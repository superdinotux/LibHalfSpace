/*
 * SOURCE_BE3D.cpp
 *
 *  Created on: 02/ott/2014
 *      Author: dino
 */


#include "SOURCE_BE3D.h"

#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>

#include "PRINT_TEMPLATE.h"


SOURCE_BE3D::SOURCE_BE3D(){
	// TODO Auto-generated constructor stub
}

SOURCE_BE3D::~SOURCE_BE3D() {
	// TODO Auto-generated destructor stub
}



extern "C" {
void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
}



void SOURCE_BE3D::DISPLACEMENT(double x,double y,double z,double U[])
{
// ! input:    x,y,z   coordinates of the point in which the displacement will be calculated
//
// ! output:   U       displacement induced by the dislocation

	double x1,y1;

	U[0]=0e0;
	U[1]=0e0;
	U[2]=0e0;

	for(int i=0; i < NBE; i++)
	{
		double U_i[3]={0,0,0};

		x1=x-x0;
		y1=y-y0;

	    Rotazione_coord(x1,y1,-phi_degree);

		BE[i].DISPLACEMENT(FLAG_COMP,MEDIUM_PAR,x1,y1,z,U_i);

		U[0] += U_i[0];
		U[1] += U_i[1];
		U[2] += U_i[2];
	}

	Rotazione_vettore(U,-phi_degree);

}




void SOURCE_BE3D::STRAIN(double x,double y,double z,double S[])
{
// ! input:    x,y,z   coordinates of the point in which the displacement will be calculated
//
// ! output:   S 	   components of the stress field induced by the dislocation

	double x1,y1;

	S[0]=0e0;
	S[1]=0e0;
	S[2]=0e0;
	S[3]=0e0;
	S[4]=0e0;
	S[5]=0e0;

	for(int i=0; i < NBE; i++)
	{
		double S_i[6]={0,0,0,0,0,0};

		x1=x-x0;
		y1=y-y0;

	    Rotazione_coord(x1,y1,-phi_degree);

		BE[i].STRAIN(FLAG_COMP,MEDIUM_PAR,x1,y1,z,S_i);

		S[0] += S_i[0];
		S[1] += S_i[1];
		S[2] += S_i[2];
		S[3] += S_i[3];
		S[4] += S_i[4];
		S[5] += S_i[5];
	}

	Rotazione_tensore(S,-phi_degree);

}




void SOURCE_BE3D::STRESS(double x,double y,double z,double S[])
{
// ! input:    x,y,z   coordinates of the point in which the displacement will be calculated
//
// ! output:   S 	   components of the stress field induced by the dislocation

	double x1,y1;

	S[0]=0e0;
	S[1]=0e0;
	S[2]=0e0;
	S[3]=0e0;
	S[4]=0e0;
	S[5]=0e0;

	for(int i=0; i < NBE; i++)
	{
		double S_i[6]={0,0,0,0,0,0};

		x1=x-x0;
		y1=y-y0;

	    Rotazione_coord(x1,y1,-phi_degree);

		BE[i].STRESS(FLAG_COMP,MEDIUM_PAR,x1,y1,z,S_i);

		S[0] += S_i[0];
		S[1] += S_i[1];
		S[2] += S_i[2];
		S[3] += S_i[3];
		S[4] += S_i[4];
		S[5] += S_i[5];
	}

	Rotazione_tensore(S,-phi_degree);

}




void SOURCE_BE3D::SOLVE()
{
	int dim=3;

	int TOT=dim*NBE;

	double *A;
	A=new double [TOT*TOT];

	double *B;
    B=new double [TOT];

	double xn,yn,z;
	double n[3],ts[3],td[3];

	double traction_S[3],traction_D[3],traction_T[3];

	BE3D_PARAM PARAM;

	int I1,I2;

	int cont=0;
	int np=pow(NBE,2);
	int c[10]={0,0,0,0,0,0,0,0,0,0};


	cout << endl << "SOLVE method" << endl;
	cout << endl << "Computation of the influence coefficients matrix" << endl;

	for(int ti=0; ti < NBE; ti++)
	{

		I1=ti*dim;

		BE[ti].GET_PARAM(PARAM);

		xn = PARAM.GEOM.pc[0];
		yn = PARAM.GEOM.pc[1];
		z  = PARAM.GEOM.pc[2];

		n[0] = PARAM.GEOM.n[0];
		n[1] = PARAM.GEOM.n[1];
		n[2] = PARAM.GEOM.n[2];

		ts[0] = PARAM.GEOM.ts[0];
		ts[1] = PARAM.GEOM.ts[1];
		ts[2] = PARAM.GEOM.ts[2];

		td[0] = PARAM.GEOM.td[0];
		td[1] = PARAM.GEOM.td[1];
		td[2] = PARAM.GEOM.td[2];


		for(int tj=0; tj < NBE; tj++)
		{

			I2=tj*dim;

			cont = ti * NBE + tj + 1;

			BE[tj].PUT_BV(1,1,1);

			double S_Stress[6],D_Stress[6],T_Stress[6];
			BE[tj].STRESS(FLAG_COMP,MEDIUM_PAR,xn,yn,z,S_Stress,D_Stress,T_Stress);


			traction_S[0] = S_Stress[0]*n[0]+S_Stress[3]*n[1]+S_Stress[4]*n[2];
			traction_S[1] = S_Stress[3]*n[0]+S_Stress[1]*n[1]+S_Stress[5]*n[2];
			traction_S[2] = S_Stress[4]*n[0]+S_Stress[5]*n[1]+S_Stress[2]*n[2];

			traction_D[0] = D_Stress[0]*n[0]+D_Stress[3]*n[1]+D_Stress[4]*n[2];
			traction_D[1] = D_Stress[3]*n[0]+D_Stress[1]*n[1]+D_Stress[5]*n[2];
			traction_D[2] = D_Stress[4]*n[0]+D_Stress[5]*n[1]+D_Stress[2]*n[2];

			traction_T[0] = T_Stress[0]*n[0]+T_Stress[3]*n[1]+T_Stress[4]*n[2];
			traction_T[1] = T_Stress[3]*n[0]+T_Stress[1]*n[1]+T_Stress[5]*n[2];
			traction_T[2] = T_Stress[4]*n[0]+T_Stress[5]*n[1]+T_Stress[2]*n[2];


			A[I2*TOT+I1]   	   = traction_S[0]*ts[0] + traction_S[1]*ts[1] + traction_S[2]*ts[2];
			A[I2*TOT+I1+1] 	   = traction_S[0]*td[0] + traction_S[1]*td[1] + traction_S[2]*td[2];
			A[I2*TOT+I1+2] 	   = traction_S[0]*n[0]  + traction_S[1]*n[1]  + traction_S[2]*n[2];

			A[(I2+1)*TOT+I1]   = traction_D[0]*ts[0] + traction_D[1]*ts[1] + traction_D[2]*ts[2];
			A[(I2+1)*TOT+I1+1] = traction_D[0]*td[0] + traction_D[1]*td[1] + traction_D[2]*td[2];
			A[(I2+1)*TOT+I1+2] = traction_D[0]*n[0]  + traction_D[1]*n[1]  + traction_D[2]*n[2];

			A[(I2+2)*TOT+I1]   = traction_T[0]*ts[0] + traction_T[1]*ts[1] + traction_T[2]*ts[2];
			A[(I2+2)*TOT+I1+1] = traction_T[0]*td[0] + traction_T[1]*td[1] + traction_T[2]*td[2];
			A[(I2+2)*TOT+I1+2] = traction_T[0]*n[0]  + traction_T[1]*n[1]  + traction_T[2]*n[2];


		}

		Cycle_counter(cont,np,c);

	}


	double BC_1,BC_2,BC_3;

	int i1,i2,i3;

	for(int i=0; i < NBE; i++)
	{
		BE[i].GET_BC(BC_1,BC_2,BC_3);

		i1 = dim*i;
		i2 = i1 + 1;
		i3 = i1 + 2;

		B[i1] = -BC_1;
		B[i2] = -BC_2;
		B[i3] = -BC_3;
	}


	int ncB=1;
	int IPIV[TOT];
	int INFO;

	cout << endl << "Try to find a solution (using DGESV routine)" << endl;

	dgesv_(&TOT, &ncB, A, &TOT, IPIV, B, &TOT, &INFO);

	if(INFO == 0)
	{
		cout << endl << "Solution identified (INFO = " << INFO << ")" << endl;
	}
	else
	{
		cout << endl << "INFO ="  << INFO << endl;

		if(INFO < 0)
		{
			cout << "If INFO = -i, the i-th argument had an illegal value." << endl;
		}

		if(INFO > 0)
		{
			cout << "The factor U(i,i) is exactly singular (A = P * L * U) and so the solution could not be computed." << endl;
			exit(0);
		}

	}

	cout << endl << endl;

	for(int i=0; i < NBE; i++)
	{
		i1 = dim*i;
		i2 = i1 + 1;
		i3 = i1 + 2;

		BE[i].PUT_BV(B[i1],B[i2],B[i3]);
	}

	delete A;
	delete B;

}




void SOURCE_BE3D::PRINT()
{

	string NAME_DATAFILE;

	string SUFFIX("_MODEL.dat");

	NAME_DATAFILE = label + SUFFIX;

	ofstream out_SOURCE_MODEL;
	out_SOURCE_MODEL.open(NAME_DATAFILE.c_str());


	double x,y,z;

	double phi = acos(-1)*phi_degree/180e0;

	double cos_phi = cos(phi);
	double sin_phi = sin(phi);

	BE3D_PARAM PARAM;

	for(int i=0; i < NBE; i++)
	{
		BE[i].GET_PARAM(PARAM);

		x = (  cos_phi * PARAM.GEOM.pc[0] + sin_phi * PARAM.GEOM.pc[1])*1e-3;
		y = (- sin_phi * PARAM.GEOM.pc[0] + cos_phi * PARAM.GEOM.pc[1])*1e-3;
		z = PARAM.GEOM.pc[2]*1e-3;

		out_SOURCE_MODEL << setw(12) << right << x << " "
						 << setw(12) << right << y << " "
						 << setw(12) << right << z << endl;
	}



	string PLOT_NAMEFILE;

	SUFFIX = "_MODEL.gnuplot";

	PLOT_NAMEFILE = label + SUFFIX;


	string IMAGENAME;

	SUFFIX = "_MODEL.eps";

	IMAGENAME = label + SUFFIX;


	ofstream PLOT_SOURCE_MODEL;
	PLOT_SOURCE_MODEL.open(PLOT_NAMEFILE.c_str());

	PLOT_SOURCE_MODEL << "# Gnuplot script											"	<< endl;
	PLOT_SOURCE_MODEL << "  														"	<< endl;
	PLOT_SOURCE_MODEL << "x0=10														"	<< endl;
	PLOT_SOURCE_MODEL << "y0=10														"	<< endl;
	PLOT_SOURCE_MODEL << "z0=10														"	<< endl;
	PLOT_SOURCE_MODEL << "  														"	<< endl;
	PLOT_SOURCE_MODEL << "set terminal postscript eps enhanced color                "	<< endl;
	PLOT_SOURCE_MODEL << "set nokey													"	<< endl;
	PLOT_SOURCE_MODEL << "set title '" << label << "'								"	<< endl;

	PLOT_SOURCE_MODEL << "set xtics -20,5,20 offset 25,9 							"	<< endl;
	PLOT_SOURCE_MODEL << "set ytics -20,5,20 offset -42.5,5 						"	<< endl;
	PLOT_SOURCE_MODEL << "set ztics 1 offset 0.25,-0.7	 							"	<< endl;

	PLOT_SOURCE_MODEL << "set xlabel 'x (km)' offset 35,12 							"	<< endl;
	PLOT_SOURCE_MODEL << "set ylabel 'y (km)' offset -56,7 							"	<< endl;
	PLOT_SOURCE_MODEL << "set zlabel 'z (km)' 										"	<< endl;

	PLOT_SOURCE_MODEL << "set output '" << IMAGENAME << "' 							"	<< endl;
	PLOT_SOURCE_MODEL << "set xyplane at 0      									"	<< endl;
	PLOT_SOURCE_MODEL << "set grid noxtics front									"	<< endl;
	PLOT_SOURCE_MODEL << "set grid xtics ytics  									"	<< endl;
	PLOT_SOURCE_MODEL << "set border 4095 lw 2 back									"	<< endl;
	PLOT_SOURCE_MODEL << "set pointsize 0.5											"	<< endl;

	PLOT_SOURCE_MODEL << "set xrange [-x0:+x0]										"	<< endl;
	PLOT_SOURCE_MODEL << "set yrange [-y0:+y0]										"	<< endl;
	PLOT_SOURCE_MODEL << "set zrange [-z0:0]										"	<< endl;

	PLOT_SOURCE_MODEL << "splot '" << NAME_DATAFILE  << "' 				  			" 	<< endl;
}




void SOURCE_BE3D::PRINT(CONSOLE &out)
{

	int int_flag,flag_print;

	string filename;

	out.GET_int_flag(int_flag);
	out.GET_flag_print(flag_print);

	ostream *stream;
	stream = out.out_stream;

	switch(int_flag)
	{
		case (0):
			PRINT_GEOM_PROP(stream);
			PRINT_STRESS_PROP(stream);
			PRINT_BV_COMP(stream);
		break;

		case (1):
			PRINT_GEOM_PROP(stream);
		break;

		case (2):
			PRINT_STRESS_PROP(stream);
		break;

		case (3):
    	    PRINT_BV_COMP(stream);
		break;

		case (10):
			PRINT_MAIN_FEATURES(stream,flag_print);
		break;

		default:
			*stream << "Not implemented for this source type" << endl;
	}

}




void SOURCE_BE3D::PRINT_GEOM_PROP(ostream *out_stream)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);


	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;
	*out_stream << setw(7)  << "i="            ;
	*out_stream << setw(15) << "x0 (m)"        ;
	*out_stream << setw(14) << "y0 (m)"        ;
	*out_stream << setw(14) << "z0 (m)"        ;
	*out_stream << setw(18) << "angle (degree)";
	*out_stream << endl;


	double fi=(acos(-1)/180e0)*phi_degree;

	double cos_fi=cos(fi);
	double sin_fi=sin(fi);

	BE3D_PARAM PARAM;

	for(int i=0; i < NBE; i++)
	{

		BE[i].GET_PARAM(PARAM);

		double x1=PARAM.GEOM.pc[0];
		double y1=PARAM.GEOM.pc[1];

		double x2=x1*cos_fi-y1*sin_fi;
		double y2=x1*sin_fi+y1*cos_fi;

		*out_stream << setw(7)  <<  i;
		*out_stream << setw(14) <<  x2+x0;
		*out_stream << setw(14) <<  y2+y0;
		*out_stream << setw(14) <<  PARAM.GEOM.pc[2];
		*out_stream << setw(15) <<  PARAM.GEOM.delta_degree;
		*out_stream << endl;

	}

	*out_stream << endl << endl;

}




void SOURCE_BE3D::PRINT_STRESS_PROP(ostream *stream)
{

	double BC1,BC2,BC3;

	stream->fill(' ');
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	*stream << label << " - Stress properties"<< endl;

	*stream << endl;
	*stream << setw(7)  << "i="        ;
	*stream << setw(17) << "BC S (MPa)" ;
	*stream << setw(14) << "BC D (MPa)" ;
	*stream << setw(14) << "BC T (MPa)" ;
	*stream << endl;


	for(int i=0; i < NBE; i++)
	{

		BE[i].GET_BC(BC1,BC2,BC3);

		stream->fill(' ');
		*stream << right << setw(8)  << i+1;

		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << BC1*1e-6;

		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << BC2*1e-6;

		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << BC3*1e-6;

		*stream << endl;

	}

	*stream << endl << endl;

}




void SOURCE_BE3D::PRINT_BV_COMP(ostream *stream)
{

	double B_1,B_2,B_3;

//	stream->fill(' ');
//	stream->setf(ios::right);
//	stream->setf(ios::showpoint);
//	stream->precision(5);

	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	*stream << label << " - Burgers vector components" << endl;

	*stream << endl;
	*stream << setw(7)  << "i="        ;
	*stream << setw(15) << "S (m)" ;
	*stream << setw(19) << "D (m)" ;
	*stream << setw(19) << "T (m)" ;
	*stream << endl;


	for(int i=0; i < NBE; i++)
	{

		BE[i].GET_BV(B_1,B_2,B_3);

		stream->fill(' ');
		*stream << right << noshowpos << setw(8)  << i+1;

		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << showpos << setw(8) << B_1 ;

		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << B_2 ;
		stream->fill(' ');
		*stream << setw(6) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << B_3 ;

		*stream << endl;

	}

	*stream << endl << endl;

}




void SOURCE_BE3D::PRINT(DATAFILE &out)
{
	int flag_datatype;
	out.GET_flag_datatype(flag_datatype);

	if(flag_datatype == 1)
	{
		int int_flag;
		out.GET_int_flag(int_flag);

		string filename;

		int N_1,N_2;
		double DX1,DX2;

		out.GET_file_parameters(filename,N_1,N_2,DX1,DX2);

		switch(int_flag)
		{
			case (0):
        		print_maps(*this,filename,N_1,N_2,DX1,DX2);
				break;

			case (1):
            	print_map_displ(*this,filename,N_1,N_2,DX1,DX2);
				break;

			case (2):
            	print_map_strain(*this,filename,N_1,N_2,DX1,DX2);
				break;

			case (3):
            	print_map_stress(*this,filename,N_1,N_2,DX1,DX2);
				break;

			default:
				cout << "Opzione " << int_flag << " non disponibile";
				break;
		}

	}
	else
	{
		cout << label << ": datafile object is not compatible with print method" << endl;
	}


}







void SOURCE_BE3D::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);

	double mu=MEDIUM_PAR.GET_mu();
	double nu=MEDIUM_PAR.GET_nu();

	double GET_nu();

	if(flag_print == 0)
	{

		*out_stream << label << " - Main features" << endl;

		*out_stream << endl;

		*out_stream << "Position of source center at:" << endl;
		*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
		*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
		*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

		*out_stream << endl;

		*out_stream << "Parameters of the boundary element grid:" << endl;
		*out_stream << "NDP = " << NDP << endl;
		*out_stream << "NBE = " << NBE << endl;

		*out_stream << endl;

		*out_stream << "The elastic parameter of the medium are:" << endl;
		*out_stream << "mu = " << mu << endl;
		*out_stream << "nu = " << nu << endl;

		*out_stream << endl;

	}
	else
	{

		*out_stream << "# " << label << " - Main features" << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "Position of source center at:" << endl;
		*out_stream << "# " << "x0 = " << x0*1e-3 << " km"     << endl;
		*out_stream << "# " << "y0 = " << y0*1e-3 << " km"     << endl;
		*out_stream << "# " << "z0 = " << z0*1e-3 << " km"     << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "Parameters of the boundary element grid:" << endl;
		*out_stream << "# " << "NDP = " << NDP << endl;
		*out_stream << "# " << "NBE = " << NBE << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;

	}

}




void SOURCE_BE3D::PUT_x0_y0_z0(double x0_i, double y0_i, double z0_i)
{

	double Dx,Dy,Dz;

	Dx = x0_i - x0;
	Dy = y0_i - y0;
	Dz = z0_i - z0;

	x0 = x0_i;
	y0 = y0_i;
	z0 = z0_i;

	for(int ti=0 ; ti < NBE; ti++)
	{
		BE[ti].SHIFT(Dx, Dy, Dz);
	}

}







