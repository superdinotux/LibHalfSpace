/*
 * OKADA.cpp
 *
 *  Created on: 24/ott/2014
 *      Author: dino
 */

#include "OKADA.h"

#include <math.h>

#include "PRINT_TEMPLATE.h"



OKADA::OKADA(){
	// TODO Auto-generated constructor stub
}

OKADA::~OKADA() {
	// TODO Auto-generated destructor stub
}




void OKADA::NEW(double L_i,double W_i,double theta_deg_i,double phi_deg_i)
{

	theta_deg = theta_deg_i;
	phi_deg = phi_deg_i;

	double pi=acos(-1);

	double angolo_piatto = 180;
	double delta=(theta_deg/angolo_piatto)*pi;

	double COS_DELTA,SIN_DELTA;

	if((theta_deg == 90) || (theta_deg == -90))
	{
		SIN_DELTA=sin(delta);
		COS_DELTA=0;
	}
	else if((theta_deg == 180) || (theta_deg == -180) || (theta_deg == 0))
	{
		SIN_DELTA=0;
		COS_DELTA=cos(delta);
	}
	else
	{
		SIN_DELTA=sin(delta);
		COS_DELTA=cos(delta);
	}


	L=L_i;
	W=W_i;


	double ccp[2],pcp[2];

	double xccI=-L*0.5;
	double yccI=-W*0.5;

	double xpcI=0;
	double ypcI=0;


	double c1=-z0;

	double xc,yc,zc,xp,yp,zp;

    double n_x,n_y,n_z;
	double ts_x,ts_y,ts_z;
	double td_x,td_y,td_z;

	ccp[0]=xccI;
	ccp[1]=yccI;

	pcp[0]=xpcI;
	pcp[1]=ypcI;

	xc =  ccp[0];
	yc =  ccp[1]*COS_DELTA;
	zc =  ccp[1]*SIN_DELTA-c1;

	xp =  pcp[0];
	yp =  pcp[1]*COS_DELTA;
	zp =  pcp[1]*SIN_DELTA-c1;

	n_x =  0;
	n_y = -SIN_DELTA;
	n_z =  COS_DELTA;

    ts_x = 1;
    ts_y = 0;
    ts_z = 0;

    td_x = 0;
    td_y = COS_DELTA;
    td_z = SIN_DELTA;

	BE.NEW(xc,yc,zc,xp,yp,zp,0.,theta_deg,L,W,n_x,n_y,n_z,ts_x,ts_y,ts_z,td_x,td_y,td_z);

}



void OKADA::SET(MEDIUM MEDIUM_PAR_i,double INPUT_1,double INPUT_2,double INPUT_3,int FLAG_COMP_i)
{

	MEDIUM_PAR = MEDIUM_PAR_i;

	BE.PUT_BV(INPUT_1,INPUT_2,INPUT_3);

	FLAG_COMP = FLAG_COMP_i;

}



void OKADA::DISPLACEMENT(double x,double y,double z,double U[])
{
	// ! input:    x,y,z   coordinates of the point in which the displacement will be calculated
	//
	// ! output:   U       displacement induced by the dislocation

	double x1,y1;

	U[0]=0e0;
	U[1]=0e0;
	U[2]=0e0;

	x1=x-x0;
	y1=y-y0;

	Rotazione_coord(x1,y1,-phi_deg);

	BE.DISPLACEMENT(FLAG_COMP,MEDIUM_PAR,x1,y1,z,U);

	Rotazione_vettore(U,-phi_deg);

}




void OKADA::STRAIN(double x,double y,double z,double S[])
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

	x1=x-x0;
	y1=y-y0;

	Rotazione_coord(x1,y1,-phi_deg);

	BE.STRAIN(FLAG_COMP,MEDIUM_PAR,x1,y1,z,S);

	Rotazione_tensore(S,-phi_deg);

}




void OKADA::STRESS(double x,double y,double z,double S[])
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

	x1=x-x0;
	y1=y-y0;

	Rotazione_coord(x1,y1,-phi_deg);

	BE.STRESS(FLAG_COMP,MEDIUM_PAR,x1,y1,z,S);

	Rotazione_tensore(S,-phi_deg);

}




void OKADA::PRINT()
{

	string NAME_DATAFILE;

	string SUFFIX("_MODEL.dat");

	NAME_DATAFILE = label + SUFFIX;

	ofstream out_SOURCE_MODEL;
	out_SOURCE_MODEL.open(NAME_DATAFILE.c_str());


	double phi = acos(-1)*phi_deg/180e0;

	double cos_phi = cos(phi);
	double sin_phi = sin(phi);

	BE3D_PARAM PARAM;
	BE.GET_PARAM(PARAM);

	double x_pc = (cos_phi * PARAM.GEOM.pc[0] - sin_phi * PARAM.GEOM.pc[1]) * 1e-3;
	double y_pc = (sin_phi * PARAM.GEOM.pc[0] + cos_phi * PARAM.GEOM.pc[1]) * 1e-3;
	double z_pc = 		   	 PARAM.GEOM.pc[2]								  * 1e-3;

	double x_cc = (cos_phi * PARAM.GEOM.cc[0] - sin_phi * PARAM.GEOM.cc[1]) * 1e-3;
	double y_cc = (sin_phi * PARAM.GEOM.cc[0] + cos_phi * PARAM.GEOM.cc[1]) * 1e-3;
	double z_cc = 		   	 PARAM.GEOM.cc[2]								* 1e-3;

	out_SOURCE_MODEL << setw(12) << right << x_pc << " "
					 << setw(12) << right << y_pc << " "
					 << setw(12) << right << z_pc << " "
					 << setw(12) << right << x_cc << " "
					 << setw(12) << right << y_cc << " "
					 << setw(12) << right << z_cc << endl;

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

	PLOT_SOURCE_MODEL << "splot '" << NAME_DATAFILE  << "' using 1:2:3,"
						           << NAME_DATAFILE  << "' using 4:5:6"             	<< endl;

}



void OKADA::PRINT(CONSOLE &out)
{

	int int_flag,flag_print;

	string filename;

	out.GET_int_flag(int_flag);
	out.GET_flag_print(flag_print);

	ostream *stream;
	stream = out.out_stream;


	switch(int_flag)
	{
		case (1):
            PRINT_GEOM_PROP(stream);
		break;

		case (10):
			PRINT_MAIN_FEATURES(stream,flag_print);
		break;

		default:
			*stream << "Choice " << int_flag << " is not available" << endl;
		break;
	}


}



void OKADA::PRINT_GEOM_PROP(ostream *out_stream)
{

	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;

	*out_stream << "Position of source center at:" << endl;
	*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
	*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
	*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

	*out_stream << endl;

	*out_stream << "Okada source dimensions:" << endl;
	*out_stream << "L = " << L*1e-3 << " km" << endl;
	*out_stream << "W = " << W*1e-3 << " km" << endl;

	*out_stream << endl;

	*out_stream << "Source orientation:"   << endl;
	*out_stream << "phi = "   << phi_deg   << endl;
	*out_stream << "theta = " << theta_deg << endl;

	*out_stream << endl;

}




void OKADA::PRINT(DATAFILE &out)
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








void OKADA::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);

	double mu=MEDIUM_PAR.GET_mu();
	double nu=MEDIUM_PAR.GET_nu();

	if(flag_print == 0)
	{
		*out_stream << label << " - Main features" << endl;

		*out_stream << endl;

		*out_stream << "Position of source center at:" << endl;
		*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
		*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
		*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

		*out_stream << endl;

		*out_stream << "Geometric properties:" << endl;
		*out_stream << "L = " << L << " m" << endl;
		*out_stream << "W = " << W << " m" << endl;

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

		*out_stream << "# " << "Geometric properties:" << endl;
		*out_stream << "# " << "L = " << L << " m" << endl;
		*out_stream << "# " << "W = " << W << " m" << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;
	}

}

