/*
 * PARAL.cpp
 *
 *  Created on: 01/lug/2014
 *      Author: dino
 */

#include <cmath>
#include <iomanip>
#include <typeinfo>


using namespace std;


#include "PARAL.h"
#include "GRID.h"
#include "PRINT_TEMPLATE.h"
#include "foo.h"
#include "MEDIUM.h"


extern "C" {
void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
}


void F(int &N_a,int &N_b,int NUM,double h_a,double h_b,double passo);

void F(int &N_a,int &N_b,int NUM,double h_a,double h_b,double passo)
{
    double cfr1[NUM],cfr2[NUM];
    double min1,min2;

	for(int tk=0; tk < NUM; tk++)
	{
		double passo_a=h_a/(tk+1);
		cfr1[tk]=abs((passo_a-passo)/passo);

		double passo_b=h_b/(tk+1);
		cfr2[tk]=abs((passo_b-passo)/passo);
	}

	min1=cfr1[0];
	min2=cfr2[0];

	for(int tk=0; tk < NUM; tk++)
	{
		if(cfr1[tk] <= min1)
		{
			min1=cfr1[tk];
			N_a=tk+1;
		}

		if(cfr2[tk] <= min2)
		{
			min2=cfr2[tk];
			N_b=tk+1;
		}
	}
}



void PARAL::NEW(double h1,double h2,double h3)
{
	int NUM=NDP;

    int N_1=0;
    int N_2=0;
    int N_3=0;

    h_1=h1;
    h_2=h2;
    h_3=h3;

    double passo_1;
    double passo_2;
    double passo_3;


   	if(((h_1 == h_2) && (h_1 == h_3)))
    {
    	N_1=NUM;
		N_2=NUM;
		N_3=NUM;

		passo_1=h_1/N_1;
		passo_2=h_2/N_2;
		passo_3=h_3/N_3;
   	}
   	else if((h_1 >= h_2) && (h_1 >= h_3))
   	{
    	N_1=NUM;
		passo_1=h_1/N_1;

		F(N_2,N_3,NUM,h_2,h_3,passo_1);

		passo_2=h_2/N_2;
		passo_3=h_3/N_3;
    }
    else if((h_2 >= h_1) && (h_2 >= h_3))
    {
    	N_2=NUM;
    	passo_2=h_2/N_2;

		F(N_1,N_3,NUM,h_1,h_3,passo_2);

    	passo_1=h_1/N_1;
    	passo_3=h_3/N_3;
    }
    else if((h_3 >= h_1) && (h_3 >= h_2))
    {
    	N_3=NUM;
		passo_3=h_3/N_3;

		F(N_1,N_2,NUM,h_1,h_2,passo_3);

		passo_1=h_1/N_1;
		passo_2=h_2/N_2;
    }
    else
    {
    	cout << "ERROR" << endl;

    	exit(1);
    }


   	GRID_PAR = GRID(N_1,N_2,N_3);

   	NBE = GRID_PAR.TOTAL + 1;
   	BE = new BE3D [NBE];


    double c_1=-z0;

    int cont1,cont2;

    double passo_x=passo_1;
    double passo_y=passo_2;
    double passo_z=passo_3;

    double xI=h1*0.5;

    double yccI=-h2*0.5;
    double zccI=-c_1+h3*0.5;

    double ypcI=-(N_2-1)*passo_y*0.5;
    double zpcI=-c_1+(N_3-1)*passo_z/2;

    double zcc,zpc;
    double ycc,ypc;

    for(int tz=0; tz < N_3; tz++)
    {
        zcc=zccI-passo_z*tz;
        zpc=zpcI-passo_z*tz;

        for(int ty=0; ty < N_2; ty++)
        {
            ycc=yccI+passo_y*ty;
            ypc=ypcI+passo_y*ty;

            cont1=tz*N_2+ty;
            cont2=cont1+N_2*N_3;

            BE[cont1].NEW( xI,ycc,zcc, xI,ypc,zpc,-90.,-90.,h2/N_2,h3/N_3,1,0,0,0,1,0,0,0,1);
            BE[cont2].NEW(-xI,ycc,zcc,-xI,ypc,zpc,-90.,-90.,h2/N_2,h3/N_3,1,0,0,0,1,0,0,0,1);
        }
    }


    double yI=h2*0.5;

    double xccI=-h1*0.5;
    zccI=-c_1+h3*0.5;

    double xpcI=-(N_1-1)*passo_x*0.5;
    zpcI=-c_1+(N_3-1)*passo_z/2;

    double xcc,xpc;

    for(int tz=0; tz < N_3; tz++)
    {
        zcc=zccI-passo_z*tz;
        zpc=zpcI-passo_z*tz;

        for(int tx=0; tx < N_1; tx++)
        {
            xcc=xccI+passo_x*tx;
            xpc=xpcI+passo_x*tx;

            cont1=N_2*N_3*2+tz*N_1+tx;
            cont2=cont1+N_1*N_3;

            BE[cont1].NEW(xcc, yI,zcc,xpc, yI,zpc,0.,-90.,h1/N_1,h3/N_3,0,1,0,1,0,0,0,0,1);
            BE[cont2].NEW(xcc,-yI,zcc,xpc,-yI,zpc,0.,-90.,h1/N_1,h3/N_3,0,1,0,1,0,0,0,0,1);
        }
    }


    double zI=-c_1+h3*0.5;

    xccI=-h1*0.5;
    yccI=-h2*0.5;

    xpcI=-(N_1-1)*passo_x*0.5;
    ypcI=-(N_2-1)*passo_y*0.5;

    for(int tx=0; tx < N_1; tx++)
    {
        xcc=xccI+passo_x*tx;
        xpc=xpcI+passo_x*tx;

        for(int ty=0; ty < N_2; ty++)
        {
            ycc=yccI+passo_y*ty;
            ypc=ypcI+passo_y*ty;

            cont1=(N_2+N_1)*N_3*2+tx*N_2+ty;
            cont2=cont1+N_1*N_2;

            BE[cont1].NEW(xcc,ycc,zI   ,xpc,ypc,zI   ,0.,0.,h1/N_1,h2/N_2,0,0,1,1,0,0,0,1,0);
            BE[cont2].NEW(xcc,ycc,zI-h3,xpc,ypc,zI-h3,0.,0.,h1/N_1,h2/N_2,0,0,1,1,0,0,0,1,0);
        }
    }


	double factor1=0.9999;

	xcc = -(h1/N_1)*0.5*factor1;
	ycc = -(h2/N_2)*0.5*factor1;

	zI  = -c_1+1e-1;

	xpc = 0;
	ypc = 0;

    BE[NBE-1].NEW(xcc,ycc,zI,xpc,ypc,zI,0.,0.,h1*factor1,h2*factor1,0,0,1,1,0,0,0,1,0);
}




void PARAL::SET(MEDIUM MEDIUM_PAR_i,int option_flag,double INPUT_1,double INPUT_2,double INPUT_3,int FLAG_COMP_i)
{
	MEDIUM_PAR = MEDIUM_PAR_i;

	switch (option_flag)
	{
		case (1):

				for(int i=0; i < NBE; i++)
				{
					if(i == NBE-1)
					{
						BE[i].PUT_BV(0,0,0);
					}
					else
					{
						BE[i].PUT_BV(INPUT_1,INPUT_2,INPUT_3);
					}
				}

				break;

		case (2):

				for(int i=0; i < NBE; i++)
				{
					if(i == NBE-1)
					{
						BE[i].PUT_BC(0,0,0);
					}
					else
					{
						BE[i].PUT_BC(INPUT_1,INPUT_2,INPUT_3);
					}
				}

				break;

		default:

			cout << "Option " << option_flag << " is not available" << endl;
	}

	FLAG_COMP = FLAG_COMP_i;
};




void PARAL::SOLVE()
{
	int dim=3;

	int TOTALE=NBE;
	int TOT=dim*NBE;

	double *A;
	A=new double [TOT*TOT];

	double *B;
    B=new double [TOT];

	double xn,yn,z;
	double n[3],ts[3],td[3];

	double traction_S[3],traction_D[3],traction_T[3];

	BE3D_PARAM PARAM;

	int start_ti = TOTALE-1;

	int I1,I2;

	int cont=0;
	int np=pow(NBE,2);
	int c[10]={0,0,0,0,0,0,0,0,0,0};



	cout << endl << "SOLVE method" << endl;
	cout << endl << "Computation of the influence coefficients matrix" << endl;

	for(int ti=0; ti < TOTALE; ti++)
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


		for(int tj=0; tj < TOTALE; tj++)
		{
			I2=tj*dim;

			cont = ti * NBE + tj + 1;

			BE[tj].PUT_BV(1,1,1);

			if(ti < start_ti)
			{
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
			else
			{
				double S_Displacement[3],D_Displacement[3],T_Displacement[3];
				BE[tj].DISPLACEMENT(FLAG_COMP,MEDIUM_PAR,xn,yn,z,S_Displacement,D_Displacement,T_Displacement);

				A[I2*TOT+I1]		= S_Displacement[0];
				A[I2*TOT+I1+1]		= S_Displacement[1];
				A[I2*TOT+I1+2]		= S_Displacement[2];

				A[(I2+1)*TOT+I1]	= D_Displacement[0];
				A[(I2+1)*TOT+I1+1]	= D_Displacement[1];
				A[(I2+1)*TOT+I1+2]	= D_Displacement[2];

				A[(I2+2)*TOT+I1]	= T_Displacement[0];
				A[(I2+2)*TOT+I1+1]	= T_Displacement[1];
				A[(I2+2)*TOT+I1+2]	= T_Displacement[2];
			}
		}

		Cycle_counter(cont,np,c);

	}


	int i1,i2,i3;

	double BC1,BC2,BC3;

	for(int i=0; i < TOTALE; i++)
	{
		BE[i].GET_BC(BC1,BC2,BC3);

		i1 = dim*i;
		i2 = i1 + 1;
		i3 = i1 + 2;

		B[i1] = -BC1;
		B[i2] = -BC2;
		B[i3] = -BC3;
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

	for(int i=0; i < TOTALE; i++)
	{
		i1 = dim*i;
		i2 = i1 + 1;
		i3 = i1 + 2;

		BE[i].PUT_BV(B[i1],B[i2],B[i3]);
	}

	delete A;
	delete B;

}




void PARAL::PRINT(CONSOLE &out)
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



void PARAL_f_GEOM_PROP(ostream *out_stream)
{
	*out_stream << endl;
	*out_stream << setw(7)  << "i="            ;
	*out_stream << setw(15) << "x0 (m)"        ;
	*out_stream << setw(14) << "y0 (m)"        ;
	*out_stream << setw(14) << "z0 (m)"        ;
	*out_stream << setw(18) << "angle (degree)";
	*out_stream << setw(6)  << "j=";
}



void PARAL::PRINT_GEOM_PROP(ostream *out_stream)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);

	*out_stream << label << " - Geometrical properties" << endl;

	BE3D_PARAM PARAM;

	int j=1;
	int fc=1;

	for(int i=0; i < NBE; i++)
	{

		if(j == 1)
		{
			*out_stream << endl;

			if(i == 0)
			{
				*out_stream << "FACE X+";
				PARAL_f_GEOM_PROP(out_stream);

				fc=1;
			}
			else if(i == GRID_PAR.NX)
			{
				*out_stream << "FACE X-";
				PARAL_f_GEOM_PROP(out_stream);

				fc=1;
			}
			else if(i == 2*GRID_PAR.NX)
			{
				*out_stream << "FACE Y+";
				PARAL_f_GEOM_PROP(out_stream);

				fc=2;
			}
			else if(i == 2*GRID_PAR.NX+GRID_PAR.NY)
			{
				*out_stream << "FACE Y-";
				PARAL_f_GEOM_PROP(out_stream);

				fc=2;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY))
			{
				*out_stream << "FACE Z+";
				PARAL_f_GEOM_PROP(out_stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY)+GRID_PAR.NZ)
			{
				*out_stream << "FACE Z-";
				PARAL_f_GEOM_PROP(out_stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY+GRID_PAR.NZ))
			{
				*out_stream << "IE";
				PARAL_f_GEOM_PROP(out_stream);

				fc=3;
			}

			*out_stream << endl;

		}


		BE[i].GET_PARAM(PARAM);

		double x1=PARAM.GEOM.pc[0];
		double y1=PARAM.GEOM.pc[1];

		double fi=(acos(-1)/180e0)*phi_degree;

		double cos_fi=cos(fi);
		double sin_fi=sin(fi);

		double x2=x1*cos_fi-y1*sin_fi;
		double y2=x1*sin_fi+y1*cos_fi;

		*out_stream << setw(8) <<   i+1;
//		*out_stream << setw(14) <<  PARAM.GEOM.pc[0];
//		*out_stream << setw(14) <<  PARAM.GEOM.pc[1];
		*out_stream << setw(14) <<  x2+x0;
		*out_stream << setw(14) <<  y2+y0;
		*out_stream << setw(14) <<  PARAM.GEOM.pc[2];
		*out_stream << setw(15) <<  PARAM.GEOM.delta_degree;
		*out_stream << setw(10) <<  j;

		*out_stream << endl;

		j++;

		if(((j == GRID_PAR.NX+1) && (fc == 1)) || ((j == GRID_PAR.NY+1) && (fc == 2)) || ((j == GRID_PAR.NZ+1) && (fc == 3)))
		{
			j=1;
		}

	}

	*out_stream << endl << endl;

}



void f_STRESS_PROP(ostream *stream)
{
	*stream << endl;
	*stream << setw(7)  << "i=";
	*stream << setw(15) << "BC S (MPa)";
	*stream << setw(15) << "BC D (MPa)";
	*stream << setw(15) << "BC T (MPa)";
	*stream << setw(7)  << "j=";
}



void f_DISL_PROP(ostream *stream)
{
	*stream << endl;
	*stream << setw(7)  << "i=";
	*stream << setw(15) << "BC S (m)";
	*stream << setw(15) << "BC D (m)";
	*stream << setw(15) << "BC T (m)";
	*stream << setw(7)  << "j=";
}


void PARAL::PRINT_STRESS_PROP(ostream *stream)
{

	double BC1,BC2,BC3;

	stream->fill(' ');
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	*stream << label << " - Stress properties"<< endl;

	int j=1;
	int fc=1;

	for(int i=0; i < NBE; i++)
	{

		if(j == 1)
		{
			*stream << endl;

			if(i == 0)
			{
				*stream << "FACE X+";
				f_STRESS_PROP(stream);

				fc=1;
			}
			else if(i == GRID_PAR.NX)
			{
				*stream << "FACE X-";
				f_STRESS_PROP(stream);

				fc=1;
			}
			else if(i == 2*GRID_PAR.NX)
			{
				*stream << "FACE Y+";
				f_STRESS_PROP(stream);

				fc=2;
			}
			else if(i == 2*GRID_PAR.NX+GRID_PAR.NY)
			{
				*stream << "FACE Y-";
				f_STRESS_PROP(stream);

				fc=2;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY))
			{
				*stream << "FACE Z+";
				f_STRESS_PROP(stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY)+GRID_PAR.NZ)
			{
				*stream << "FACE Z-";
				f_STRESS_PROP(stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY+GRID_PAR.NZ))
			{
				*stream << "IE";
				f_DISL_PROP(stream);

				fc=3;
			}

			*stream << endl;

		}

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

		*stream << setw(10) <<  j;

		*stream << endl;

		j++;

		if(((j == GRID_PAR.NX+1) && (fc == 1)) || ((j == GRID_PAR.NY+1) && (fc == 2)) || ((j == GRID_PAR.NZ+1) && (fc == 3)))
		{
			j=1;
		}

	}

	*stream << endl << endl;

}



void f_BV_PROP(ostream *stream)
{
	*stream << endl;

	*stream << setw(7)  << "i="     ;
	*stream << setw(17) << "S (m)" ;
	*stream << setw(20) << "D (m)" ;
	*stream << setw(18) << "T (m)" ;
	*stream << setw(9)  << "j=";
}



void PARAL::PRINT_BV_COMP(ostream *stream)
{

	double B_1,B_2,B_3;

	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	GET_label(label);

	*stream << label << " - Burgers vector components" << endl;

	int j=1;
	int fc=1;

	for(int i=0; i < NBE; i++)
	{

		if(j == 1)
		{
			*stream << endl;

			if(i == 0)
			{
				*stream << "FACE X+";
				f_BV_PROP(stream);

				fc=1;
			}
			else if(i == GRID_PAR.NX)
			{
				*stream << "FACE X-";
				f_BV_PROP(stream);

				fc=1;
			}
			else if(i == 2*GRID_PAR.NX)
			{
				*stream << "FACE Y+";
				f_BV_PROP(stream);

				fc=2;
			}
			else if(i == 2*GRID_PAR.NX+GRID_PAR.NY)
			{
				*stream << "FACE Y-";
				f_BV_PROP(stream);

				fc=2;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY))
			{
				*stream << "FACE Z+";
				f_BV_PROP(stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY)+GRID_PAR.NZ)
			{
				*stream << "FACE Z-";
				f_BV_PROP(stream);

				fc=3;
			}
			else if(i == 2*(GRID_PAR.NX+GRID_PAR.NY+GRID_PAR.NZ))
			{
				*stream << "IE";
				f_BV_PROP(stream);

				fc=3;
			}

			*stream << endl;

		}

		BE[i].GET_BV(B_1,B_2,B_3);

		stream->fill(' ');
		*stream << right << setw(8)  << i+1;

		stream->fill(' ');
		*stream << setw(7) << " ";

		stream->fill(' ');
		*stream << right << showpos << setw(8) << B_1 ;

		stream->fill(' ');
		*stream << setw(7) << " ";

		stream->fill(' ');
		*stream << right << setw(8) << B_2 ;

		stream->fill(' ');
		*stream << setw(7) << " ";

		stream->fill(' ');
		*stream << right << noshowpos << setw(8) << B_3 ;

		*stream << setw(10) <<  j;

		*stream << endl;

		j++;

		if(((j == GRID_PAR.NX+1) && (fc == 1)) || ((j == GRID_PAR.NY+1) && (fc == 2)) || ((j == GRID_PAR.NZ+1) && (fc == 3)))
		{
			j=1;
		}

	}

	*stream << endl << endl;

}




void PARAL::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
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
		*out_stream << "h_1 = " << h_1 << " m" << endl;
		*out_stream << "h_2 = " << h_2 << " m" << endl;
		*out_stream << "h_3 = " << h_3 << " m" << endl;

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

		*out_stream << "# " << "Geometric properties:" << endl;
		*out_stream << "# " << "h_1 = " << h_1 << " m" << endl;
		*out_stream << "# " << "h_2 = " << h_2 << " m" << endl;
		*out_stream << "# " << "h_3 = " << h_3 << " m" << endl;

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



void PARAL::PRINT(DATAFILE &out)
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





