/*
 * SILL.cpp
 *
 *  Created on: 07/nov/2014
 *      Author: dino
 */

#include "SILL.h"

#include <iomanip>

#include "PRINT_TEMPLATE.h"


SILL::SILL() {
	// TODO Auto-generated constructor stub

}

SILL::~SILL() {
	// TODO Auto-generated destructor stub
}



void SILL::NEW(double a_i,double b_i,double delta_gradi)
{

	theta_degree = delta_gradi;

	double pi=acos(-1);

	double angolo_piatto = 180;
	double delta=(delta_gradi/angolo_piatto)*pi;

	double COS_DELTA,SIN_DELTA;

	if((delta_gradi == 90) || (delta_gradi == -90))
	{
		SIN_DELTA=sin(delta);
		COS_DELTA=0;
	}
	else if((delta_gradi == 180) || (delta_gradi == -180) || (delta_gradi == 0))
	{
		SIN_DELTA=0;
		COS_DELTA=cos(delta);
	}
	else
	{
		SIN_DELTA=sin(delta);
		COS_DELTA=cos(delta);
	}


	int NUM=NDP;

	int N_1=0;
	int N_2=0;

	a=a_i;
	b=b_i;

	double LATO_h_1=a*2;
	double LATO_h_2=b*2;

	double passo_1;
	double passo_2;

	double cfr[NUM];
	double min;


	if(a == b)
	{
		N_1=NUM;
		N_2=NUM;

		passo_1=LATO_h_1/N_1;
		passo_2=LATO_h_2/N_2;
	}
	else if(a > b)
	{
		N_1=NUM;

		passo_1=LATO_h_1/N_1;

		for(int tk=0; tk < N_1; tk++)
		{
			passo_2=LATO_h_2/(tk+1);
			cfr[tk]=abs((passo_2-passo_1)/passo_1);
		}

		min=cfr[0];

		for(int tk=0; tk < N_1; tk++)
		{

			if(cfr[tk] <= min)
			{
				min=cfr[tk];
				N_2=tk+1;
			}

		}

		passo_2=LATO_h_2/N_2;
	}
	else
	{
		N_2=NUM;
		passo_2=LATO_h_2/N_2;

		for(int tk=0; tk < N_2; tk++)
		{
			passo_1=LATO_h_1/(tk+1);
			cfr[tk]=abs((passo_1-passo_2)/passo_2);
		}


		min=cfr[0];

		for(int tk=0; tk < N_2; tk++)
		{
			if(cfr[tk] <= min)
			{
				min=cfr[tk];
				N_1=tk+1;
			}
		}

		passo_1=LATO_h_1/N_1;
	}


	double x_est,x_cen,y_est,y_cen;

	double vector_elemets_x_cen[N_1];
	double vector_elemets_x_est[N_1];
	double vector_elemets_y_cen[N_2];
	double vector_elemets_y_est[N_2];


	double x_est_0 = -a;
	double y_est_0 =  b - passo_2;

	for (int tk=0; tk < N_1; tk++)
	{
		x_est=x_est_0+tk*passo_1;
		x_cen=x_est+passo_1*0.5;

		vector_elemets_x_est[tk]=x_est;
		vector_elemets_x_cen[tk]=x_cen;
	}


	for (int tk=0; tk < N_2; tk++)
	{
		y_est=y_est_0-tk*passo_2;
		y_cen=y_est+passo_2*0.5;

		vector_elemets_y_est[tk]=y_est;
		vector_elemets_y_cen[tk]=y_cen;
	}



	int tj=0;

	double x_element,y_element;
	double dist_quad;

	for (int tl=0; tl < N_2; tl++)
	{
		y_element=vector_elemets_y_cen[tl]/b;

		for (int tk=0; tk < N_1; tk++)
		{
			x_element=vector_elemets_x_cen[tk]/a;

			dist_quad=pow(x_element,2)+pow(y_element,2);

			if(dist_quad < 1)
			{
				tj=++tj;
			}
		}
	}


	NBE=tj;


	double ccp[NBE][2],pcp[NBE][2];

	tj=0;

	for (int tl=0; tl < N_2; tl++)
	{
		y_element=vector_elemets_y_cen[tl]/b;

		for (int tk=0; tk < N_1; tk++)
		{
			x_element=vector_elemets_x_cen[tk]/a;

			dist_quad=pow(x_element,2)+pow(y_element,2);

			if(dist_quad < 1)
			{
				ccp[tj][0]=vector_elemets_x_est[tk];
				ccp[tj][1]=vector_elemets_y_est[tl];

				pcp[tj][0]=vector_elemets_x_cen[tk];
				pcp[tj][1]=vector_elemets_y_cen[tl];

				tj=++tj;
			}
		}
	}


	BE = new BE3D [NBE];

	double c1=-z0;

	double xc,yc,zc,xp,yp,zp;

	double BE_L=passo_1;
	double BE_W=passo_2;

	double n_x,n_y,n_z;
	double ts_x,ts_y,ts_z;
	double td_x,td_y,td_z;

	double delta_phi=0e0;

	for(int i=0; i < NBE; i++)
	{
		xc =  ccp[i][0];
		yc =  ccp[i][1]*COS_DELTA;
		zc =  ccp[i][1]*SIN_DELTA-c1;

		xp =  pcp[i][0];
		yp =  pcp[i][1]*COS_DELTA;
		zp =  pcp[i][1]*SIN_DELTA-c1;

		n_x =  0;
		n_y = -SIN_DELTA;
		n_z =  COS_DELTA;

        ts_x = 1;
        ts_y = 0;
        ts_z = 0;

        td_x = 0;
        td_y = COS_DELTA;
        td_z = SIN_DELTA;

		BE[i].NEW(xc,yc,zc,xp,yp,zp,delta_phi,delta_gradi,BE_L,BE_W,n_x,n_y,n_z,ts_x,ts_y,ts_z,td_x,td_y,td_z);
	}

}





void SILL::PRINT(CONSOLE &out)
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





void SILL::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
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
		*out_stream << "a (major axis) = " << a << " m" << endl;
		*out_stream << "b (minor axis) = " << b << " m" << endl;

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
		*out_stream << "# " << "a (major axis) = " << a << " m" << endl;
		*out_stream << "# " << "b (minor axis) = " << b << " m" << endl;

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



void SILL::PRINT(DATAFILE &out)
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












