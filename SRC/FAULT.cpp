
#include "FAULT.h"

#include <cmath>
#include <string>
#include <iomanip>

#include "BE3D.h"
#include "PRINT_TEMPLATE.h"



FAULT::FAULT(){
	// TODO Auto-generated constructor stub
}

FAULT::~FAULT() {
	// TODO Auto-generated destructor stub
}




void FAULT::NEW(double Li,double Wi,double delta_gradi)
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

	int N_1;
	int N_2;

	L=Li;
	W=Wi;

	double LATO_h_1=L;
	double LATO_h_2=W;

	double passo_1;
	double passo_2;

	double cfr[NUM];
	double min;


	if((L == W) || (NUM == 1))
	{
		N_1=NUM;
		N_2=NUM;

		passo_1=LATO_h_1/N_1;
		passo_2=LATO_h_2/N_2;
	}
	else if(L > W)
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


	NBE=N_1*N_2;
	BE = new BE3D [NBE];

	GRID_PAR.N_1=N_1;
	GRID_PAR.N_2=N_2;

	double ccp[2];
	double pcp[2];

    double xccI=-L*0.5;
    double yccI=-W*0.5;

    double xpcI=-(N_1-1)*passo_1*0.5;
    double ypcI=-(N_2-1)*passo_2*0.5;

    double xcc,ycc;
    double xpc,ypc;

	double c1=-z0;

    double xc,yc,zc;
    double xp,yp,zp;

    double BE_L=L/N_1;
    double BE_W=W/N_2;

    double n_x,n_y,n_z;
	double ts_x,ts_y,ts_z;
	double td_x,td_y,td_z;

	double phi_angle=0e0;


    int tj;

    for(int tx=0; tx < N_1; tx++)
    {

        xcc=xccI+passo_1*tx;
        xpc=xpcI+passo_1*tx;

        for(int ty=0; ty < N_2; ty++)
        {

            ycc=yccI+passo_2*ty;
            ypc=ypcI+passo_2*ty;

            tj=tx*N_2+ty;

            ccp[0]=xcc;
            ccp[1]=ycc;

            pcp[0]=xpc;
            pcp[1]=ypc;

    		xc =  ccp[0];
    		yc =  ccp[1]*COS_DELTA;
    		zc =  ccp[1]*SIN_DELTA-c1;

    		xp =  pcp[0];
    		yp =  pcp[1]*COS_DELTA;
    		zp =  pcp[1]*SIN_DELTA-c1;

            n_x =  0;
            n_y = -SIN_DELTA;
            n_z =  COS_DELTA;

            ts_x =  1;
            ts_y =  0;
            ts_z =  0;

            td_x =  0;
            td_y =  COS_DELTA;
            td_z =  SIN_DELTA;

    	   	BE[tj].NEW(xc,yc,zc,xp,yp,zp,phi_angle,delta_gradi,BE_L,BE_W,n_x,n_y,n_z,ts_x,ts_y,ts_z,td_x,td_y,td_z);

        }

    }

}



void FAULT::PRINT(CONSOLE &out)
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

		case (2):
            	PRINT_GEOM_PROP(stream);
				break;

		case (3):
				PRINT_STRESS_PROP(stream);
				break;

		case (4):
				PRINT_BV_COMP(stream);
				break;

		case (10):
			    PRINT_MAIN_FEATURES(stream,flag_print);
				break;

		default:
				*stream << "Not implemented for this source type" << endl;
	}


}



void FAULT::PRINT_GEOM_PROP(ostream *out_stream)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);

	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;
	*out_stream << setw(7)  << "i="            ;
	*out_stream << setw(6)  << "j="            ;
	*out_stream << setw(15) << "x0 (m)"        ;
	*out_stream << setw(14) << "y0 (m)"        ;
	*out_stream << setw(14) << "z0 (m)"        ;
	*out_stream << setw(18) << "angle (degree)";
	*out_stream << endl;

	BE3D_PARAM PARAM;

	int cont;

	for(int i=0; i < GRID_PAR.N_1; i++)
	{
		for(int j=0; j < GRID_PAR.N_2; j++)
		{
			cont=i*(GRID_PAR.N_2)+j;

			if(j == 0)
			{
				*out_stream << setw(8) << i+1;
				*out_stream << setw(6) << j+1;
			}
			else
			{
				*out_stream << setw(14) << j+1;
			}


			BE[cont].GET_PARAM(PARAM);

			double x1=PARAM.GEOM.pc[0];
			double y1=PARAM.GEOM.pc[1];

			double fi=(acos(-1)/180e0)*phi_degree;

			double cos_fi=cos(fi);
			double sin_fi=sin(fi);

			double x2=x1*cos_fi-y1*sin_fi;
			double y2=x1*sin_fi+y1*cos_fi;

			*out_stream << setw(14) <<  x2+x0;
			*out_stream << setw(14) <<  y2+y0;
			*out_stream << setw(14) <<  PARAM.GEOM.pc[2];
			*out_stream << setw(15) <<  PARAM.GEOM.delta_degree;
			*out_stream << endl;
		}

		*out_stream << endl;
	}

	*out_stream << endl << endl;

}



void FAULT::PRINT_STRESS_PROP(ostream *stream)
{

	double BC1,BC2,BC3;

	stream->fill(' ');
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	*stream << label << " - Stress properties"<< endl << endl;

	*stream << endl;
	*stream << setw(7)  << "i="            ;
	*stream << setw(6)  << "j="            ;
	*stream << setw(15) << "BC1 (MPa)"     ;
	*stream << setw(14) << "BC2 (MPa)"     ;
	*stream << setw(14) << "BC3 (MPa)"     ;
	*stream << endl;

	int cont;

	for(int i=0; i < GRID_PAR.N_1; i++)
	{
		for(int j=0; j < GRID_PAR.N_2; j++)
		{
			cont=i*(GRID_PAR.N_2)+j;

			if(j == 0)
			{
				*stream << setw(8) << i+1;
				*stream << setw(6) << j+1;
			}
			else
			{
				*stream << setw(14) << j+1;
			}

			BE[cont].GET_BC(BC1,BC2,BC3);

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

		*stream << endl;
	}

	*stream << endl << endl;

}



void FAULT::PRINT_BV_COMP(ostream *stream)
{

	double B_1,B_2,B_3;

	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->precision(5);

	*stream << label << " - Burgers vector components" << endl << endl;

	*stream << endl;
	*stream << setw(7)  << "i="            ;
	*stream << setw(6)  << "j="            ;
	*stream << setw(16) << "S (m)" ;
	*stream << setw(16) << "D (m)" ;
	*stream << setw(16) << "T (m)" ;
	*stream << endl;

	int cont;

	for(int i=0; i < GRID_PAR.N_1; i++)
	{
		for(int j=0; j < GRID_PAR.N_2; j++)
		{
			cont=i*(GRID_PAR.N_2)+j;

			if(j == 0)
			{
				*stream << noshowpos << setw(8) << i+1;
				*stream << setw(6) << j+1;
			}
			else
			{
				*stream << noshowpos << setw(14) << j+1;
			}

			BE[cont].GET_BV(B_1,B_2,B_3);

			stream->fill(' ');
			*stream << setw(6) << " ";

			stream->fill(' ');
			*stream << right << showpos << setw(8) << B_1 ;

			stream->fill(' ');
			*stream << setw(4) << " ";

			stream->fill(' ');
			*stream << right << setw(8) << B_2 ;
			stream->fill(' ');
			*stream << setw(4) << " ";

			stream->fill(' ');
			*stream << right << setw(8) << B_3 ;

			*stream << endl;

		}

		*stream << endl;

	}

	*stream << endl << endl;

}






void FAULT::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
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
		*out_stream << "# " << "L = " << L << " m" << endl;
		*out_stream << "# " << "W = " << W << " m" << endl;

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




void FAULT::PRINT(DATAFILE &out)
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



