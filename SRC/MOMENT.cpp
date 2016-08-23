/*
 * MOMENT.cpp
 *
 *  Created on: 24/ott/2014
 *      Author: dino
 */

#include "MOMENT.h"

#include "cmath"


#include "PRINT_TEMPLATE.h"


MOMENT::MOMENT() {
	// TODO Auto-generated constructor stub

}

MOMENT::~MOMENT() {
	// TODO Auto-generated destructor stub
}





double G_hs(int tj,int tk,int ti,double X,double Y,double xi,double mu,double nu);

void MOMENT::DISPLACEMENT(double x,double y,double z,double U[])
{

	// M=[
	//    8.6835e+16,   6.7210e+02,   9.5627e+10;
	//    6.7210e+02,   8.6835e+16,   9.5627e+10;
	//    9.5627e+10,   9.5627e+10,  2.4692e+17
	// ]

	if(z != 0)
	{
		cout << "z != 0: we can compute only ground displacement" << endl;
		exit(1);
	}

	double X = x-x0;
	double Y = y-y0;
	double xi = -z0;

	double mu = MEDIUM_PAR.GET_mu();
	double nu = MEDIUM_PAR.GET_nu();


	double W_12_1 = G_hs(0,1,0,X,Y,xi,mu,nu);
	double W_13_1 = G_hs(0,2,0,X,Y,xi,mu,nu);
	double W_23_1 = G_hs(1,2,0,X,Y,xi,mu,nu);

	U[0] = M[0]*G_hs(0,0,0,X,Y,xi,mu,nu) + M[1]*G_hs(1,1,0,X,Y,xi,mu,nu) + M[2]*G_hs(2,2,0,X,Y,xi,mu,nu)
    				+ M[3]*W_12_1 + M[4]*W_13_1 + M[5]*W_23_1;


	double W_12_2 = G_hs(0,1,1,X,Y,xi,mu,nu);
	double W_13_2 = G_hs(0,2,1,X,Y,xi,mu,nu);
	double W_23_2 = G_hs(1,2,1,X,Y,xi,mu,nu);

	U[1] = M[0]*G_hs(0,0,1,X,Y,xi,mu,nu) + M[1]*G_hs(1,1,1,X,Y,xi,mu,nu)  +M[2]*G_hs(2,2,1,X,Y,xi,mu,nu)
    				+ M[3]*W_12_2 + M[4]*W_13_2 + M[5]*W_23_2;


	double W_12_3 = G_hs(0,1,2,X,Y,xi,mu,nu);
	double W_13_3 = G_hs(0,2,2,X,Y,xi,mu,nu);
	double W_23_3 = G_hs(1,2,2,X,Y,xi,mu,nu);

	U[2] = -( M[0]*G_hs(0,0,2,X,Y,xi,mu,nu) + M[1]*G_hs(1,1,2,X,Y,xi,mu,nu) + M[2]*G_hs(2,2,2,X,Y,xi,mu,nu)
    				+ M[3]*W_12_3 + M[4]*W_13_3 + M[5]*W_23_3 );


}



double G_hs(int tj,int tk,int ti,double X,double Y,double xi,double mu,double nu)
{

	double x[2];

	x[0]=X;
	x[1]=Y;

	double R=sqrt(pow(X,2)+pow(Y,2)+pow(xi,2));
	double alpha=(3*R+xi)/(pow(R,3)*pow(R+xi,3));
	double beta=(2*R+xi)/(pow(R,3)*pow(R+xi,2));
	double eta=1e0/(R*pow(R+xi,2));
	double zeta=1e0/(R*(R+xi));

	double delta_i1,delta_i2,delta_ij;

	double U_jj_i,U_jj_3;
	double U_33_i,U_33_3;

	double W_12_i,W_12_3;
	double W_j3_i,W_j3_3;

	double result;

	double pi = acos(-1);


	if((ti <= 2) && (tj <= 2) && (tk <= 2))
	{

		if(tj == tk)
		{

			if((ti != 2) && (tj != 2))
			{

				if(ti == tj)
				{
					delta_ij=1e0;
				}
				else
				{
					delta_ij=0e0;
				}

				U_jj_i=(x[ti]/(4e0*pi*mu))*(-1e0/pow(R,3)+3e0*pow(x[tj],2)/pow(R,5)+(1e0-2e0*nu)*((1e0+2e0*delta_ij)*eta-pow(x[tj],2)*alpha));
				result=U_jj_i;
			}
			else if((ti == 2) && (tj != 2))
			{
				U_jj_3=(1e0/(4e0*pi*mu))*(xi/pow(R,3)-3*pow(x[tj],2)*xi/pow(R,5)-(1e0-2e0*nu)*(zeta-pow(x[tj],2)*beta));
				result=U_jj_3;
			}
			else if((ti != 2) && (tj == 2))
			{
				U_33_i=(x[ti]/(4e0*pi*mu))*(3e0*pow(xi,2)/pow(R,5)-2e0*nu/pow(R,3));
				result=U_33_i;
			}
			else if((ti == 2) && (tj == 2))
			{
				U_33_3=(xi/(4e0*pi*mu))*(-3e0*pow(xi,2)/pow(R,5)+2e0*nu/pow(R,3));
				result=U_33_3;
			}

		}
		else
		{

			if((tj == 0) && (tk == 1))
			{

				if(ti != 2)
				{

					if(ti == 1)
					{
						delta_i2=1e0;
					}
					else
					{
						delta_i2=0e0;
					}

					if(ti == 0)
					{
						delta_i1=1e0;
					}
					else
					{
						delta_i1=0e0;
					}

					W_12_i=((x[1]*(1e0-delta_i2)+x[0]*(1e0-delta_i1))/(2e0*pi*mu))*(3e0*pow(x[ti],2)/pow(R,5)+(1e0-2e0*nu)*(eta-pow(x[ti],2)*alpha));
					result=W_12_i;
				}
				else
				{
					W_12_3=(x[0]*x[1]/(2e0*pi*mu))*(-3e0*xi/pow(R,5)+(1e0-2e0*nu)*beta);
					result=W_12_3;
				}
			}
			else if(tk == 2)
			{
				if(ti != 2)
				{
					W_j3_i=(x[ti]*x[tj]/(2e0*pi*mu))*(-3e0*xi/pow(R,5));
					result=W_j3_i;
				}
				else
				{
					W_j3_3=(x[tj]*xi/(2e0*pi*mu))*(3e0*xi/pow(R,5));
					result=W_j3_3;
				}

			}

		}

	}
	else
	{
		cout << "error" << endl;
	}

	return result;

}





void MOMENT::PRINT(CONSOLE &out)
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
		break;

		case (1):
				PRINT_GEOM_PROP(stream);
		break;

		case (2):
				PRINT_STRESS_PROP(stream);
		break;

		case (10):
			    PRINT_MAIN_FEATURES(stream,flag_print);
				break;

		default:
				*stream << "Not implemented for this source type" << endl;
	}

}



void MOMENT::PRINT_GEOM_PROP(ostream *out_stream)
{
	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;

	*out_stream << "Position of source center at:" << endl;
	*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
	*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
	*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

	*out_stream << endl;

}




void MOMENT::PRINT_STRESS_PROP(ostream *out_stream)
{
	*out_stream << label << " - Moment tensor components" << endl;

	*out_stream << endl;

	*out_stream << "Mxx = " << M[0] << endl;
	*out_stream << "Myy = " << M[1] << endl;
	*out_stream << "Mzz = " << M[2] << endl;
	*out_stream << "Mxy = " << M[3] << endl;
	*out_stream << "Mxz = " << M[4] << endl;
	*out_stream << "Myz = " << M[5] << endl;

	*out_stream << endl;
}




void MOMENT::PRINT(DATAFILE &out)
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
			case (1):
            	print_map_displ(*this,filename,N_1,N_2,DX1,DX2);
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






void MOMENT::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
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

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;
	}
}





