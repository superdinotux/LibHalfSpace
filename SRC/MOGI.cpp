/*
 * MOGI.cpp
 *
 *  Created on: 23/ott/2014
 *      Author: dino
 */

#include "MOGI.h"

#include <cmath>


#include "PRINT_TEMPLATE.h"


MOGI::MOGI() {
	// TODO Auto-generated constructor stub

}

MOGI::~MOGI() {
	// TODO Auto-generated destructor stub
}



void MOGI::SET(MEDIUM MEDIUM_PAR_i,double P_i)
{

	DP = P_i;

	MEDIUM_PAR = MEDIUM_PAR_i;

}



void MOGI::DISPLACEMENT(double r,double z,double &Ur,double &Uz)
//double Strength,double mu,double K)
{

	double mu = MEDIUM_PAR.GET_mu();
	double K = MEDIUM_PAR.GET_K();

	double d=-z0;

	double R1=sqrt(pow(r,2)+pow(z-d,2));
	double R1_ast=sqrt(pow(r,2)+pow(z+d,2));

	double Strength = DP * pow(R,3);

	double f_rz=((3*K+7*mu)/(3*K+mu))*(r/pow(R1_ast,3))-(6*z*r*(z+d))/pow(R1_ast,5);

	Ur =  (Strength/(4*mu))*((r/pow(R1,3))+f_rz);

	double g_rz=(-(3*K+7*mu)/(3*K+mu))*((z+d)/pow(R1_ast,3))+(2*z)/pow(R1_ast,3)-(6*z*pow(z+d,2))/pow(R1_ast,5);

	Uz = -(Strength/(4*mu))*(((z-d)/pow(R1,3))+g_rz);

}



void MOGI::DISPLACEMENT(double x,double y,double z,double U[])
{
	double xp=x-x0;
	double yp=y-y0;

	double r = sqrt(pow(xp,2)+pow(yp,2));

	double Uz,Ur;

	DISPLACEMENT(r,z,Uz,Ur);

	double teta=atan2(yp,xp);

	U[0]=Ur*cos(teta);
	U[1]=Ur*sin(teta);
	U[2]=Uz;
}




void MOGI::PRINT(CONSOLE &out)
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



void MOGI::PRINT_GEOM_PROP(ostream *out_stream)
{
	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;

	*out_stream << "Position of source center at:" << endl;
	*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
	*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
	*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

	*out_stream << endl;

	*out_stream << "Source radius = " << R << " m" << endl;

	*out_stream << endl;

}




void MOGI::PRINT_STRESS_PROP(ostream *out_stream)
{
	*out_stream << label << " - Stress properties" << endl;

	*out_stream << endl;

	*out_stream << "Source overpressure:" << endl;
	*out_stream << "DP = " << DP*1e-6 << " MPa" << endl;

	*out_stream << endl;
}





void MOGI::PRINT(DATAFILE &out)
{

	int flag_datatype;
	out.GET_flag_datatype(flag_datatype);

	if(flag_datatype == 0)
	{

		int int_flag;
		out.GET_int_flag(int_flag);

		string filename;

		int Nr;
		double ri,rf;

		out.GET_file_parameters(filename,Nr,ri,rf);

		filename = label + filename + ".dat";

		switch(int_flag)
		{
		case (1):
			print_map_radial_displ(*this,filename,Nr,ri,rf);

		break;

		default:
			cout << "Opzione " << int_flag << " non disponibile";
			break;
		}

	}
	else if(flag_datatype == 1)
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

}





void MOGI::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
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

		*out_stream << "# " << endl;

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;
	}
}






