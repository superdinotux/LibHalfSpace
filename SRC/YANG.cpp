/*
 * YANG.cpp
 *
 *  Created on: 17/ott/2014
 *      Author: dino
 */

#include "YANG.h"

#include <cmath>
#include <iostream>


#include "PRINT_TEMPLATE.h"


YANG::YANG() {
	// TODO Auto-generated constructor stub

}

YANG::~YANG() {
	// TODO Auto-generated destructor stub
}



void YANG::NEW(double a_i,double b_i,double phi_i,double theta_i)
{
	if(a_i > b_i)
	{
		a = a_i*1e-3;
		b = b_i*1e-3;
	}
	else
	{
		cout << "Axis a must be greater than axis b." << endl;
		exit(1);
	}

	c = sqrt(pow(a,2)-pow(b,2));

	phi= phi_i;
	theta = theta_i;

	double factor = acos(-1)/180e0;

	phi_rad = factor*phi_i;
	theta_rad = factor*theta_i;
};



void YANG::SET(MEDIUM MEDIUM_PAR_i,double DP_i)
{

	DP = DP_i*1e-6;

	MEDIUM_PAR = MEDIUM_PAR_i;

	double lamda = MEDIUM_PAR.GET_lambda()*1e-7;
	double mu = MEDIUM_PAR.GET_mu()*1e-7;
	double nu = MEDIUM_PAR.GET_nu();

	double pi=acos(-1);

	// Model expressions
	double ac=(a-c)/(a+c);
	double L1=log(ac);
	//L1=log((a-c)/(a+c));

	double iia=2/a/pow(c,2) + L1/pow(c,3);
	double iiaa=2e0/3e0/pow(a,3)/pow(c,2) + 2/a/pow(c,4) + L1/pow(c,5);
	double coef1=-2*pi*a*pow(b,2);

	double Ia=coef1*iia;
	double Iaa=coef1*iiaa;

	double u=8*pi*(1-nu);
	double Q=3/u;
	double R=(1-2*nu)/u;

	double a11=2*R*(Ia-4*pi);
	double a12=-2*R*(Ia+4*pi);
	double a21=Q*pow(a,2)*Iaa + R*Ia - 1;
	double a22=-(Q*pow(a,2)*Iaa + Ia*(2*R-Q));

	double coef2=3*lamda+2*mu;
	double w=1/(a11*a22-a12*a21);

	double e11=(3*a22-a12)*DP*w/coef2;
	double e22=(a11-3*a21)*DP*w/coef2;

	double Pdila_i=2*mu*(e11-e22);
	double Pstar_i=lamda*e11 + 2*(lamda+mu)*e22;

	Pstar = Pstar_i;
	Pdila = Pdila_i;

	double a1_i=-2*pow(b,2)*Pdila;
	double b1_i= 3*pow(b,2)*Pdila/pow(c,2) + 2*(1-2*nu)*Pstar; // !PL version had (1-nu) in the 2nd term!

	a1 = a1_i;
	b1 = b1_i;

}





template <typename T> T sign(T t);


void YANG::DISPL(double x,double y,double z,double U[],double xi)
{

	// Calculate the double force (star) and dilatation (dila) displacements U
	// for a SPHEROIDAL pressure source in an elastic halfspace
	// (based on Yang et al., vol 93, JGR, 4249-4257, 1988) with arbitrary plunge (theta)
	// of the long axis of the spheroid (theta = 90, vertical; theta = 0, horizontal).
	// Evaluate at for xi.
	//
	// Inputs: theta: dip angle of source
	//             P: pressure change in magma chamber
	//             a: semimajor axis of spheroid
	//             b: semiminor axis of spheriod
	//            xi: evaluate integrals at +- c
	//            z0: depth of source (allowed to vary with topo)
	//             x: x location of point on surface
	//             y: y location of point on surface
	// Output: rd: calculated range displacement
	// NOTE: the x, y locations assume a source at origin
	// ALSO: the units need to be in mks units so input x, y, and z0
	//       in km will be changed into meters
	// NOTE: In the expressions of Yang et al. the spheroid dips parallel to the y axis
	//       at x=0. We will assume (initially) that it is also centered at y=0.

	//double epsn=1e-15;

	//	int dim = x.size();

	double pi=acos(-1);

	double tp=0;

	//Poisson's ratio, Young's modulus, and the Lame coeffiecents mu and lamda
	double nu = MEDIUM_PAR.GET_nu();

	double nu4 = 3-4*nu;
	double nu1 = 1-nu;

	double mu = MEDIUM_PAR.GET_mu()*1e-7;

	double coeffs_0 = 1e0/(16*mu*(1-nu));
	double coeffs_1 = nu4;
	double coeffs_2 = 4*(1-nu)*(1-2*nu);

	double coeff=a*pow(b,2)/pow(c,3)*coeffs_0;


	double sinth = sin(theta_rad);
	double costh = cos(theta_rad);


	// Introduce new coordinates and parameters (Yang et al., 1988, page 4251):
	double xi2=xi*costh;
	double xi3=xi*sinth;
	double y0=0;


	double z00;

	double x1 = x;
	double y1 = x1;

	double x2;
	double x3;

	double y2;
	double y3;

	double xbar3;
	double ybar3;

	double Z_0=-z0*1e-3;

//	z00 = tp+z0;
	z00 = tp+Z_0;

	x2 = y-y0;
	x3 = z-z00;

	y2 = x2-xi2;
	y3 = x3-xi3;

	xbar3 = z+z00;
	ybar3 = xbar3+xi3;

	double r2 =  x2*sinth - x3*costh;
	double q2 =  x2*sinth + xbar3*costh;
	double r3 =  x2*costh + x3*sinth;
	double q3 = -x2*costh + xbar3*sinth;


	double rbar3;
	double qbar3;

	double R1,R2;

	double C0;

	double betatop;
	double betabottom;

	rbar3 = r3-xi;
	qbar3 = q3+xi;

	R1 = sqrt(pow(y1,2)+pow(y2,2)+pow(y3,2));
	R2 = sqrt(pow(y1,2)+pow(y2,2)+pow(ybar3,2));

	C0 = y0*costh + z00*sinth;  // check this

	betatop = costh*q2+(1+sinth)*(R2+qbar3);
	betabottom = costh*y1;

	double atnbeta;

	if(abs(betabottom == 0))
	{
		atnbeta=pi/2*sign(betatop);
	}
	else
	{
		atnbeta=atan(betatop/betabottom);
	}


	// Set up other parameters for dipping spheroid (Yang et al., 1988, page 4252):
	// precalculate some repeatedly used natural logs:

	double Rr = R1 + rbar3;
	double Rq = R2 + qbar3;
	double Ry = R2 + ybar3;

	double lRr = log(Rr);
	double lRq = log(Rq);
	double lRy = log(Ry);


	double A1star = a1/(R1*Rr) + b1*(lRr+(r3+xi)/Rr);

	double Abar1star = -a1/(R2*Rq) - b1*(lRq+(q3-xi)/Rq);

	double A1 = xi/R1 + lRr;

	double Abar1 = xi/R2 - lRq;

	double A2 = R1-r3*lRr;

	double Abar2 = R2 - q3*lRq;

	double A3 = xi*rbar3/R1 + R1;

	double Abar3 = xi*qbar3/R2 -R2;


	double B = xi*(xi+C0)/R2 - Abar2 - C0*lRq;

	double Bstar = a1/R1 + 2*b1*A2 + coeffs_1*(a1/R2 + 2*b1*Abar2);


	double F1=0;
	double F1star=0;
	double F2=0;
	double F2star=0;

	// Skip if displacement calculated at surface (z=0)
	if(z != 0)
	{
		F1 = -2*sinth*z*(xi*(xi+C0)/pow(R2,3) + (R2+xi+C0)/(R2*(Rq)) + 4*(1-nu)*(R2+xi)/(R2*(Rq)));

		F1star = 2*z*(costh*q2*(a1*(2*Rq)/(pow(R2,3)*pow(Rq,2))
				- b1*(R2+2*xi)/(R2*pow(Rq,2)))+ sinth*(a1/pow(R2,3) -2*b1*(R2+xi)/(R2*(Rq))));

		F2 = -2*sinth*z*(xi*(xi+C0)*qbar3/pow(R2,3) + C0/R2 + (5-4*nu)*Abar1);

		F2star = 2*z*(a1*ybar3/pow(R2,3) - 2*b1*(sinth*Abar1 + costh*q2*(R2+xi)/(R2*Rq)));
	}


	// calculate little f's

	double ff1 = xi*y1/Ry + 3/pow(costh,2)*(y1*lRy*sinth - y1*lRq
			+ 2*q2*atnbeta) +	2*y1*lRq - 4*xbar3*atnbeta/costh;

	double ff2 = xi*y2/Ry + 3/pow(costh,2)*(q2*lRq*sinth - q2*lRy + 2*y1*atnbeta*sinth
			+ costh*(R2-ybar3)) - 2*costh*Abar2 + 2/costh*(xbar3*lRy - q3*lRq);

	double ff3 = (q2*lRq - q2*lRy*sinth + 2*y1*atnbeta)/costh + 2*sinth*Abar2 + q3*lRy - xi;


	// Assemble into x, y, z displacements (1,2,3):

	U[0] = coeff*(A1star + nu4*Abar1star + F1star)*y1;

	U[1] = coeff*(sinth*(A1star*r2+(nu4*Abar1star+F1star)*q2)
			+ costh*(Bstar-F2star) + 2*sinth*costh*z*Abar1star);

	U[2] = coeff*(-costh*(Abar1star*r2+(nu4*Abar1star-F1star)*q2)
			+ sinth*(Bstar+F2star) + 2*pow(costh,2)*z*Abar1star);


	U[0] += 2*coeff*Pdila*((A1 + nu4*Abar1 + F1)*y1 - coeffs_2*ff1);

	U[1] += 2*coeff*Pdila*(sinth*(A1*r2+(nu4*Abar1+F1)*q2)
			- coeffs_2*ff2 + 4*nu1*costh*(A2+Abar2) + costh*(A3 - nu4*Abar3 - F2));

	U[2] += 2*coeff*Pdila*(costh*(-A1*r2 + (nu4*Abar1 + F1)*q2)
			+ coeffs_2*ff3 + 4*nu1*sinth*(A2+Abar2) + sinth*(A3 + nu4*Abar3 + F2 - 2*nu4*B));

}




template <typename T>
T sign(T t)
{
	if(t == 0)
	{
		return T(0);
	}
	else if(t <0)
	{
		return T(-1);
	}
	else
	{
		return T(+1);
	}
}





void YANG::DISPLACEMENT(double x,double y,double z,double U[])
{

	// Calculate range displacement

	double xn = (x-x0)*1e-3;
	double yn = (y-y0)*1e-3;

	double cosp = cos(phi_rad);
	double sinp = sin(phi_rad);


	// Rotate points

	double xp = xn*cosp + yn*sinp;
	double yp = yn*cosp - xn*sinp;


	// Calculate model at integration limits

	double xi;

	double Up[3];

	xi=c;

	DISPL(xp,yp,z,Up,xi);


	double Um[3];

	xi=-xi;

	DISPL(xp,yp,z,Um,xi);


	// Sum

	double Ur[3];

	for(int i=0; i < 2; i++)
	{
		Ur[i] = -Up[i] + Um[i];
	}

	// Rotate horiz. displacements back to the orig. coordinate system:

	U[0] = (Ur[0]*cosp - Ur[1]*sinp)*1e2;
	U[1] = (Ur[0]*sinp + Ur[1]*cosp)*1e2;
	U[2] = (Up[2]      - Um[2])*1e2;

}





void YANG::PRINT(CONSOLE &out)
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



void YANG::PRINT_GEOM_PROP(ostream *out_stream)
{
	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;

	*out_stream << "Position of source center at:" << endl;
	*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
	*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
	*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

	*out_stream << endl;

	*out_stream << "Axis dimensions:" << endl;
	*out_stream << "a = " << a << " km" << endl;
	*out_stream << "b = " << b << " km" << endl;
	*out_stream << "c = " << c << " km" << endl;

	*out_stream << endl;

	*out_stream << "Source orientation:"   << endl;
	*out_stream << "phi = "   << phi   << endl;
	*out_stream << "theta = " << theta << endl;

	*out_stream << endl << endl;

}




void YANG::PRINT_STRESS_PROP(ostream *out_stream)
{
	*out_stream << label << " - Stress properties" << endl;

	*out_stream << endl;

	*out_stream << "Source overpressure:" << endl;
	*out_stream << "DP = " << DP << " MPa" << endl;

	*out_stream << endl << endl;
}




void YANG::PRINT(DATAFILE &out)
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
		cout << endl << label << ": unsymmetric source, datafile object is not compatible with print method" << endl;
	}
}






void YANG::PRINT_MAIN_FEATURES(ostream *out_stream, int flag_print)
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

		*out_stream << "Axis dimensions:" << endl;
		*out_stream << "a = " << a << " km" << endl;
		*out_stream << "b = " << b << " km" << endl;
		*out_stream << "c = " << c << " km" << endl;

		*out_stream << endl;

		*out_stream << "Source orientation:"   << endl;
		*out_stream << "phi = "   << phi   << endl;
		*out_stream << "theta = " << theta << endl;

		*out_stream << endl << endl;

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


		*out_stream << "# " << "Axis dimensions:" << endl;
		*out_stream << "# " << "a = " << a << " km" << endl;
		*out_stream << "# " << "b = " << b << " km" << endl;
		*out_stream << "# " << "c = " << c << " km" << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "Source orientation:"   << endl;
		*out_stream << "# " << "phi = "   << phi   << endl;
		*out_stream << "# " << "theta = " << theta << endl;

		*out_stream << "# " << endl << endl;

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;

	}

}


