/*
 * Utilities.h
 *
 *  Created on: 21/ott/2014
 *      Author: dino
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_


#include "cmath"


template <class type>
void Rotazione_coord(type &x,type &y,int flag)
{
      float cv;

      if(flag == 1)
      {
		  cv=x;
		  x=y;
		  y=-cv;
      }
};


template <class type>
void Rotazione_coord(type &x,type &y,double phi_degree)
{
	double pi=acos(-1);

	double phi = pi*phi_degree/180e0;

	double x1 =  cos(phi)*x - sin(phi)*y;
	double y1 = +sin(phi)*x + cos(phi)*y;

	x=x1;
	y=y1;
};



template <class type>
void Rotazione_vettore(type V[],double phi_degree)
{
	double pi=acos(-1);

	double phi = pi*phi_degree/180e0;

	double V1[2];

	V1[0] =  cos(phi)*V[0] + sin(phi)*V[1];
	V1[1] = -sin(phi)*V[0] + cos(phi)*V[1];

	V[0] = V1[0];
	V[1] = V1[1];
};



template <class type>
void Rotazione_tensore(type S[],double phi_degree)
{
	double pi=acos(-1);

	double S1[6];

	double phi = pi*phi_degree/180e0;

	double cos_phi = cos(phi);
	double sin_phi = sin(phi);

	double cos_phi_2 = pow(cos_phi,2);
	double sin_phi_2 = pow(sin_phi,2);

	double cos_2_phi = cos(2*phi);
	double sin_2_phi = sin(2*phi);

	S1[0] =  S[0]*cos_phi_2 + S[1]*sin_phi_2 + S[3]*sin_2_phi;
	S1[1] =  S[0]*sin_phi_2 + S[1]*cos_phi_2 - S[3]*sin_2_phi;
	S1[2] =  S[2];
	S1[3] =  S[3]*cos_2_phi + 0.5*(S[1]-S[0])*sin_2_phi;
	S1[4] =  S[4]*cos_phi + S[5]*sin_phi;
	S1[5] = -S[4]*sin_phi + S[5]*cos_phi;

	for(int i=0; i < 6; i++)
	{
		S[i] = S1[i];
	}
};


#define ERRTOL 0.05
#define TINY 1.0e-25
#define BIG 4.5e21
#define C1 (3.0/14.0)
#define C2 (1.0/6.0)
#define C3 (9.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.25*C3)
#define C6 (1.5*C4)


static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


float rd(float x, float y, float z);


void Cycle_counter(int jf,int nsg,int c[]);

void gsl_eigenvalueproblem(int dimension,double eigenvalues[],double Matrix[]);

void ref_func(double &x,double &y,double &z,double X1,double X2,double X3_plane,int flag_plane);


void DC3D(float alpha,float X,float Y,float Z,float DEPTH,float DIP,float AL1,float AL2,float AW1,float AW2,float DISL1,float DISL2,float DISL3,float &UX,float &UY,float &UZ,float &UXX,float &UYX,float &UZX,float &UXY,float &UYY,float &UZY,float &UXZ,float &UYZ,float &UZZ,int &IRET,int CASE);






#endif /* UTILITIES_H_ */
