/*
 * SOURCE.h
 *
 *  Created on: 07/lug/2014
 *      Author: dino
 */

#ifndef SOURCE_H_
#define SOURCE_H_

///*! \brief Brief description.
// *         Brief description continued.
// *
// *  Detailed description starts here.
// */

#include <string>
#include <iostream>
#include <stdlib.h>
#include <iomanip>


#include "CONSOLE.h"
#include "DATAFILE.h"
#include "MEDIUM.h"

using namespace std;

class SOURCE{

protected: 			// All PROTECTED MEMBERS of the CLASS are DATA

	string label;
	double x0,y0,z0;
	MEDIUM MEDIUM_PAR;

public:				// All  PUBLIC MEMBERS of the CLASS are METHODS

	SOURCE();

	SOURCE(string label_i, double x0_i, double y0_i, double z0_i)
	{
		label = label_i; x0 = x0_i; y0 = y0_i;

		if(z0_i <= 0)
			{z0 = z0_i;}
		else
			{cout << "Source depth cannot be positive." << endl; exit(1);}
	};

	virtual ~SOURCE();

/*! Association of the medium param with the source buried in the Half-space  */

	virtual void SET(MEDIUM MEDIUM_PAR_i){MEDIUM_PAR = MEDIUM_PAR_i;};

/*! Computation of the displacement components (Ux=U[0], Uy=U[1], Uz=U[2]) at the observation point (x,y,z) */

	virtual void DISPLACEMENT(double x,double y,double z,double U[])
						{cout << "Method not implemented for this type of source" << endl;};

/*! Computation of the displacement components (Ur radial component, Uz vertical component]) at the observation point (x,y,z) (available only for axially symmetric sources) */

	virtual void DISPLACEMENT(double r,double z,double &Ur,double &Uz)
						{cout << "Method not implemented for this type of source" << endl;};

/*! Computation of the strain tensor components
 * (Sxx=Strain[0], Syy=Strain[1], Szz=Strain[2],Sxy=Strain[3],Sxz=Strain[4], Syz=Strain[5])
 *  at the observation point (x,y,z) */

	virtual void STRAIN(double x,double y,double z,double Strain[])
						{cout << "Method not implemented for this type of source" << endl;};

/*! Computation of the stress tensor components
 * (Sxx=Stress0], Syy=Stress[1], Szz=Stress[2],Sxy=Stress[3],Sxz=Stress[4], Syz=Stress[5])
 *  at the observation point (x,y,z) */

	virtual void STRESS(double x,double y,double z,double Stress[])
						{cout << "Method not implemented for this type of source" << endl;};


/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

	virtual void PRINT(){cout << "Method not implemented for this type of source" << endl;};

/*! Print to standard output useful information about the source */

	virtual void PRINT(CONSOLE &out)
		{cout << "Source center is located at point " << x0*1e-3 << "," << y0*1e-3 << "," << z0*1e-3 << " (km)" << endl;};

/*! Generation of a datafile and a gnuplot script.
 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

	virtual void PRINT(DATAFILE &out) = 0;

/*! We retrieve the label associated to the source */

	virtual void GET_label(string &label_o){label_o = label;};

/*! We retrieve the coordinates of the source center (x0,y0,z0) */

	virtual void GET_x0_y0_z0(double &x0_o, double &y0_o, double &z0_o){x0_o = x0; y0_o = y0; z0_o = z0;};

/*! We retrieve the MEDIUM_PAR object which describes the properties of the medium where the source is buried */

	virtual void GET_MEDIUM_PAR(MEDIUM &MEDIUM_PAR_o){MEDIUM_PAR_o = MEDIUM_PAR;};

/*! We change the label associated to the source */

	virtual void PUT_label(string label_i){label = label_i;};

/*! We change the coordinates of the source center (x0,y0,z0) */

	virtual void PUT_x0_y0_z0(double x0_i, double y0_i, double z0_i){x0 = x0_i; y0 = y0_i; z0 = z0_i;};

/*! We change the MEDIUM_PAR object which describes the properties of the medium where the source is buried */

	virtual void PUT_MEDIUM_PAR(MEDIUM MEDIUM_PAR_i){MEDIUM_PAR = MEDIUM_PAR_i;};

};


#endif
