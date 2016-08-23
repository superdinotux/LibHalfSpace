/*
 * YANG.h
 *
 *  Created on: 17/ott/2014
 *      Author: dino
 */

#ifndef YANG_H_
#define YANG_H_

#include <string>
#include <cmath>


#include "PRESSURIZED.h"


class YANG: public PRESSURIZED
{

	double a;
	double b;

	double phi;
	double theta;


	double c;

	double phi_rad;
	double theta_rad;

	double Pstar;
	double Pdila;
	double a1;
	double b1;

    void NEW(double a_i,double b_i,double phi_i,double theta_i);

	void DISPL(double x,double y,double z,double U[],double xi);

    void PRINT_GEOM_PROP(ostream *out_stream);
    void PRINT_STRESS_PROP(ostream *stream);

    void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print=0);

public:

    YANG();

	YANG(string label_i,double x0_i,double y0_i,double z0_i,double a_i,double b_i,double phi_i,double theta_i):
				PRESSURIZED(label_i,x0_i,y0_i,z0_i)
				{NEW(a_i,b_i,phi_i,theta_i);}

	virtual ~YANG();

	/*! Association of the medium param with the source buried in the Half-space. It also assigns the overpressure.  */

	void SET(MEDIUM MEDIUM_PAR_i,double P_i);

	/*! Computation of the displacement components (Ux=U[0], Uy=U[1], Uz=U[2]) at the observation point (x,y,z) */

	void DISPLACEMENT(double x,double y,double z,double U[]);

	/*! Print to standard output useful information about the source */

	void PRINT(CONSOLE  &out);

	/*! Generation of a datafile and a gnuplot script.
	 *  According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface. */

	void PRINT(DATAFILE &out);


};



#endif /* YANG_H_ */
