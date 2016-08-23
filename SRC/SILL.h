/*
 * SILL.h
 *
 *  Created on: 07/nov/2014
 *      Author: dino
 */

#ifndef SILL_H_
#define SILL_H_

#include "SOURCE_BE3D.h"

class SILL: public SOURCE_BE3D
{

	double a;				// major axis
	double b;				// minor axis

	double theta_degree;	// inclination with respect to the interface

	void NEW(double a,double b,double delta_gradi);

    void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print = 0);

public:

	SILL();

	SILL(string s_label,double x0,double y0,double z0,double a,double b,double phi_degree,double theta_degree,int N)
			:SOURCE_BE3D(s_label,x0,y0,z0,phi_degree,N)
			{NEW(a,b,theta_degree);};

	virtual ~SILL();


	/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

    void PRINT(){SOURCE_BE3D::PRINT();};

	/*! Print to standard output useful information about the source */

    void PRINT(CONSOLE &out);

	/*! Generation of a datafile and a gnuplot script.
	 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

    void PRINT(DATAFILE &out);

};

#endif /* SILL_H_ */
