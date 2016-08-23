/*
 * FAULT.h
 *
 *  Created on: 18/lug/2014
 *      Author: dino
 */

#ifndef FAULT_H_
#define FAULT_H_


#include <typeinfo>
#include "BE3D.h"
#include "CONSOLE.h"
#include "DATAFILE.h"
#include "GRID.h"
#include "SOURCE_BE3D.h"


class FAULT: public SOURCE_BE3D
{
protected:

	double L;
	double W;

	double theta_degree;


	void NEW(double L,double W,double delta_gradi);

    void PRINT_GEOM_PROP(ostream *out_stream);
    void PRINT_STRESS_PROP(ostream *stream);
    void PRINT_BV_COMP(ostream *stream);

    void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print = 0);

public:

    FAULT();

	FAULT(string s_label,double x0,double y0,double z0,double L,double W,double phi_degree,double theta_degree,int N)
			:SOURCE_BE3D(s_label,x0,y0,z0,phi_degree,N)
			{NEW(L,W,theta_degree);};

	virtual ~FAULT();

	/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

    void PRINT(){SOURCE_BE3D::PRINT();};

	/*! Print to standard output useful information about the source */

    void PRINT(CONSOLE &out);

	/*! Generation of a datafile and a gnuplot script.
	 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

    void PRINT(DATAFILE &out);
};


#endif
