/*
 * FIALKO.h
 *
 *  Created on: 20/ott/2014
 *      Author: dino
 */

#ifndef FIALKO_H_
#define FIALKO_H_

#include <vector>


using namespace std;


#include "PRESSURIZED.h"


class FIALKO: public PRESSURIZED{

	double R;
	double h;

	vector<double> fi,psi;
	vector<double> t,Wt;

    void PRINT_GEOM_PROP(ostream *out_stream);
    void PRINT_STRESS_PROP(ostream *stream);

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print=0);

public:

	FIALKO();

	FIALKO(string label_i,double x0_i,double y0_i,double z0_i,double R_i):
		PRESSURIZED(label_i,x0_i,y0_i,z0_i)
		{R = R_i; h = (-z0_i)/R;};

	virtual ~FIALKO();

	/**< Association of the medium param with the source buried in the Half-space
	  *< We also fix the assigned overpressure (DP_i).
	  */

	void SET(MEDIUM MEDIUM_PAR_i,double DP_i);

	/**< Computation of the displacement components (Ur, radial component and Uz, vertical component) at the observation point of polar coordinates (r,z) */

	void DISPLACEMENT(double r,double z,double &Ur,double &Uz);

	/**< Computation of the displacement components (Ux=U[0], Uy=U[1], Uz=U[2]) at the observation point (x,y,z) */

	void DISPLACEMENT(double x,double y,double z,double U[]);

	/**< Print to standard output useful information about the source */

	void PRINT(CONSOLE  &out);

	/**< Generation of a datafile and a gnuplot script.
	 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

	void PRINT(DATAFILE &out);

};

#endif /* FIALKO_H_ */
