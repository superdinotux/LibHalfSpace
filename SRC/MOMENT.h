/*
 * MOMENT.h
 *
 *  Created on: 24/ott/2014
 *      Author: dino
 */

#ifndef MOMENT_H_
#define MOMENT_H_


#include "SOURCE_3D.h"


class MOMENT: public SOURCE_3D
{
private:

	void PRINT_GEOM_PROP(ostream *out_stream);
	void PRINT_STRESS_PROP(ostream *out_stream);

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print = 0);

protected:

	double M[6];		// M[0]=MXX,M[1]=MYY,M[2]=MZZ,M[3]=MXY,M[4]=MXZ,M[5]=MYZ

public:

	MOMENT();

	MOMENT(string label,double x0_i,double y0_i,double z0_i):
				SOURCE_3D(label,x0_i,y0_i,z0_i){};

	virtual ~MOMENT();

	/*! The medium parameters are assigned to the source. We also assigns the moment tensor associated to the source.  */

	void SET(MEDIUM medium_par_i,double MXX,double MYY,double MZZ,double MXY,double MXZ,double MYZ)
	{
		MEDIUM_PAR = medium_par_i;

		M[0]=MXX;
		M[1]=MYY;
		M[2]=MZZ;
		M[3]=MXY;
		M[4]=MXZ;
		M[5]=MYZ;
	}

	void DISPLACEMENT(double x,double y,double z,double U[]);

	void PRINT(CONSOLE &out);
	void PRINT(DATAFILE &out);

};






#endif /* MOMENT_H_ */
