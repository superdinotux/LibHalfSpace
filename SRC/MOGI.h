/*
 * MOGI.h
 *
 *  Created on: 23/ott/2014
 *      Author: dino
 */

#ifndef MOGI_H_
#define MOGI_H_


#include "PRESSURIZED.h"


class MOGI: public PRESSURIZED
{

	double R;

    void PRINT_GEOM_PROP(ostream *out_stream);
    void PRINT_STRESS_PROP(ostream *stream);

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print=0);

public:

	MOGI();

	MOGI(string label_i,double x0_i,double y0_i,double z0_i,double R_i):
			PRESSURIZED(label_i,x0_i,y0_i,z0_i)
			{R = R_i;};

	virtual ~MOGI();


	void SET(MEDIUM MEDIUM_PAR_i,double P_i);

	void DISPLACEMENT(double r,double z,double &Ur,double &Uz);
	void DISPLACEMENT(double x,double y,double z,double U[]);

	void PRINT(CONSOLE  &out);
	void PRINT(DATAFILE &out);

};

#endif /* MOGI_H_ */
