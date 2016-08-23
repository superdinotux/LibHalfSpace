/*
 * PRESSURIZED.h
 *
 *  Created on: 25/nov/2014
 *      Author: dino
 */

#ifndef PRESSURIZED_H_
#define PRESSURIZED_H_

#include "SOURCE_3D.h"

class PRESSURIZED: public SOURCE_3D {

protected:

	double DP;

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print=0);

public:

	PRESSURIZED();

	PRESSURIZED(string s_label,double x0,double y0,double z0):SOURCE_3D(s_label,x0,y0,z0){};

	virtual ~PRESSURIZED();

	virtual void SET(MEDIUM MEDIUM_PAR_i,double P_i) =0;

};

#endif /* PRESSURIZED_H_ */
