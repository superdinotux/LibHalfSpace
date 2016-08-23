/*
 * NOT_PLANAR.h
 *
 *  Created on: 29/gen/2015
 *      Author: dino
 */



#ifndef NOT_PLANAR_H_
#define NOT_PLANAR_H_


#include "HalfSpace.h"


class NOT_PLANAR: public SOURCE_BE3D {

private: 			// All PRIVATE MEMBERS of the CLASS are DATA

	double L;
	double W_1,W_2;
	double teta_1,teta_2;

public: 			// All PUBLIC MEMBERS of the CLASS are METHODS

	NOT_PLANAR();

	virtual ~NOT_PLANAR();

	NOT_PLANAR(string label,double x0,double y0,double z0,
				double L_i,double W1,double W2,double phi,double teta1,double teta2,int NDP);

};


#endif /* NOT_PLANAR_H_ */



