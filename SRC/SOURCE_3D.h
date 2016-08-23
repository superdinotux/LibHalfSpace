/*
 * SOURCE_3D.h
 *
 *  Created on: 18/nov/2014
 *      Author: dino
 */

#ifndef SOURCE_3D_H_
#define SOURCE_3D_H_

#include "SOURCE.h"

class SOURCE_3D: public SOURCE {

public:

	SOURCE_3D();

	SOURCE_3D(string s_label,double x0,double y0,double z0):SOURCE(s_label,x0,y0,z0){};

	virtual ~SOURCE_3D();

};

#endif /* SOURCE_3D_H_ */
