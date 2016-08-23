/*
 * BE3DPARAM.h
 *
 *  Created on: 09/lug/2014
 *      Author: dino
 */

#ifndef BE3DPARAM_H_
#define BE3DPARAM_H_

#include "BE3D_GEOM_PARAM.h"
#include "BE3D_STRESS_PARAM.h"

class BE3D_PARAM {

public:

	BE3D_GEOM_PARAM GEOM;
	BE3D_STRESS_PARAM STRESS;

	BE3D_PARAM()
	{
		BE3D_GEOM_PARAM GEOM;
		BE3D_STRESS_PARAM STRESS;
	};

	virtual ~BE3D_PARAM(){};

};

#endif /* BE3DPARAM_H_ */
