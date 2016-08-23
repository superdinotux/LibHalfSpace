/*
 * BE3D_STRESS_PARAM.h
 *
 *  Created on: 09/lug/2014
 *      Author: dino
 */

#ifndef BE3D_STRESS_PARAM_H_
#define BE3D_STRESS_PARAM_H_

class BE3D_STRESS_PARAM{

	double Delta_P;

public:

	BE3D_STRESS_PARAM(){};
	BE3D_STRESS_PARAM(double DP){Delta_P = DP;};

//	virtual ~BE3D_STRESS_PARAM(){};

};

#endif /* BE3D_STRESS_PARAM_H_ */
