/*
 * BE3D_GEOM_PARAM.h
 *
 *  Created on: 07/lug/2014
 *      Author: dino
 */

#ifndef BE3D_GEOM_PARAM_H_
#define BE3D_GEOM_PARAM_H_


#include <cmath>



class BE3D_GEOM_PARAM {

public:

    double cc[3];
    double pc[3];

    double n[3];
    double ts[3];
    double td[3];

    double phi_degree;
    double delta_degree;

    double L;
    double W;

    int flag;

	BE3D_GEOM_PARAM(){};
//	virtual ~BE3D_GEOM_PARAM();


	BE3D_GEOM_PARAM(double xc,double yc,double zc,
					double xp,double yp,double zp,
					double BE_phi,double BE_delta,
					double BE_L,double BE_W,
					double n_x,double n_y,double n_z,
					double ts_x,double ts_y,double ts_z,
					double td_x,double td_y,double td_z)
	{
		cc[0] = xc;
		cc[1] = yc;
		cc[2] = zc;

		pc[0] = xp;
		pc[1] = yp;
		pc[2] = zp;

		phi_degree = BE_phi;
		delta_degree = BE_delta;

		double pi=acos(-1);

		double phi_rad = pi*phi_degree/180e0;
		double delta_rad = pi*delta_degree/180e0;

		L = BE_L;
		W = BE_W;

		n[0]=n_x;
		n[1]=n_y;
		n[2]=n_z;

		ts[0]=ts_x;
		ts[1]=ts_y;
		ts[2]=ts_z;

		td[0]=td_x;
		td[1]=td_y;
		td[2]=td_z;
	};


};






#endif /* BE3D_GEOM_PARAM_H_ */




