/*
 * MEDIUM.h
 *
 *  Created on: 07/lug/2014
 *      Author: dino
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

class MEDIUM {

//! Elastic properties of the medium

	double mu;
	double nu;
	double lambda;
	double K;
	double E;

//! Medium properties in terms of density

	int num_density_layers;
	double depth_density_trans[];   // [m] 		array of num_density_layer elements filled with depth density transitions
	double layer_density[];      	// [kg/m^3] array of num_density_layer+1 elements filled with layer densities

//! Medium properties in terms of tectonic stress

	double stress_depth;  // [m]  depth under which the tectonic stress does not act any more
	double Sxx;      	  // [Pa] XX stress component
	double Syy;      	  // [Pa] YY stress component
	double Szz;      	  // [Pa] ZZ stress component
	double Sxy;      	  // [Pa] XY stress component
	double Sxz;           // [Pa] XZ stress component
	double Syz;      	  // [Pa] YZ stress component

//! Medium properties in terms of friction

	double f_coef;        // coefficient of friction

public:

	MEDIUM(){mu=0; nu=0; lambda=0; K=0; E=0;};

	MEDIUM(double mu,double nu);

	virtual ~MEDIUM(){};

	double GET_mu();
	double GET_nu();
	double GET_lambda();
	double GET_K();
	double GET_E();

};

#endif /* MEDIUM_H_ */
