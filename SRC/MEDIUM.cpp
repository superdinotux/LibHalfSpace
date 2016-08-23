/*
 * MEDIUM.cpp
 *
 *  Created on: 07/lug/2014
 *      Author: dino
 */


#include "MEDIUM.h"


MEDIUM::MEDIUM(double mu_i,double nu_i){
	// TODO Auto-generated constructor stub

	mu=mu_i;

	nu=nu_i;

	lambda=2.e0*mu*nu/(1.e0-2.e0*nu);

	K=2.e0*mu*(1.e0+nu)/(3.e0*(1.e0-2.e0*nu));

	E=2.e0*mu*(1.e0+nu);

}



double MEDIUM::GET_mu()
{
	double mu_o=mu;
	return mu_o;
}

double MEDIUM::GET_nu()
{
	double nu_o=nu;
	return nu_o;
}

double MEDIUM::GET_lambda()
{
	double lambda_o=lambda;
	return lambda_o;
}

double MEDIUM::GET_K()
{
	double K_o=K;
	return K_o;
}

double MEDIUM::GET_E()
{
	double E_o=E;
	return E_o;
}











