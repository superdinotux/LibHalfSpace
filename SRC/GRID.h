/*
 * GRID.h
 *
 *  Created on: 02/lug/2014
 *      Author: dino
 */

#ifndef GRID_H_
#define GRID_H_

class GRID
{
public:

	 int N_1;
	 int N_2;
	 int N_3;

	 int NX;
	 int NY;
	 int NZ;

	 int TOTAL;

	 int NUM_COND;

	 GRID(){};
	 GRID(int N_1_i,int N_2_i,int N_3_i,int N_COND = 0)
	 	 {N_1=N_1_i; N_2=N_2_i; N_3=N_3_i; NUM_COND=N_COND; NX=N_2*N_3; NY=N_1*N_3; NZ=N_1*N_2; TOTAL=(NX+NY+NZ)*2;};
	 ~GRID(){};
};

#endif /* GRID_H_ */
