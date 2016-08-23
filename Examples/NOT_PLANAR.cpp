/*
 * NOT_PLANAR.cpp
 *
 *  Created on: 29/gen/2015
 *      Author: dino
 */

#include "NOT_PLANAR.h"


struct grid{
	int Nx,Ny;
	double step_x,step_y;
	double cc[10000][2],pc[10000][2];
};


void retrieve_grid_param(double L_i,double W_i,grid &grid_o,int NDP_i);


NOT_PLANAR::NOT_PLANAR() {
	// TODO Auto-generated constructor stub

}

NOT_PLANAR::~NOT_PLANAR() {
	// TODO Auto-generated destructor stub
}


NOT_PLANAR::NOT_PLANAR(string label,double x0,double y0,double z0,
		double L_i,double W1,double W2,double phi,double teta1,double teta2,int NDP):SOURCE_BE3D(label,x0,y0,z0,phi,NDP)
{

	L=L_i;
	W_1=W1;
	W_2=W2;
	teta_1=teta1;
	teta_2=teta2;

	double teta_1_rad=(teta_1/180e0)*acos(-1e0),teta_2_rad=(teta_2/180e0)*acos(-1e0);

	double y01=y0-W_1*cos(teta_1_rad)*0.5,y02=y0+W_2*cos(teta_2_rad)*0.5;
	double z01=z0-W_1*sin(teta_1_rad)*0.5,z02=z0+W_2*sin(teta_2_rad)*0.5;

	grid grid_1,grid_2;

	int NDP_1=NDP;
	int NDP_2=NDP;

	double area_1,area_2;

	retrieve_grid_param(L,W_1,grid_1,NDP_1);
	retrieve_grid_param(L,W_2,grid_2,NDP_2);

	area_1=grid_1.step_x*grid_1.step_y;
	area_2=grid_2.step_x*grid_2.step_y;

	if(area_2 > area_1)
	{
		while((area_2-area_1)/area_1 > 0.01)
		{
			NDP_2++;
			retrieve_grid_param(L,W_2,grid_2,NDP_2);
			area_2=grid_2.step_x*grid_2.step_y;
		}
	}
	else if(area_1 > area_2)
	{
		while((area_1-area_2)/area_2 > 0.01)
		{
			NDP_1++;
			retrieve_grid_param(L,W_1,grid_1,NDP_1);
			area_1=grid_1.step_x*grid_1.step_y;
		}
	}

	int NBE_1=grid_1.Nx*grid_1.Ny;
	int NBE_2=grid_2.Nx*grid_2.Ny;

	NBE=NBE_1+NBE_2;

	BE = new BE3D [NBE];

	double xc,yc,zc,xp,yp,zp,angle,n_x,n_y,n_z,td_x,td_y,td_z,BE_L,BE_W;

    double ts_x =  1;
    double ts_y =  0;
    double ts_z =  0;


	for(int tj=0; tj < NBE; tj++)
	{
		int tk=tj-NBE_1;

		if(tj < NBE_1)
		{
			xc = grid_1.cc[tj][0];
			yc = grid_1.cc[tj][1]*cos(teta_1_rad) + y01;
			zc = grid_1.cc[tj][1]*sin(teta_1_rad) + z01;
			xp = grid_1.pc[tj][0];
			yp = grid_1.pc[tj][1]*cos(teta_1_rad) + y01;
			zp = grid_1.pc[tj][1]*sin(teta_1_rad) + z01;

			angle = teta_1;

			n_x = 0e0;
			n_y = -sin(teta_1_rad);
			n_z =  cos(teta_1_rad);

		    td_x = 0e0;
		    td_y = cos(teta_1_rad);
		    td_z = sin(teta_1_rad);

			BE_L = grid_1.step_x;
			BE_W = grid_1.step_y;
		}
		else
		{
			xc = grid_2.cc[tk][0];
			yc = grid_2.cc[tk][1]*cos(teta_2_rad) + y02;
			zc = grid_2.cc[tk][1]*sin(teta_2_rad) + z02;
			xp = grid_2.pc[tk][0];
			yp = grid_2.pc[tk][1]*cos(teta_2_rad) + y02;
			zp = grid_2.pc[tk][1]*sin(teta_2_rad) + z02;

			angle = teta_2;

			n_x = 0e0;
			n_y = -sin(teta_2_rad);
			n_z =  cos(teta_2_rad);

		    td_x = 0e0;
		    td_y = cos(teta_2_rad);
		    td_z = sin(teta_2_rad);

			BE_L = grid_2.step_x;
			BE_W = grid_2.step_y;
		}

		BE[tj].NEW(xc,yc,zc,xp,yp,zp,0e0,angle,BE_L,BE_W,n_x,n_y,n_z,ts_x,ts_y,ts_z,td_x,td_y,td_z);

	}

}



//void NOT_PLANAR::retrieve_grid_param(double W_i,grid &grid_o,int NDP_i)
void retrieve_grid_param(double L,double W_i,grid &grid_o,int NDP_i)
{
	int NUM=NDP_i;

	int N_1;
	int N_2;

	double LATO_h_1=L;
	double LATO_h_2=W_i;

	double step_1;
	double step_2;

	double cfr[NUM];
	double min;


	if((L == W_i) || (NUM == 1))
	{
		N_1=NUM;
		N_2=NUM;

		step_1=LATO_h_1/N_1;
		step_2=LATO_h_2/N_2;
	}
	else if(L > W_i)
	{
		N_1=NUM;

		step_1=LATO_h_1/N_1;

		for(int tk=0; tk < N_1; tk++)
		{
			step_2=LATO_h_2/(tk+1);
			cfr[tk]=abs((step_2-step_1)/step_1);
		}

		min=cfr[0];

		for(int tk=0; tk < N_1; tk++)
		{

			if(cfr[tk] <= min)
			{
				min=cfr[tk];
				N_2=tk+1;
			}

		}

		step_2=LATO_h_2/N_2;
	}
	else
	{
		N_2=NUM;
		step_2=LATO_h_2/N_2;

		for(int tk=0; tk < N_2; tk++)
		{
			step_1=LATO_h_1/(tk+1);
			cfr[tk]=abs((step_1-step_2)/step_2);
		}


		min=cfr[0];

		for(int tk=0; tk < N_2; tk++)
		{
			if(cfr[tk] <= min)
			{
				min=cfr[tk];
				N_1=tk+1;
			}
		}

		step_1=LATO_h_1/N_1;
	}

	grid_o.Nx=N_1;
	grid_o.Ny=N_2;

	grid_o.step_x=step_1;
	grid_o.step_y=step_2;


	double xccI=-L*0.5;
	double yccI=-W_i*0.5;

	int N_x=grid_o.Nx;
	int N_y=grid_o.Ny;

	double step_x=grid_o.step_x;
	double step_y=grid_o.step_y;

	double xpcI=-(N_x-1)*step_x*0.5;
	double ypcI=-(N_y-1)*step_y*0.5;

	double xcc,xpc;
	double ycc,ypc;

	int tj;

	for(int tx=0; tx < N_x; tx++)
	{

		xcc=xccI+step_x*tx;
		xpc=xpcI+step_x*tx;

		for(int ty=0; ty < N_y; ty++)
		{

			ycc=yccI+step_y*ty;
			ypc=ypcI+step_y*ty;

			tj=tx*N_y+ty;

			grid_o.cc[tj][0]=xcc;
			grid_o.cc[tj][1]=ycc;

			grid_o.pc[tj][0]=xpc;
			grid_o.pc[tj][1]=ypc;
		}
	}

}





