//
// C++ class definition: BE3D.h
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef BE3D_H_
#define BE3D_H_


#include "BE3D_PARAM.h"
#include "foo.h"
#include "MEDIUM.h"



class BE3D
{

	double B1,B2,B3;		// Burgers vector components
	double BC1,BC2,BC3;		// Boundary conditions

	BE3D_PARAM PARAM;


	void Parametri(double Xi,double Eta,double QQ,int KXi,int KEta);

	void Par(double &SIN_DELTA,double &COS_DELTA,double &alpha,double &coef_alpha_1,double &coef_alpha_2,
			double delta_gradi,double mu,double lambda);

	void Singular_cases(double x,double y,double z,double c,double AL1,double AL2,double AW1,double AW2,
			double Vett_Xi[],double Vett_1_Eta[],double Vett_2_Eta[],int &flag_QM,int &flag_QP,
			int KXiM[],int KEtaM[],int KXiP[],int KEtaP[],double &QM,double &QP);

	void initialize(struct_var &struttura);


	void Strain_tensor(double Strain_components[],struct_var FF);
	void Stress_tensor(double S_Stress_components[],struct_var FF,double mu,double lambda);

	void IS_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void IS_B(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void IS_C(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void IS_u_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void IS_u_B(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void IS_u_C(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);

	void RS_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);
	void RS_u_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);

	void IS_A_Z0(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z);

	void Derivates(double burgers_vector[],struct_var FT[],double z);
	void CDisplacement(double burgers_vector[],struct_var FT[],double z);

	double fSX1A(double Xi,double Eta,double z);
	double fSX1B(double Xi,double Eta,double z);
	double fSX1C(double Xi,double Eta,double z);
	double fSX2A(double Xi,double Eta,double z);
	double fSX2B(double Xi,double Eta,double z);
	double fSX2C(double Xi,double Eta,double z);
	double fSX3A(double Xi,double Eta,double z);
	double fSX3B(double Xi,double Eta,double z);
	double fSX3C(double Xi,double Eta,double z);
	double fSY1A(double Xi,double Eta,double z);
	double fSY1B(double Xi,double Eta,double z);
	double fSY1C(double Xi,double Eta,double z);
	double fSY2A(double Xi,double Eta,double z);
	double fSY2B(double Xi,double Eta,double z);
	double fSY2C(double Xi,double Eta,double z);
	double fSY3A(double Xi,double Eta,double z);
	double fSY3B(double Xi,double Eta,double z);
	double fSY3C(double Xi,double Eta,double z);
	double fSZ1A(double Xi,double Eta,double z);
	double fSZ1B(double Xi,double Eta,double z);
	double fSu1C(double Xi,double Eta,double z);
	double fSZ1C(double Xi,double Eta,double z);
	double fSZ2A(double Xi,double Eta,double z);
	double fSZ2B(double Xi,double Eta,double z);
	double fSu2C(double Xi,double Eta,double z);
	double fSZ2C(double Xi,double Eta,double z);
	double fSZ3A(double Xi,double Eta,double z);
	double fSZ3B(double Xi,double Eta,double z);
	double fSu3C(double Xi,double Eta,double z);
	double fSZ3C(double Xi,double Eta,double z);

	double fSu1A(double Xi,double Eta,double z);
	double fSu1B(double Xi,double Eta,double z);
	double fSu2A(double Xi,double Eta,double z);
	double fSu3A(double Xi,double Eta,double z);
	double fSu2B(double Xi,double Eta,double z);
	double fSu3B(double Xi,double Eta,double z);


	double fDX1A(double Xi,double Eta,double z);
	double fDX1B(double Xi,double Eta,double z);
	double fDX1C(double Xi,double Eta,double z);
	double fDX2A(double Xi,double Eta,double z);
	double fDX2B(double Xi,double Eta,double z);
	double fDX2C(double Xi,double Eta,double z);
	double fDX3A(double Xi,double Eta,double z);
	double fDX3B(double Xi,double Eta,double z);
	double fDX3C(double Xi,double Eta,double z);
	double fDY1A(double Xi,double Eta,double z);
	double fDY1B(double Xi,double Eta,double z);
	double fDY1C(double Xi,double Eta,double z);
	double fDY2A(double Xi,double Eta,double z);
	double fDY2B(double Xi,double Eta,double z);
	double fDY2C(double Xi,double Eta,double z);
	double fDY3A(double Xi,double Eta,double z);
	double fDY3B(double Xi,double Eta,double z);
	double fDY3C(double Xi,double Eta,double z);
	double fDZ1A(double Xi,double Eta,double z);
	double fDZ1B(double Xi,double Eta,double z);
	double fDu1C(double Xi,double Eta,double z);
	double fDZ1C(double Xi,double Eta,double z);
	double fDZ2A(double Xi,double Eta,double z);
	double fDZ2B(double Xi,double Eta,double z);
	double fDu2C(double Xi,double Eta,double z);
	double fDZ2C(double Xi,double Eta,double z);
	double fDZ3A(double Xi,double Eta,double z);
	double fDZ3B(double Xi,double Eta,double z);
	double fDu3C(double Xi,double Eta,double z);
	double fDZ3C(double Xi,double Eta,double z);

	double fDu1A(double Xi,double Eta,double z);
	double fDu1B(double Xi,double Eta,double z);
	double fDu2A(double Xi,double Eta,double z);
	double fDu3A(double Xi,double Eta,double z);
	double fDu2B(double Xi,double Eta,double z);
	double fDu3B(double Xi,double Eta,double z);


	double fTX1A(double Xi,double Eta,double z);
	double fTX1B(double Xi,double Eta,double z);
	double fTX1C(double Xi,double Eta,double z);
	double fTX2A(double Xi,double Eta,double z);
	double fTX2B(double Xi,double Eta,double z);
	double fTX2C(double Xi,double Eta,double z);
	double fTX3A(double Xi,double Eta,double z);
	double fTX3B(double Xi,double Eta,double z);
	double fTX3C(double Xi,double Eta,double z);
	double fTY1A(double Xi,double Eta,double z);
	double fTY1B(double Xi,double Eta,double z);
	double fTY1C(double Xi,double Eta,double z);
	double fTY2A(double Xi,double Eta,double z);
	double fTY2B(double Xi,double Eta,double z);
	double fTY2C(double Xi,double Eta,double z);
	double fTY3A(double Xi,double Eta,double z);
	double fTY3B(double Xi,double Eta,double z);
	double fTY3C(double Xi,double Eta,double z);
	double fTZ1A(double Xi,double Eta,double z);
	double fTZ1B(double Xi,double Eta,double z);
	double fTu1C(double Xi,double Eta,double z);
	double fTZ1C(double Xi,double Eta,double z);
	double fTZ2A(double Xi,double Eta,double z);
	double fTZ2B(double Xi,double Eta,double z);
	double fTu2C(double Xi,double Eta,double z);
	double fTZ2C(double Xi,double Eta,double z);
	double fTZ3A(double Xi,double Eta,double z);
	double fTZ3B(double Xi,double Eta,double z);
	double fTu3C(double Xi,double Eta,double z);
	double fTZ3C(double Xi,double Eta,double z);

	double fTu1A(double Xi,double Eta,double z);
	double fTu1B(double Xi,double Eta,double z);
	double fTu2A(double Xi,double Eta,double z);
	double fTu3A(double Xi,double Eta,double z);
	double fTu2B(double Xi,double Eta,double z);
	double fTu3B(double Xi,double Eta,double z);


	void Displacement_Okada(double S_Displacement[],double D_Displacement[],double T_Displacement[],
			double burgers_vector[],double x,double y,double z,double c,double delta_gradi,
			double AL1,double AL2,double AW1,double AW2,
			double mu,double lambda);

	void Strain_Okada(double Strain[],double burgers_vector[],double x,double y,double z,double c,
			double delta_gradi,double AL1,double AL2,double AW1,double AW2,
			double mu,double lambda);

	void Stress_Okada(double S_Stress_components[],double D_Stress_components[],double T_Stress_components[],
			double burgers_vector[],double x,double y,double z,double c,
			double delta_gradi,double AL1,double AL2,double AW1,double AW2,
			double mu,double lambda);

public:

	BE3D(){};

	BE3D(double xc,double yc,double zc,double xp,double yp,double zp,
			double BE_phi,double BE_delta,
			double BE_L,double BE_W,
			double n_x,double n_y,double n_z,
	  	  	double ts_x,double ts_y,double ts_z,
	  	  	double td_x,double td_y,double td_z)
	//			  int BE_flag,double nx,double ny,double nz)
	{
		B1=0; B2=0; B3=0;
		BC1=0; BC2=0; BC3=0;

		BE3D_GEOM_PARAM(xc,yc,zc,xp,yp,zp,BE_phi,BE_delta,BE_L,BE_W,n_x,n_y,n_z,ts_x,ts_y,ts_z,td_x,td_y,td_z);
	};

	virtual ~BE3D(){};


//	void NEW(double xc,double yc,double zc,double xp,double yp,double zp,
//				double BE_phi,double BE_delta,double BE_L,double BE_W,double nx,double ny,double nz);

	void NEW(double xc,double yc,double zc,
				double xp,double yp,double zp,
			  	double BE_phi,double BE_delta,
			  	double BE_L,double BE_W,
			  	double n_x,double n_y,double n_z,
			  	double ts_x,double ts_y,double ts_z,
			  	double td_x,double td_y,double td_z);


	void DISPLACEMENT(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,
						double S_Displacement[],double D_Displacement[],double T_Displacement[]);

	void DISPLACEMENT(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Displacement[]);


	void STRESS(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,
					double S_Stress[],double D_Stress[],double T_Stress[]);

	void STRESS(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Stress[]);


	void STRAIN(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Strain[]);


	void GET_PARAM(BE3D_PARAM &PARAM_o){PARAM_o = PARAM;};

	void GET_BV(double &B1_o,double &B2_o,double &B3_o){B1_o = B1; B2_o = B2; B3_o = B3;};
	void PUT_BV(double B1_i,double B2_i,double B3_i){B1=B1_i;B2=B2_i;B3=B3_i;};

	void GET_BC(double &BC1_o,double &BC2_o,double &BC3_o){BC1_o = BC1; BC2_o = BC2; BC3_o = BC3;};
	void PUT_BC(double BC1_i,double BC2_i,double BC3_i){BC1=BC1_i;BC2=BC2_i;BC3=BC3_i;};

	void SHIFT(double Dx,double Dy,double Dz);

};



#endif  /* BE3D_H_ */

