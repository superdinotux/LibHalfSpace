//
// C++ Implementation: Okada
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <cmath>
#include <iostream>


using namespace std;


#include "BE3D.h"
#include "MEDIUM.h"
#include "Utilities.h"


const double EPS=1e-6;

double SIN_DELTA,COS_DELTA,SIN_DELTA_2,COS_DELTA_2;
double alpha;
double coef_alpha_1,coef_alpha_2;

double p,q;
double Y,D;
double RpD,RpD2;
double Xi2,Eta2,q2,q3,Y2,D2;
double R,R2,R3,R5;

double D_11;
double X_11,X_32,X_53;
double Y_11,Y_32,Y_53;
double Y_0;

double E,F,G,H;
double Ep,Fp,Gp,Hp;

double Teta,X;

double Sing_iii,Sing_iv;

double I_1,I_2,I_3,I_4; 
double J_1,J_2,J_3,J_4,J_5,J_6;  
double K_1,K_2,K_3,K_4;

double C,h;
double Z_32,Z_53,Z_0;
double CR;
double P,Q,Pp,Qp;
double qR,CDR,YY0;
	
const float pi=acos(-1);
const float coefpi=1.E0/(2.E0*pi);





void DC3D(float alpha,float X,float Y,float Z,float DEPTH,float DIP,float AL1,float AL2,float AW1,float AW2,
			float DISL1,float DISL2,float DISL3,float &UX,float &UY,float &UZ,
			float &UXX,float &UYX,float &UZX,float &UXY,float &UYY,float &UZY,float &UXZ,float &UYZ,float &UZZ,int &IRET);



//void BE3D::NEW(double xc,double yc,double zc,double xp,double yp,double zp,
//		  double BE_phi,double BE_delta,
//		  double BE_L,double BE_W,
//		  double nx,double ny,double nz)
//{
//	B1=0; B2=0; B3=0;
//	BC1=0; BC2=0; BC3=0;
//
//	PARAM.GEOM.cc[0] = xc;
//	PARAM.GEOM.cc[1] = yc;
//	PARAM.GEOM.cc[2] = zc;
//
//	PARAM.GEOM.pc[0] = xp;
//	PARAM.GEOM.pc[1] = yp;
//	PARAM.GEOM.pc[2] = zp;
//
//	PARAM.GEOM.phi_degree = BE_phi;
//	PARAM.GEOM.delta_degree = BE_delta;
//
//	PARAM.GEOM.L = BE_L;
//	PARAM.GEOM.W = BE_W;
//
//	PARAM.GEOM.n[0]=nx;
//	PARAM.GEOM.n[1]=ny;
//	PARAM.GEOM.n[2]=nz;
//};



void BE3D::NEW(double xc,double yc,double zc,
				  double xp,double yp,double zp,
		  	  	  double BE_phi,double BE_delta,
		  	  	  double BE_L,double BE_W,
		  	  	  double n_x,double n_y,double n_z,
		  	  	  double ts_x,double ts_y,double ts_z,
		  	  	  double td_x,double td_y,double td_z)
{
	B1=0; B2=0; B3=0;
	BC1=0; BC2=0; BC3=0;

	PARAM.GEOM.cc[0] = xc;
	PARAM.GEOM.cc[1] = yc;
	PARAM.GEOM.cc[2] = zc;

	PARAM.GEOM.pc[0] = xp;
	PARAM.GEOM.pc[1] = yp;
	PARAM.GEOM.pc[2] = zp;

	PARAM.GEOM.phi_degree = BE_phi;
	PARAM.GEOM.delta_degree = BE_delta;

	PARAM.GEOM.L = BE_L;
	PARAM.GEOM.W = BE_W;

	PARAM.GEOM.n[0]=n_x;
	PARAM.GEOM.n[1]=n_y;
	PARAM.GEOM.n[2]=n_z;

	PARAM.GEOM.ts[0]=ts_x;
	PARAM.GEOM.ts[1]=ts_y;
	PARAM.GEOM.ts[2]=ts_z;

	PARAM.GEOM.td[0]=td_x;
	PARAM.GEOM.td[1]=td_y;
	PARAM.GEOM.td[2]=td_z;
};






void BE3D::STRAIN(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Strain[])
{

	double AL1 = -PARAM.GEOM.L*0.5;
	double AL2 =  PARAM.GEOM.L*0.5;
	double AW1 = -PARAM.GEOM.W*0.5;
	double AW2 =  PARAM.GEOM.W*0.5;

	double x_0 =  PARAM.GEOM.pc[0];
	double y_0 =  PARAM.GEOM.pc[1];
	double   c = -PARAM.GEOM.pc[2];


    if(FLAG_COMP == 1)
	{

	    double mu=M_PAR.GET_mu();
	    double lambda=M_PAR.GET_lambda();

	    double delta_gradi=PARAM.GEOM.delta_degree;

	    double x = xn-x_0;
	    double y = yn-y_0;

	    Rotazione_coord(x,y,PARAM.GEOM.phi_degree);

	    double burgers_vector[3];

	    burgers_vector[0] = B1;
	    burgers_vector[1] = B2;
	    burgers_vector[2] = B3;

	    Strain_Okada(Strain,burgers_vector,x,y,z,c,delta_gradi,AL1,AL2,AW1,AW2,mu,lambda);

	}
	else if(FLAG_COMP == 2)
	{

	    float f_AL1=(float) AL1;
	    float f_AL2=(float) AL2;
	    float f_AW1=(float) AW1;
	    float f_AW2=(float) AW2;

	    float f_mu=(float) M_PAR.GET_mu();
	    float f_lambda=(float) M_PAR.GET_lambda();

	    float f_z=(float) z;

	    float f_delta_gradi=(float) PARAM.GEOM.delta_degree;

	    float f_x = (float) xn-x_0;
	    float f_y = (float) yn-y_0;
	    float f_c = (float) c;

	    Rotazione_coord(f_x,f_y,PARAM.GEOM.phi_degree);

	    float UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ;
	    int IRET;
	    float Slip_x,Slip_y,Slip_z;

	    float alpha=(f_lambda+f_mu)/(f_lambda+2*f_mu);


	    float DISL1a=(float) B1;
	    float DISL2a=(float) B2;
	    float DISL3a=(float) B3;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1a,DISL2a,DISL3a,
	    		Slip_x,Slip_y,Slip_z,
	    		UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    Strain[0] = (double) UXX;
	    Strain[1] = (double) UYY;
	    Strain[2] = (double) UZZ;
	    Strain[3] = (double) 0.5*(UXY+UYX);
	    Strain[4] = (double) 0.5*(UXZ+UZX);
	    Strain[5] = (double) 0.5*(UYZ+UZY);

	}

	Rotazione_tensore(Strain,PARAM.GEOM.phi_degree);

}



void BE3D::STRESS(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,
					double S_Stress[],double D_Stress[],double T_Stress[])
{

	double AL1 = -PARAM.GEOM.L*0.5;
	double AL2 =  PARAM.GEOM.L*0.5;
	double AW1 = -PARAM.GEOM.W*0.5;
	double AW2 =  PARAM.GEOM.W*0.5;

	double x_0 =  PARAM.GEOM.pc[0];
	double y_0 =  PARAM.GEOM.pc[1];
	double   c = -PARAM.GEOM.pc[2];


	if(FLAG_COMP == 1)
	{

	    double mu=M_PAR.GET_mu();
	    double lambda=M_PAR.GET_lambda();

	    double delta_gradi=PARAM.GEOM.delta_degree;

	    double x = xn-x_0;
	    double y = yn-y_0;

	    Rotazione_coord(x,y,PARAM.GEOM.phi_degree);


	    double burgers_vector[3];

	    burgers_vector[0] = B1;
	    burgers_vector[1] = B2;
	    burgers_vector[2] = B3;

	    Stress_Okada(S_Stress,D_Stress,T_Stress,burgers_vector,x,y,z,c,delta_gradi,AL1,AL2,AW1,AW2,mu,lambda);
	}
	else if(FLAG_COMP == 2)
	{

	    float f_AL1=(float) AL1;
	    float f_AL2=(float) AL2;
	    float f_AW1=(float) AW1;
	    float f_AW2=(float) AW2;

	    float f_mu=(float) M_PAR.GET_mu();
	    float f_lambda=(float) M_PAR.GET_lambda();

	    float f_z=(float) z;

	    float f_delta_gradi=(float) PARAM.GEOM.delta_degree;

	    float f_x = (float) xn-x_0;
	    float f_y = (float) yn-y_0;
	    float f_c = (float) c;

	    Rotazione_coord(f_x,f_y,PARAM.GEOM.phi_degree);


	    float UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ;
	    float e_kk;
	    int IRET;
	    float Slip_x,Slip_y,Slip_z;

	    float alpha=(f_lambda+f_mu)/(f_lambda+2*f_mu);


	    float DISL1a=1;
	    float DISL2a=0;
	    float DISL3a=0;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1a,DISL2a,DISL3a,
	    		Slip_x,Slip_y,Slip_z,
	    		UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    e_kk=UXX+UYY+UZZ;

	    S_Stress[0] = (double) f_lambda*e_kk+2*f_mu*UXX;
	    S_Stress[1] = (double) f_lambda*e_kk+2*f_mu*UYY;
	    S_Stress[2] = (double) f_lambda*e_kk+2*f_mu*UZZ;
	    S_Stress[3] = (double) f_mu*(UXY+UYX);
	    S_Stress[4] = (double) f_mu*(UXZ+UZX);
	    S_Stress[5] = (double) f_mu*(UYZ+UZY);


	    float DISL1b=0;
	    float DISL2b=1;
	    float DISL3b=0;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1b,DISL2b,DISL3b,
	    		Slip_x,Slip_y,Slip_z,
	    		UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    e_kk=UXX+UYY+UZZ;

	    D_Stress[0] = (double) f_lambda*e_kk+2*f_mu*UXX;
	    D_Stress[1] = (double) f_lambda*e_kk+2*f_mu*UYY;
	    D_Stress[2] = (double) f_lambda*e_kk+2*f_mu*UZZ;
	    D_Stress[3] = (double) f_mu*(UXY+UYX);
	    D_Stress[4] = (double) f_mu*(UXZ+UZX);
	    D_Stress[5] = (double) f_mu*(UYZ+UZY);


	    float DISL1c=0;
	    float DISL2c=0;
	    float DISL3c=1;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1c,DISL2c,DISL3c,
	    		Slip_x,Slip_y,Slip_z,
	    		UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    e_kk=UXX+UYY+UZZ;

	    T_Stress[0]=(double) f_lambda*e_kk+2*f_mu*UXX;
	    T_Stress[1]=(double) f_lambda*e_kk+2*f_mu*UYY;
	    T_Stress[2]=(double) f_lambda*e_kk+2*f_mu*UZZ;
	    T_Stress[3]=(double) f_mu*(UXY+UYX);
	    T_Stress[4]=(double) f_mu*(UXZ+UZX);
	    T_Stress[5]=(double) f_mu*(UYZ+UZY);
	}


	Rotazione_tensore(S_Stress,PARAM.GEOM.phi_degree);
	Rotazione_tensore(D_Stress,PARAM.GEOM.phi_degree);
	Rotazione_tensore(T_Stress,PARAM.GEOM.phi_degree);

}



void BE3D::STRESS(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Stress[])
{

	double AL1 = -PARAM.GEOM.L*0.5;
	double AL2 =  PARAM.GEOM.L*0.5;
	double AW1 = -PARAM.GEOM.W*0.5;
	double AW2 =  PARAM.GEOM.W*0.5;

	double x_0 =  PARAM.GEOM.pc[0];
	double y_0 =  PARAM.GEOM.pc[1];
	double   c = -PARAM.GEOM.pc[2];

    double S_Stress[6],D_Stress[6],T_Stress[6];

	if(FLAG_COMP == 1)
	{

	    double mu=M_PAR.GET_mu();
	    double lambda=M_PAR.GET_lambda();

	    double delta_gradi=PARAM.GEOM.delta_degree;

	    double x = xn-x_0;
	    double y = yn-y_0;

	    Rotazione_coord(x,y,PARAM.GEOM.phi_degree);

	    double burgers_vector[3];

	    burgers_vector[0] = B1;
	    burgers_vector[1] = B2;
	    burgers_vector[2] = B3;

	    Stress_Okada(S_Stress,D_Stress,T_Stress,burgers_vector,x,y,z,c,delta_gradi,AL1,AL2,AW1,AW2,mu,lambda);

		Stress[0]=B1*S_Stress[0]+B2*D_Stress[0]+B3*T_Stress[0];
		Stress[1]=B1*S_Stress[1]+B2*D_Stress[1]+B3*T_Stress[1];
		Stress[2]=B1*S_Stress[2]+B2*D_Stress[2]+B3*T_Stress[2];
		Stress[3]=B1*S_Stress[3]+B2*D_Stress[3]+B3*T_Stress[3];
		Stress[4]=B1*S_Stress[4]+B2*D_Stress[4]+B3*T_Stress[4];
		Stress[5]=B1*S_Stress[5]+B2*D_Stress[5]+B3*T_Stress[5];

	}
	else if(FLAG_COMP == 2)
	{

	    float f_AL1=(float) AL1;
	    float f_AL2=(float) AL2;
	    float f_AW1=(float) AW1;
	    float f_AW2=(float) AW2;

	    float f_mu=(float) M_PAR.GET_mu();
	    float f_lambda=(float) M_PAR.GET_lambda();

	    float f_z=(float) z;

	    float f_delta_gradi=(float) PARAM.GEOM.delta_degree;

	    float f_x = (float) xn-x_0;
	    float f_y = (float) yn-y_0;
	    float f_c = (float) c;

	    Rotazione_coord(f_x,f_y,PARAM.GEOM.phi_degree);

	    float UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ;
	    float e_kk;
	    int IRET;
	    float Slip_x,Slip_y,Slip_z;

	    float alpha=(f_lambda+f_mu)/(f_lambda+2*f_mu);


	    float DISL1a=(float) B1;
	    float DISL2a=(float) B2;
	    float DISL3a=(float) B3;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1a,DISL2a,DISL3a,
	    		Slip_x,Slip_y,Slip_z,
	    		UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    e_kk=UXX+UYY+UZZ;

	    Stress[0] = (double) f_lambda*e_kk+2*f_mu*UXX;
	    Stress[1] = (double) f_lambda*e_kk+2*f_mu*UYY;
	    Stress[2] = (double) f_lambda*e_kk+2*f_mu*UZZ;
	    Stress[3] = (double) f_mu*(UXY+UYX);
	    Stress[4] = (double) f_mu*(UXZ+UZX);
	    Stress[5] = (double) f_mu*(UYZ+UZY);

	}

	Rotazione_tensore(Stress,PARAM.GEOM.phi_degree);

}



// # Componenti della deformazione generate da una dislocazione rettangolare (Okada 1992)

void BE3D::Strain_Okada(double Strain[],double burgers_vector[],double x,double y,double z,double c,
							double delta_gradi,double AL1,double AL2,double AW1,double AW2,double mu,double lambda)
{

	Par(SIN_DELTA,COS_DELTA,alpha,coef_alpha_1,coef_alpha_2,delta_gradi,mu,lambda);

	struct_var FT[3];

	for(int ti=0; ti < 3; ti++)
	{
	    initialize(FT[ti]);
	}

	double Xi,Eta;
	double Vett_Xi[2];
	double Vett_1_Eta[2],Vett_2_Eta[2];
	int KXiM[2],KEtaM[2];
	int KXiP[2],KEtaP[2];
	double QM=0,QP=0;
	int flag_QM=0,flag_QP=0;

	Singular_cases(x,y,z,c,AL1,AL2,AW1,AW2,Vett_Xi,Vett_1_Eta,Vett_2_Eta,flag_QM,flag_QP,KXiM,KEtaM,KXiP,KEtaP,QM,QP);


	if((flag_QM == 0) && (flag_QP == 0))
	{

	    double sign;

	    for(int t_Xi=0; t_Xi < 2; t_Xi++)
	    {

		Xi=Vett_Xi[t_Xi];


		for(int t_Eta=0; t_Eta < 2; t_Eta++)
		{

		    sign=pow(-1,t_Xi)*pow(-1,t_Eta);


		    Eta=Vett_1_Eta[t_Eta];

		    Parametri(Xi,Eta,QM,KXiM[t_Eta],KEtaM[t_Xi]);

		    if(z < 0)
		    {

			IS_A(burgers_vector,FT,sign,Xi,Eta,z);
			IS_B(burgers_vector,FT,sign,Xi,Eta,z);
			IS_C(burgers_vector,FT,sign,Xi,Eta,z);

			Eta=Vett_2_Eta[t_Eta];

			Parametri(Xi,Eta,QP,KXiP[t_Eta],KEtaP[t_Xi]);

			RS_A(burgers_vector,FT,sign,Xi,Eta,z);

		    }
		    else if(z == 0)
		    {

			IS_A_Z0(burgers_vector,FT,sign,Xi,Eta,z);
			IS_B(burgers_vector,FT,sign,Xi,Eta,z);
			IS_u_C(burgers_vector,FT,sign,Xi,Eta,z);

		    }
		    else if(z > 0)
		    {

			cout << "Error: z must be <=0" << endl;

			return;

		    }

		}

	    }

	    Derivates(burgers_vector,FT,z);

	}


	double S_Strain_components[6],D_Strain_components[6],T_Strain_components[6];

	Strain_tensor(S_Strain_components,FT[0]);
	Strain_tensor(D_Strain_components,FT[1]);
	Strain_tensor(T_Strain_components,FT[2]);

	Strain[0] = coefpi*(burgers_vector[0]*S_Strain_components[0]+burgers_vector[1]*D_Strain_components[0]+burgers_vector[2]*T_Strain_components[0]);
	Strain[1] = coefpi*(burgers_vector[0]*S_Strain_components[1]+burgers_vector[1]*D_Strain_components[1]+burgers_vector[2]*T_Strain_components[1]);
	Strain[2] = coefpi*(burgers_vector[0]*S_Strain_components[2]+burgers_vector[1]*D_Strain_components[2]+burgers_vector[2]*T_Strain_components[2]);
	Strain[3] = coefpi*(burgers_vector[0]*S_Strain_components[3]+burgers_vector[1]*D_Strain_components[3]+burgers_vector[2]*T_Strain_components[3]);
	Strain[4] = coefpi*(burgers_vector[0]*S_Strain_components[4]+burgers_vector[1]*D_Strain_components[4]+burgers_vector[2]*T_Strain_components[4]);
	Strain[5] = coefpi*(burgers_vector[0]*S_Strain_components[5]+burgers_vector[1]*D_Strain_components[5]+burgers_vector[2]*T_Strain_components[5]);

}



// # Componenti dello sforzo generate da una dislocazione rettangolare (Okada 1992)

void BE3D::Stress_Okada(double S_Stress_components[],double D_Stress_components[],double T_Stress_components[],
						double burgers_vector[],double x,double y,double z,double c,double delta_gradi,
						double AL1,double AL2,double AW1,double AW2,double mu,double lambda)
{
	
	Par(SIN_DELTA,COS_DELTA,alpha,coef_alpha_1,coef_alpha_2,delta_gradi,mu,lambda);

	struct_var FT[3];
	
	for(int ti=0; ti < 3; ti++)
	{
	    initialize(FT[ti]);
	}
	
	double Xi,Eta;	
	double Vett_Xi[2];
	double Vett_1_Eta[2],Vett_2_Eta[2];
	int KXiM[2],KEtaM[2]; 
	int KXiP[2],KEtaP[2]; 
	double QM=0,QP=0;
	int flag_QM=0,flag_QP=0;
	Singular_cases(x,y,z,c,AL1,AL2,AW1,AW2,Vett_Xi,Vett_1_Eta,Vett_2_Eta,flag_QM,flag_QP,KXiM,KEtaM,KXiP,KEtaP,QM,QP);
 
	  
	if((flag_QM == 0) && (flag_QP == 0))
	{
	  
	    double sign;	  
		  
	    for(int t_Xi=0; t_Xi < 2; t_Xi++)
	    {
	      
	    	Xi=Vett_Xi[t_Xi];

			for(int t_Eta=0; t_Eta < 2; t_Eta++)
			{

				sign=pow(-1,t_Xi)*pow(-1,t_Eta);
	

//				if(CASE == 0)
//				{
		  
					Eta=Vett_1_Eta[t_Eta];

					Parametri(Xi,Eta,QM,KXiM[t_Eta],KEtaM[t_Xi]);

					if(z < 0)
					{

						IS_A(burgers_vector,FT,sign,Xi,Eta,z);
						IS_B(burgers_vector,FT,sign,Xi,Eta,z);
						IS_C(burgers_vector,FT,sign,Xi,Eta,z);


						Eta=Vett_2_Eta[t_Eta];

						Parametri(Xi,Eta,QP,KXiP[t_Eta],KEtaP[t_Xi]);

						RS_A(burgers_vector,FT,sign,Xi,Eta,z);

					}
					else if(z == 0)
					{

						IS_A_Z0(burgers_vector,FT,sign,Xi,Eta,z);
						IS_B(burgers_vector,FT,sign,Xi,Eta,z);
						IS_u_C(burgers_vector,FT,sign,Xi,Eta,z);

					}
					else if(z > 0)
					{

						cout << "Error: z must be <=0" << endl;

						return;

					}

//				}
//				else if(CASE == 1)
//				{
//
//					Eta=Vett_2_Eta[t_Eta];
//
//					Parametri(Xi,Eta,QP,KXiP[t_Eta],KEtaP[t_Xi]);
//
//					RS_A(burgers_vector,FT,sign,Xi,Eta,z);
//
//				}

			}
		
	    }  
		    
	    Derivates(burgers_vector,FT,z);

	}

	Stress_tensor(S_Stress_components,FT[0],mu,lambda);
	Stress_tensor(D_Stress_components,FT[1],mu,lambda);
	Stress_tensor(T_Stress_components,FT[2],mu,lambda);

}





void BE3D::DISPLACEMENT(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,
							double S_Displacement[],double D_Displacement[],double T_Displacement[])
{

	double AL1 = -PARAM.GEOM.L*0.5;
	double AL2 =  PARAM.GEOM.L*0.5;
	double AW1 = -PARAM.GEOM.W*0.5;
	double AW2 =  PARAM.GEOM.W*0.5;

	double x_0 =  PARAM.GEOM.pc[0];
	double y_0 =  PARAM.GEOM.pc[1];
	double   c = -PARAM.GEOM.pc[2];


	if(FLAG_COMP == 1)
	{

	    double mu=M_PAR.GET_mu();
	    double lambda=M_PAR.GET_lambda();

	    double delta_gradi=PARAM.GEOM.delta_degree;

	    double x = xn-x_0;
	    double y = yn-y_0;

	    Rotazione_coord(x,y,PARAM.GEOM.flag);

	    double burgers_vector[3];

	    burgers_vector[0] = B1;
	    burgers_vector[1] = B2;
	    burgers_vector[2] = B3;

	    Displacement_Okada(S_Displacement,D_Displacement,T_Displacement,
	    						burgers_vector,x,y,z,c,delta_gradi,AL1,AL2,AW1,AW2,mu,lambda);

	}
	else if(FLAG_COMP == 2)
	{

	    float f_AL1=(float) AL1;
	    float f_AL2=(float) AL2;
	    float f_AW1=(float) AW1;
	    float f_AW2=(float) AW2;

	    float f_mu=(float) M_PAR.GET_mu();
	    float f_lambda=(float) M_PAR.GET_lambda();

	    float f_z=(float) z;

	    float f_delta_gradi=(float) PARAM.GEOM.delta_degree;

	    float f_x = (float) xn-x_0;
	    float f_y = (float) yn-y_0;
	    float f_c = (float) c;

	    Rotazione_coord(f_x,f_y,PARAM.GEOM.phi_degree);

	    float UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ;
	    int IRET;
	    float Slip_x,Slip_y,Slip_z;

	    float alpha=(f_lambda+f_mu)/(f_lambda+2*f_mu);

	    float DISL1a=1;
	    float DISL2a=0;
	    float DISL3a=0;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1a,DISL2a,DISL3a,
	    			Slip_x,Slip_y,Slip_z,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    S_Displacement[0] = (double) Slip_x;
	    S_Displacement[1] = (double) Slip_y;
	    S_Displacement[2] = (double) Slip_z;


	    float DISL1b=0;
	    float DISL2b=1;
	    float DISL3b=0;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1b,DISL2b,DISL3b,
	    			Slip_x,Slip_y,Slip_z,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    D_Displacement[0] = (double) Slip_x;
	    D_Displacement[1] = (double) Slip_y;
	    D_Displacement[2] = (double) Slip_z;


	    float DISL1c=0;
	    float DISL2c=0;
	    float DISL3c=1;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1c,DISL2c,DISL3c,
	    			Slip_x,Slip_y,Slip_z,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    T_Displacement[0] = (double) Slip_x;
	    T_Displacement[1] = (double) Slip_y;
	    T_Displacement[2] = (double) Slip_z;

	}

	Rotazione_vettore(S_Displacement,PARAM.GEOM.phi_degree);
	Rotazione_vettore(D_Displacement,PARAM.GEOM.phi_degree);
	Rotazione_vettore(T_Displacement,PARAM.GEOM.phi_degree);

}




void BE3D::DISPLACEMENT(int FLAG_COMP,MEDIUM M_PAR,double xn,double yn,double z,double Displacement[])
{

	double AL1 = -PARAM.GEOM.L*0.5;
	double AL2 =  PARAM.GEOM.L*0.5;
	double AW1 = -PARAM.GEOM.W*0.5;
	double AW2 =  PARAM.GEOM.W*0.5;

	double x_0 =  PARAM.GEOM.pc[0];
	double y_0 =  PARAM.GEOM.pc[1];
	double   c = -PARAM.GEOM.pc[2];


    double S_Displacement_1[3],D_Displacement_1[3],T_Displacement_1[3];

	if(FLAG_COMP == 1)
	{

	    double mu=M_PAR.GET_mu();
	    double lambda=M_PAR.GET_lambda();

	    double delta_gradi=PARAM.GEOM.delta_degree;

	    double x = xn-x_0;
	    double y = yn-y_0;

	    Rotazione_coord(x,y,PARAM.GEOM.phi_degree);


	    double burgers_vector[3];

	    burgers_vector[0] = B1;
	    burgers_vector[1] = B2;
	    burgers_vector[2] = B3;

	    Displacement_Okada(S_Displacement_1,D_Displacement_1,T_Displacement_1,
	    					burgers_vector,x,y,z,c,delta_gradi,AL1,AL2,AW1,AW2,mu,lambda);

// # Spostamento generato dalla dislocazione con vettore di burger assegnato

	    Displacement[0]=S_Displacement_1[0]*B1+D_Displacement_1[0]*B2+T_Displacement_1[0]*B3;
	    Displacement[1]=S_Displacement_1[1]*B1+D_Displacement_1[1]*B2+T_Displacement_1[1]*B3;
	    Displacement[2]=S_Displacement_1[2]*B1+D_Displacement_1[2]*B2+T_Displacement_1[2]*B3;

	}
	else if(FLAG_COMP == 2)
	{

	    float f_AL1=(float) AL1;
	    float f_AL2=(float) AL2;
	    float f_AW1=(float) AW1;
	    float f_AW2=(float) AW2;

	    float f_mu=(float) M_PAR.GET_mu();
	    float f_lambda=(float) M_PAR.GET_lambda();

	    float f_z=(float) z;

	    float f_delta_gradi=(float) PARAM.GEOM.delta_degree;

	    float f_x = (float) xn-x_0;
	    float f_y = (float) yn-y_0;
	    float f_c = (float) c;

	    Rotazione_coord(f_x,f_y,PARAM.GEOM.phi_degree);


	    float alpha=(f_lambda+f_mu)/(f_lambda+2*f_mu);

	    float UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ;
	    int IRET;
	    float Slip_x,Slip_y,Slip_z;

	    float DISL1a = (float) B1;
	    float DISL2a = (float) B2;
	    float DISL3a = (float) B3;

	    DC3D(alpha,f_x,f_y,f_z,f_c,f_delta_gradi,f_AL1,f_AL2,f_AW1,f_AW2,DISL1a,DISL2a,DISL3a,
	    		Slip_x,Slip_y,Slip_z,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET);

	    Displacement[0] = (double) Slip_x;
	    Displacement[1] = (double) Slip_y;
	    Displacement[2] = (double) Slip_z;

	}

	Rotazione_vettore(Displacement,PARAM.GEOM.phi_degree);

}




// # Componenti dello spostamento generate da una dislocazione rettangolare (Okada 1992)

void BE3D::Displacement_Okada(double S_Displacement[],double D_Displacement[],double T_Displacement[],
								double burgers_vector[],double x,double y,double z,double c,double delta_gradi,
								double AL1,double AL2,double AW1,double AW2,double mu,double lambda)
{
	
	Par(SIN_DELTA,COS_DELTA,alpha,coef_alpha_1,coef_alpha_2,delta_gradi,mu,lambda);	
	
	struct_var FT[3];
	
	for(int ti=0; ti < 3; ti++)
	{
	    initialize(FT[ti]);
	}
	
	int t_Xi,t_Eta;
	double Xi,Eta;	
	double Vett_Xi[2];
	double Vett_1_Eta[2],Vett_2_Eta[2];
	int KXiM[2],KEtaM[2]; 
	int KXiP[2],KEtaP[2]; 
	double QM=0,QP=0;
	int flag_QM=0,flag_QP=0;	
	Singular_cases(x,y,z,c,AL1,AL2,AW1,AW2,Vett_Xi,Vett_1_Eta,Vett_2_Eta,flag_QM,flag_QP,KXiM,KEtaM,KXiP,KEtaP,QM,QP);

	
	if((flag_QM == 0) && (flag_QP == 0))
	{
	  		
	    double sign;
	  
	  
	    for (t_Xi=0; t_Xi < 2; t_Xi++)
	    {

	    	Xi=Vett_Xi[t_Xi];


	    	for (t_Eta=0; t_Eta < 2; t_Eta++)
	    	{

	    		sign=pow(-1,t_Xi)*pow(-1,t_Eta);


//	    		if(CASE == 0)
//	    		{

	    			// # Real source contribution

	    			Eta=Vett_1_Eta[t_Eta];

	    			Parametri(Xi,Eta,QM,KXiM[t_Eta],KEtaM[t_Xi]);

	    			if(z < 0)
	    			{

	    				IS_u_A(burgers_vector,FT,sign,Xi,Eta,z);
	    				IS_u_B(burgers_vector,FT,sign,Xi,Eta,z);
	    				IS_u_C(burgers_vector,FT,sign,Xi,Eta,z);

	    				// # Image source contribution

	    				Eta=Vett_2_Eta[t_Eta];

	    				Parametri(Xi,Eta,QP,KXiP[t_Eta],KEtaP[t_Xi]);

	    				RS_u_A(burgers_vector,FT,sign,Xi,Eta,z);

	    			}
	    			else if(z == 0)
	    			{

	    				IS_u_B(burgers_vector,FT,sign,Xi,Eta,z);

	    			}
	    			else if(z > 0)
	    			{

	    				cout << "Error: z must be <=0" << endl;

	    				return;

	    			}

//	    		}
//	    		else if(CASE == 1)
//	    		{
//
//	    			// # Image source contribution
//
//	    			Eta=Vett_2_Eta[t_Eta];
//
//	    			Parametri(Xi,Eta,QP,KXiP[t_Eta],KEtaP[t_Xi]);
//
//	    			RS_u_A(burgers_vector,FT,sign,Xi,Eta,z);
//
//	    		}

	    	}

	    }

	    CDisplacement(burgers_vector,FT,z);

	}

	
// # Spostamento generato dalla dislocazioni di tipo S, D e T	
	
	S_Displacement[0]=coefpi*FT[0].Ux;
	S_Displacement[1]=coefpi*FT[0].Uy;
	S_Displacement[2]=coefpi*FT[0].Uz;	
	
	D_Displacement[0]=coefpi*FT[1].Ux;
	D_Displacement[1]=coefpi*FT[1].Uy;
	D_Displacement[2]=coefpi*FT[1].Uz;	
	
	T_Displacement[0]=coefpi*FT[2].Ux;
	T_Displacement[1]=coefpi*FT[2].Uy;
	T_Displacement[2]=coefpi*FT[2].Uz;	
	
}




void BE3D::Parametri(double Xi,double Eta,double QQ,int KXi,int KEta)
{

	q=QQ;
	
	if(abs(q) < EPS)
	{
	    q=0e0;
	}	
	
	q2=q*q;
	q3=q2*q;
	
	if(abs(Xi) < EPS)
	{
	    Xi=0e0;
	}
	
	if(abs(Eta) < EPS)
	{
	    Eta=0e0;
	}

	Xi2=Xi*Xi;
	Eta2=Eta*Eta;	
	
	R2=Xi2+Eta2+q2;
	R=sqrt(R2);
	
	if(R == 0e0)
	{
	    return; 
	}

	R3=R*R2;
	R5=R3*R2;
	
	X=sqrt(Xi2+q2);
	
	if(q != 0e0)
	{
	    Teta=atan(Xi*Eta/(q*R));
	}
	else
	{
	    Teta=0;
	}
	

	Y=Eta*COS_DELTA+q*SIN_DELTA;
	Y2=Y*Y;	

	D=Eta*SIN_DELTA-q*COS_DELTA;
	D2=D*D;

	
	RpD=R+D;
	RpD2=(R+D)*(R+D);

	D_11=1.E0/(R*RpD);

	double RpXi;

	if(KXi == 0)	
	{
	    RpXi=R+Xi;
	    
	    X_11=1.E0/(R*RpXi);
	    X_32=(2.E0*R+Xi)*X_11*X_11/R;     
	    X_53=(8.E0*R2+9.E0*R*Xi+3.E0*Xi2)*X_11*X_11*X_11/R2;    
	    
	    Sing_iii=log(RpXi);	    
	}
	else
	{
	    X_11=0e0;	  
	    X_32=0e0;
	    X_53=0e0;		
	    
	    Sing_iii=-log(R-Xi);	    
	}
	
	
	double RpEta;	
	
	if(KEta == 0)
	{
	    RpEta=R+Eta;
	    
	    Y_11=1.E0/(R*RpEta);
	    Y_32=(2.E0*R+Eta)*Y_11*Y_11/R;  
	    Y_53=(8.E0*R2+9.E0*R*Eta+3.E0*Eta2)*Y_11*Y_11*Y_11/R2;	
	    
	    Y_0=Y_11-Xi2*Y_32;
	    
	    Sing_iv=log(RpEta);	        
	}
	else
	{
	    Y_11=0e0;
	    Y_32=0e0;
	    Y_53=0e0;
	    
	    Y_0=0e0;
	    
	    Sing_iv=-log(R-Eta);	    
	}

	
	E=SIN_DELTA/R-Y*q/R3;
	F=D/R3+Xi2*Y_32*SIN_DELTA;
	G=2.E0*X_11*SIN_DELTA-Y*q*X_32;
	H=D*q*X_32+Xi*q*Y_32*SIN_DELTA;
	
	Ep=COS_DELTA/R+D*q/R3;
	Fp=Y/R3+Xi2*Y_32*COS_DELTA;
	Gp=2.E0*X_11*COS_DELTA+D*q*X_32;
	Hp=Y*q*X_32+Xi*q*Y_32*COS_DELTA;

}




void BE3D::Par(double &SIN_DELTA,double &COS_DELTA,double &alpha,double &coef_alpha_1,double &coef_alpha_2,
					double delta_gradi,double mu,double lambda)
{
  
	float delta;
	delta=pi*(delta_gradi/180e0);

	if((delta_gradi == 90) || (delta_gradi == -90))
	{
	    SIN_DELTA=sin(delta);
	    COS_DELTA=0;
	}
	else if((delta_gradi == 180) || (delta_gradi == -180) || (delta_gradi == 0))
	{
	    SIN_DELTA=0;
	    COS_DELTA=cos(delta);
	}
	else
	{
	    SIN_DELTA=sin(delta);
	    COS_DELTA=cos(delta);
	}

	SIN_DELTA_2=pow(SIN_DELTA,2);
	COS_DELTA_2=pow(COS_DELTA,2);

	alpha=(lambda+mu)/(lambda+2*mu);
	coef_alpha_1=1.E0-alpha;
	coef_alpha_2=coef_alpha_1/alpha;

}	
	
	
void BE3D::initialize(struct_var &struttura)
{
  
	struttura.j1A=0,struttura.j2A=0,struttura.j3A=0;
	struttura.j1B=0,struttura.j2B=0,struttura.j3B=0;
	struttura.j1C=0,struttura.j2C=0,struttura.j3C=0;
	struttura.J1A=0,struttura.J2A=0,struttura.J3A=0;
	
	struttura.k1A=0,struttura.k2A=0,struttura.k3A=0;
	struttura.k1B=0,struttura.k2B=0,struttura.k3B=0;
	struttura.k1C=0,struttura.k2C=0,struttura.k3C=0;
	struttura.K1A=0,struttura.K2A=0,struttura.K3A=0;	
	
	struttura.l1A=0,struttura.l2A=0,struttura.l3A=0;
	struttura.l1B=0,struttura.l2B=0,struttura.l3B=0;	
	struttura.l1C=0,struttura.l2C=0,struttura.l3C=0;
	struttura.L1A=0,struttura.L2A=0,struttura.L3A=0;	
	

	struttura.Ux_x=0,struttura.Uy_y=0,struttura.Uz_z=0;
	struttura.Ux_y=0,struttura.Uy_x=0,struttura.Ux_z=0;
	struttura.Uz_x=0,struttura.Uy_z=0,struttura.Uz_y=0;  
	
	struttura.u1A=0,struttura.U1A=0,struttura.u1B=0,struttura.u1C=0;
	struttura.u2A=0,struttura.U2A=0,struttura.u2B=0,struttura.u2C=0;
	struttura.u3A=0,struttura.U3A=0,struttura.u3B=0,struttura.u3C=0;
	
	struttura.Ux=0,struttura.Uy=0,struttura.Uz=0;

}	
	


void BE3D::Singular_cases(double x,double y,double z,double c,double AL1,double AL2,double AW1,double AW2,
							double Vett_Xi[],double Vett_1_Eta[],double Vett_2_Eta[],
							int &flag_QM,int &flag_QP,int KXiM[],int KEtaM[],int KXiP[],int KEtaP[],double &QM,double &QP)
{
	
	Vett_Xi[0]=x-AL1;
	Vett_Xi[1]=x-AL2;
	
// ----------------------------------------------                                       
// ----- REJECT SINGULAR CASE ON FAULT EDGE -----                                   
// ----------------------------------------------                                      
	
	const double F0=0; 
	
	if(abs(Vett_Xi[0]) < EPS)
	{
	  Vett_Xi[0]=F0;
	}
	if(abs(Vett_Xi[1]) < EPS)
	{
	  Vett_Xi[1]=F0;
	}	
	

	double PM=0;

	PM=y*COS_DELTA+(c-z)*SIN_DELTA;	
	QM=y*SIN_DELTA-(c-z)*COS_DELTA;	
	
	Vett_1_Eta[0]=PM-AW1;		
	Vett_1_Eta[1]=PM-AW2;
	

	if(abs(QM) < EPS)
	{
	  QM=F0;
	}	
	if(abs(Vett_1_Eta[0]) < EPS)
	{
	  Vett_1_Eta[0]=F0;
	}
	if(abs(Vett_1_Eta[1]) < EPS)
	{
	  Vett_1_Eta[1]=F0;
	}


	double PP=0;
	double zm=-z;
	
	PP=y*COS_DELTA+(c-zm)*SIN_DELTA;
	QP=y*SIN_DELTA-(c-zm)*COS_DELTA;

	Vett_2_Eta[0]=PP-AW1;		
	Vett_2_Eta[1]=PP-AW2;
	
	
	if(abs(QP) < EPS)
	{
	  QP=F0;
	}
	if(abs(Vett_2_Eta[0]) < EPS)
	{
	  Vett_2_Eta[0]=F0;
	}
	if(abs(Vett_2_Eta[1]) < EPS)
	{
	  Vett_2_Eta[1]=F0;
	}	
	

	if(QM == F0 && ((Vett_Xi[0]*Vett_Xi[1] <= F0 && Vett_1_Eta[0]*Vett_1_Eta[1] == F0) || (Vett_1_Eta[0]*Vett_1_Eta[1] <= F0 && Vett_Xi[0]*Vett_Xi[1] == F0)))
	{
	  flag_QM=1;
	}

	if(QP == F0 && ((Vett_Xi[0]*Vett_Xi[1] <= F0 && Vett_2_Eta[0]*Vett_2_Eta[1] == F0) || (Vett_2_Eta[0]*Vett_2_Eta[1] <= F0 && Vett_Xi[0]*Vett_Xi[1] == F0)))
	{
	  flag_QP=1;
	}	
	

// ----- ON NEGATIVE EXTENSION OF FAULT EDGE

	double R12,R21,R22;
	
	KXiM[0]=0;                                                          
	KXiM[1]=0;                                                          
	KEtaM[0]=0;                                                          
	KEtaM[1]=0;    
	
	R12=sqrt(Vett_Xi[0]*Vett_Xi[0]+Vett_1_Eta[1]*Vett_1_Eta[1]+QM*QM);                            
	R21=sqrt(Vett_Xi[1]*Vett_Xi[1]+Vett_1_Eta[0]*Vett_1_Eta[0]+QM*QM);                            
	R22=sqrt(Vett_Xi[1]*Vett_Xi[1]+Vett_1_Eta[1]*Vett_1_Eta[1]+QM*QM); 
	
	if(Vett_Xi[0] < F0 && R21+Vett_Xi[1] < EPS)
	{
	    KXiM[0]=1;
	}
	if(Vett_Xi[0] < F0 && R22+Vett_Xi[1] < EPS)
	{
	    KXiM[1]=1;
	}
	
	if(Vett_1_Eta[0] < F0 && R12+Vett_1_Eta[1] < EPS)
	{
	    KEtaM[0]=1;
	}
	if(Vett_1_Eta[0] < F0 && R22+Vett_1_Eta[1] < EPS)
	{
	    KEtaM[1]=1;
	}
	
	
	KXiP[0]=0;                                                          
	KXiP[1]=0;                                                          
	KEtaP[0]=0;                                                          
	KEtaP[1]=0;    
	
	R12=sqrt(Vett_Xi[0]*Vett_Xi[0]+Vett_2_Eta[1]*Vett_2_Eta[1]+QP*QP);                            
	R21=sqrt(Vett_Xi[1]*Vett_Xi[1]+Vett_2_Eta[0]*Vett_2_Eta[0]+QP*QP);                            
	R22=sqrt(Vett_Xi[1]*Vett_Xi[1]+Vett_2_Eta[1]*Vett_2_Eta[1]+QP*QP); 
	
	if(Vett_Xi[0] < F0 && R21+Vett_Xi[1] < EPS)
	{
	    KXiP[0]=1;
	}
	if(Vett_Xi[0] < F0 && R22+Vett_Xi[1] < EPS)
	{
	    KXiP[1]=1;
	}
	
	if(Vett_2_Eta[0] < F0 && R12+Vett_2_Eta[1] < EPS)
	{
	    KEtaP[0]=1;
	}
	if(Vett_2_Eta[0] < F0 && R22+Vett_2_Eta[1] < EPS)
	{
	    KEtaP[1]=1;
	}	
	
}



// # Funzioni

void BE3D::IS_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	if(burgers_vector[0] != 0)
	{
	    FT[0].j1A += sign * fSX1A(Xi,Eta,z);
	    FT[0].j2A += sign * fSX2A(Xi,Eta,z);
	    FT[0].j3A += sign * fSX3A(Xi,Eta,z);
	    FT[0].k1A += sign * fSY1A(Xi,Eta,z);
	    FT[0].k2A += sign * fSY2A(Xi,Eta,z);
	    FT[0].k3A += sign * fSY3A(Xi,Eta,z);
	    FT[0].l1A += sign * fSZ1A(Xi,Eta,z);
	    FT[0].l2A += sign * fSZ2A(Xi,Eta,z);
	    FT[0].l3A += sign * fSZ3A(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0)
	{    
	    FT[1].j1A += sign * fDX1A(Xi,Eta,z);
	    FT[1].j2A += sign * fDX2A(Xi,Eta,z);
	    FT[1].j3A += sign * fDX3A(Xi,Eta,z);
	    FT[1].k1A += sign * fDY1A(Xi,Eta,z);
	    FT[1].k2A += sign * fDY2A(Xi,Eta,z);
	    FT[1].k3A += sign * fDY3A(Xi,Eta,z);
	    FT[1].l1A += sign * fDZ1A(Xi,Eta,z);
	    FT[1].l2A += sign * fDZ2A(Xi,Eta,z);
	    FT[1].l3A += sign * fDZ3A(Xi,Eta,z);
	}
	if(burgers_vector[2] != 0)
	{    
	    FT[2].j1A += sign * fTX1A(Xi,Eta,z);
	    FT[2].j2A += sign * fTX2A(Xi,Eta,z);
	    FT[2].j3A += sign * fTX3A(Xi,Eta,z);
	    FT[2].k1A += sign * fTY1A(Xi,Eta,z);
	    FT[2].k2A += sign * fTY2A(Xi,Eta,z);
	    FT[2].k3A += sign * fTY3A(Xi,Eta,z);
	    FT[2].l1A += sign * fTZ1A(Xi,Eta,z);
	    FT[2].l2A += sign * fTZ2A(Xi,Eta,z);
	    FT[2].l3A += sign * fTZ3A(Xi,Eta,z);
	}    
}



void BE3D::IS_A_Z0(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	if(burgers_vector[0] != 0)
	{
	    FT[0].l1A += sign * fSZ1A(Xi,Eta,z);
	    FT[0].l2A += sign * fSZ2A(Xi,Eta,z);
	    FT[0].l3A += sign * fSZ3A(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0)
	{    
	    FT[1].l1A += sign * fDZ1A(Xi,Eta,z);
	    FT[1].l2A += sign * fDZ2A(Xi,Eta,z);
	    FT[1].l3A += sign * fDZ3A(Xi,Eta,z);
	}
	if(burgers_vector[2] != 0)
	{    
	    FT[2].l1A += sign * fTZ1A(Xi,Eta,z);
	    FT[2].l2A += sign * fTZ2A(Xi,Eta,z);
	    FT[2].l3A += sign * fTZ3A(Xi,Eta,z);
	}    
}



void BE3D::RS_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	if(burgers_vector[0] != 0)
	{
	    FT[0].J1A += sign * fSX1A(Xi,Eta,-z);
	    FT[0].J2A += sign * fSX2A(Xi,Eta,-z);
	    FT[0].J3A += sign * fSX3A(Xi,Eta,-z);
	    FT[0].K1A += sign * fSY1A(Xi,Eta,-z);
	    FT[0].K2A += sign * fSY2A(Xi,Eta,-z);
	    FT[0].K3A += sign * fSY3A(Xi,Eta,-z);
	    FT[0].L1A += sign * fSZ1A(Xi,Eta,-z);
	    FT[0].L2A += sign * fSZ2A(Xi,Eta,-z);
	    FT[0].L3A += sign * fSZ3A(Xi,Eta,-z);
	}
	if(burgers_vector[1] != 0)
	{    
	    FT[1].J1A += sign * fDX1A(Xi,Eta,-z);
	    FT[1].J2A += sign * fDX2A(Xi,Eta,-z);
	    FT[1].J3A += sign * fDX3A(Xi,Eta,-z);
	    FT[1].K1A += sign * fDY1A(Xi,Eta,-z);
	    FT[1].K2A += sign * fDY2A(Xi,Eta,-z);
	    FT[1].K3A += sign * fDY3A(Xi,Eta,-z);
	    FT[1].L1A += sign * fDZ1A(Xi,Eta,-z);
	    FT[1].L2A += sign * fDZ2A(Xi,Eta,-z);
	    FT[1].L3A += sign * fDZ3A(Xi,Eta,-z);
	}
	if(burgers_vector[2] != 0)
	{    
	    FT[2].J1A += sign * fTX1A(Xi,Eta,-z);
	    FT[2].J2A += sign * fTX2A(Xi,Eta,-z);
	    FT[2].J3A += sign * fTX3A(Xi,Eta,-z);
	    FT[2].K1A += sign * fTY1A(Xi,Eta,-z);
	    FT[2].K2A += sign * fTY2A(Xi,Eta,-z);
	    FT[2].K3A += sign * fTY3A(Xi,Eta,-z);
	    FT[2].L1A += sign * fTZ1A(Xi,Eta,-z);
	    FT[2].L2A += sign * fTZ2A(Xi,Eta,-z);
	    FT[2].L3A += sign * fTZ3A(Xi,Eta,-z);
	}    
}



void BE3D::IS_B(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{

	if(COS_DELTA != 0e0)
	{
	    K_1=(Xi/COS_DELTA)*(D_11-Y_11*SIN_DELTA);
	    K_3=(q*Y_11-Y*D_11)/COS_DELTA;		
	}
	else
	{
	    K_1=Xi*q*D_11/RpD;
	    K_3=(SIN_DELTA/RpD)*(Xi2*D_11-1);		
	}

	K_2=1.E0/R+K_3*SIN_DELTA;	
	K_4=Xi*Y_11*COS_DELTA-K_1*SIN_DELTA;	
	

	J_5=-(D+Y2/RpD)*D_11;
	
	if(COS_DELTA != 0e0)
	{
	    J_6=(K_3-J_5*SIN_DELTA)/COS_DELTA;
	}
	else
	{
	    J_6=-Y*(Xi2*D_11-0.5)/RpD2;
	}

	J_1=J_5*COS_DELTA-J_6*SIN_DELTA;
	J_2=(Xi*Y/RpD)*D_11;

	if(COS_DELTA != 0e0)
	{
	    J_3=(K_1-J_2*SIN_DELTA)/COS_DELTA;
	}
	else
	{
	    J_3=-Xi*(q2*D_11-0.5)/RpD2;
	}

	J_4=-Xi*Y_11-J_2*COS_DELTA+J_3*SIN_DELTA;

	
	if(burgers_vector[0] != 0)
	{
	    FT[0].j1B += sign * fSX1B(Xi,Eta,z);
	    FT[0].j2B += sign * fSX2B(Xi,Eta,z);
	    FT[0].j3B += sign * fSX3B(Xi,Eta,z);
	    FT[0].k1B += sign * fSY1B(Xi,Eta,z);
	    FT[0].k2B += sign * fSY2B(Xi,Eta,z);
	    FT[0].k3B += sign * fSY3B(Xi,Eta,z);
	    FT[0].l1B += sign * fSZ1B(Xi,Eta,z);
	    FT[0].l2B += sign * fSZ2B(Xi,Eta,z);
	    FT[0].l3B += sign * fSZ3B(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0)
	{    
	    FT[1].j1B += sign * fDX1B(Xi,Eta,z);
	    FT[1].j2B += sign * fDX2B(Xi,Eta,z);
	    FT[1].j3B += sign * fDX3B(Xi,Eta,z);
	    FT[1].k1B += sign * fDY1B(Xi,Eta,z);
	    FT[1].k2B += sign * fDY2B(Xi,Eta,z);
	    FT[1].k3B += sign * fDY3B(Xi,Eta,z);
	    FT[1].l1B += sign * fDZ1B(Xi,Eta,z);
	    FT[1].l2B += sign * fDZ2B(Xi,Eta,z);
	    FT[1].l3B += sign * fDZ3B(Xi,Eta,z);
	}
	if(burgers_vector[2] != 0)
	{    
	    FT[2].j1B += sign * fTX1B(Xi,Eta,z);
	    FT[2].j2B += sign * fTX2B(Xi,Eta,z);
	    FT[2].j3B += sign * fTX3B(Xi,Eta,z);
	    FT[2].k1B += sign * fTY1B(Xi,Eta,z);
	    FT[2].k2B += sign * fTY2B(Xi,Eta,z);
	    FT[2].k3B += sign * fTY3B(Xi,Eta,z);
	    FT[2].l1B += sign * fTZ1B(Xi,Eta,z);
	    FT[2].l2B += sign * fTZ2B(Xi,Eta,z);
	    FT[2].l3B += sign * fTZ3B(Xi,Eta,z);
	}    

	
}



void BE3D::IS_C(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
  
	C=D+z;  

	h=q*COS_DELTA-z;

	Z_32=SIN_DELTA/R3-h*Y_32;
	Z_53=3.E0*SIN_DELTA/R5-h*Y_53;	
	Z_0=Z_32-Xi2*Z_53;  
	
	CR=3.E0*C/R5;
	
	P=COS_DELTA/R3+q*Y_32*SIN_DELTA;
	Q=CR*D-(z*Y_32+Z_32+Z_0)*SIN_DELTA;
	
	Pp=SIN_DELTA/R3-q*Y_32*COS_DELTA;
	Qp=CR*Y+q*Y_32-(z*Y_32+Z_32+Z_0)*COS_DELTA;
	
	qR=3.E0*q/R5;      
	CDR=(C+D)/R3;            
	YY0=Y/R3-Y_0*COS_DELTA;	
	
	if(burgers_vector[0] != 0)
	{   
	    FT[0].j1C += sign * fSX1C(Xi,Eta,z);
	    FT[0].j2C += sign * fSX2C(Xi,Eta,z);
	    FT[0].j3C += sign * fSX3C(Xi,Eta,z);
	    FT[0].k1C += sign * fSY1C(Xi,Eta,z);
	    FT[0].k2C += sign * fSY2C(Xi,Eta,z);
	    FT[0].k3C += sign * fSY3C(Xi,Eta,z);
	    FT[0].l1C += sign * fSZ1C(Xi,Eta,z);
	    FT[0].l2C += sign * fSZ2C(Xi,Eta,z);
	    FT[0].l3C += sign * fSZ3C(Xi,Eta,z);
	    FT[0].u1C += sign * fSu1C(Xi,Eta,z);
	    FT[0].u2C += sign * fSu2C(Xi,Eta,z);
	    FT[0].u3C += sign * fSu3C(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0)
	{  
	    FT[1].j1C += sign * fDX1C(Xi,Eta,z);
	    FT[1].j2C += sign * fDX2C(Xi,Eta,z);
	    FT[1].j3C += sign * fDX3C(Xi,Eta,z);
	    FT[1].k1C += sign * fDY1C(Xi,Eta,z);
	    FT[1].k2C += sign * fDY2C(Xi,Eta,z);
	    FT[1].k3C += sign * fDY3C(Xi,Eta,z);
	    FT[1].l1C += sign * fDZ1C(Xi,Eta,z);
	    FT[1].l2C += sign * fDZ2C(Xi,Eta,z);
	    FT[1].l3C += sign * fDZ3C(Xi,Eta,z);
	    FT[1].u1C += sign * fDu1C(Xi,Eta,z);
	    FT[1].u2C += sign * fDu2C(Xi,Eta,z);
	    FT[1].u3C += sign * fDu3C(Xi,Eta,z);
	}
	if(burgers_vector[2] != 0)
	{
	    FT[2].j1C += sign * fTX1C(Xi,Eta,z);
	    FT[2].j2C += sign * fTX2C(Xi,Eta,z);
	    FT[2].j3C += sign * fTX3C(Xi,Eta,z);
	    FT[2].k1C += sign * fTY1C(Xi,Eta,z);
	    FT[2].k2C += sign * fTY2C(Xi,Eta,z);
	    FT[2].k3C += sign * fTY3C(Xi,Eta,z);
	    FT[2].l1C += sign * fTZ1C(Xi,Eta,z);
	    FT[2].l2C += sign * fTZ2C(Xi,Eta,z);
	    FT[2].l3C += sign * fTZ3C(Xi,Eta,z);
	    FT[2].u1C += sign * fTu1C(Xi,Eta,z);
	    FT[2].u2C += sign * fTu2C(Xi,Eta,z);
	    FT[2].u3C += sign * fTu3C(Xi,Eta,z);    
	}
    
}




// # Spostamento generato da una dislocazione di tipo SF-slip

void BE3D::IS_u_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	if(burgers_vector[0] != 0)
	{  
	    FT[0].u1A += sign * fSu1A(Xi,Eta,z);
	    FT[0].u2A += sign * fSu2A(Xi,Eta,z);
	    FT[0].u3A += sign * fSu3A(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0) 
	{
	    FT[1].u1A += sign * fDu1A(Xi,Eta,z);
	    FT[1].u2A += sign * fDu2A(Xi,Eta,z);
	    FT[1].u3A += sign * fDu3A(Xi,Eta,z);     
	}
	if(burgers_vector[2] != 0) 
	{
	    FT[2].u1A += sign * fTu1A(Xi,Eta,z);
	    FT[2].u2A += sign * fTu2A(Xi,Eta,z);
	    FT[2].u3A += sign * fTu3A(Xi,Eta,z); 
	}
}



void BE3D::RS_u_A(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	if(burgers_vector[0] != 0)
	{  
	    FT[0].U1A += sign * fSu1A(Xi,Eta,-z);
	    FT[0].U2A += sign * fSu2A(Xi,Eta,-z);
	    FT[0].U3A += sign * fSu3A(Xi,Eta,-z);
	}
	if(burgers_vector[1] != 0)
	{
	    FT[1].U1A += sign * fDu1A(Xi,Eta,-z);
	    FT[1].U2A += sign * fDu2A(Xi,Eta,-z);
	    FT[1].U3A += sign * fDu3A(Xi,Eta,-z);      
	}
	if(burgers_vector[2] != 0)
	{
	    FT[2].U1A += sign * fTu1A(Xi,Eta,-z);
	    FT[2].U2A += sign * fTu2A(Xi,Eta,-z);
	    FT[2].U3A += sign * fTu3A(Xi,Eta,-z); 
	}
}



void BE3D::IS_u_B(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
  
	if(COS_DELTA != 0e0)
	{
	    I_3=Y/(RpD*COS_DELTA)-(Sing_iv-SIN_DELTA*log(RpD))/COS_DELTA_2;
	    
	    if(Xi == 0e0)
	    {
	    	I_4=0e0;
	    }
	    else
	    {
	    	I_4=(1.E0/COS_DELTA_2)*(SIN_DELTA*COS_DELTA*Xi/RpD+2.E0*atan((Eta*(X+q*COS_DELTA)+X*(R+X)*SIN_DELTA)/(Xi*(R+X)*COS_DELTA)));
	    }	    
	}
	else
	{
	    I_3=0.5*(Eta/RpD+Y*q/RpD2-Sing_iv);
		    
	    I_4=0.5*Xi*Y/RpD2;	    
	}
	 
	I_1=-(Xi/RpD)*COS_DELTA-I_4*SIN_DELTA;	
	I_2=log(RpD)+I_3*SIN_DELTA;	
	    
      
	if(burgers_vector[0] != 0) 
	{    
	    FT[0].u1B += sign * fSu1B(Xi,Eta,z); 
	    FT[0].u2B += sign * fSu2B(Xi,Eta,z); 
	    FT[0].u3B += sign * fSu3B(Xi,Eta,z); 
	}
	if(burgers_vector[1] != 0) 
	{
	    FT[1].u1B += sign * fDu1B(Xi,Eta,z); 
	    FT[1].u2B += sign * fDu2B(Xi,Eta,z); 
	    FT[1].u3B += sign * fDu3B(Xi,Eta,z);  
	}
	if(burgers_vector[2] != 0) 
	{   
	    FT[2].u1B += sign * fTu1B(Xi,Eta,z); 
	    FT[2].u2B += sign * fTu2B(Xi,Eta,z); 
	    FT[2].u3B += sign * fTu3B(Xi,Eta,z);  
	}
	
}



void BE3D::IS_u_C(double burgers_vector[],struct_var FT[],double sign,double Xi,double Eta,double z)
{
	C=D+z;   
	h=q*COS_DELTA-z;  
	Z_32=SIN_DELTA/R3-h*Y_32;  
  
	if(burgers_vector[0] != 0) 
	{  
	    FT[0].u1C += sign * fSu1C(Xi,Eta,z);
	    FT[0].u2C += sign * fSu2C(Xi,Eta,z);
	    FT[0].u3C += sign * fSu3C(Xi,Eta,z);
	}
	if(burgers_vector[1] != 0) 
	{
	    FT[1].u1C += sign * fDu1C(Xi,Eta,z);
	    FT[1].u2C += sign * fDu2C(Xi,Eta,z);
	    FT[1].u3C += sign * fDu3C(Xi,Eta,z);
	}
	if(burgers_vector[2] != 0) 
	{
	    FT[2].u1C += sign * fTu1C(Xi,Eta,z);
	    FT[2].u2C += sign * fTu2C(Xi,Eta,z);
	    FT[2].u3C += sign * fTu3C(Xi,Eta,z);      
	}
}



// # Funzioni per una faglia di tipo strike-slip

double BE3D::fSX1A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*q*Y_11-0.5*alpha*Xi2*q*Y_32;
	return result;
}
double BE3D::fSX2A(double Xi,double Eta,double z)
{
	double result=-0.5*alpha*Xi*q/R3;
	return result;
}
double BE3D::fSX3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Xi*Y_11+0.5*alpha*Xi*q2*Y_32;
	return result;
}
double BE3D::fSY1A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Xi*Y_11*SIN_DELTA+0.5*D*X_11+0.5*alpha*Xi*F;
	return result;
}
double BE3D::fSY2A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*E;
	return result;
}
double BE3D::fSY3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*(COS_DELTA/R+q*Y_11*SIN_DELTA)-0.5*alpha*q*F;
	return result;
}
double BE3D::fSZ1A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Xi*Y_11*COS_DELTA+0.5*Y*X_11+0.5*alpha*Xi*Fp;
	return result;
}
double BE3D::fSZ2A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*Ep;
	return result;
}
double BE3D::fSZ3A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*(SIN_DELTA/R-q*Y_11*COS_DELTA)-0.5*alpha*q*Fp;
	return result;
}


double BE3D::fSX1B(double Xi,double Eta,double z)
{
	double result=Xi2*q*Y_32-coef_alpha_2*J_1*SIN_DELTA;
	return result;
}
double BE3D::fSX2B(double Xi,double Eta,double z)
{
	double result=Xi*q/R3-coef_alpha_2*J_2*SIN_DELTA;
	return result;
}
double BE3D::fSX3B(double Xi,double Eta,double z)
{
	double result=-Xi*q2*Y_32-coef_alpha_2*J_3*SIN_DELTA;
	return result;
}
double BE3D::fSY1B(double Xi,double Eta,double z)
{
	double result=-Xi*F-D*X_11+coef_alpha_2*(Xi*Y_11+J_4)*SIN_DELTA;
	return result;
}
double BE3D::fSY2B(double Xi,double Eta,double z)
{
	double result=-E+coef_alpha_2*(1.E0/R+J_5)*SIN_DELTA;
	return result;
}
double BE3D::fSY3B(double Xi,double Eta,double z)
{
	double result=q*F-coef_alpha_2*(q*Y_11-J_6)*SIN_DELTA;
	return result;
}
double BE3D::fSZ1B(double Xi,double Eta,double z)
{
	double result=-Xi*Fp-Y*X_11+coef_alpha_2*K_1*SIN_DELTA;
	return result;
}
double BE3D::fSZ2B(double Xi,double Eta,double z)
{
	double result=-Ep+coef_alpha_2*Y*D_11*SIN_DELTA;
	return result;
}
double BE3D::fSZ3B(double Xi,double Eta,double z)
{
	double result=q*Fp+coef_alpha_2*K_2*SIN_DELTA;
	return result;
}


double BE3D::fSX1C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*Y_0*COS_DELTA-alpha*q*Z_0;
	return result;
}
double BE3D::fSX2C(double Xi,double Eta,double z)
{
// 	double result=-coef_alpha_1*Xi*(COS_DELTA/R3+2*q*Y_32*SIN_DELTA)+alpha*3.E0*C*Xi*q/R5;
	double result=-coef_alpha_1*Xi*(COS_DELTA/R3+2.E0*q*Y_32*SIN_DELTA)+alpha*C*Xi*qR;
	return result;
}
double BE3D::fSX3C(double Xi,double Eta,double z)
{
	double result=-coef_alpha_1*Xi*q*Y_32*COS_DELTA+alpha*Xi*(3.E0*C*Eta/R5-z*Y_32-Z_32-Z_0);
	return result;
}

double BE3D::fSY1C(double Xi,double Eta,double z)
{
	double result=-coef_alpha_1*Xi*P*COS_DELTA-alpha*Xi*Q;
	return result;
}
double BE3D::fSY2C(double Xi,double Eta,double z)
{
// 	double result=2.E0*coef_alpha_1*(D/R3-Y_0*SIN_DELTA)*SIN_DELTA-Y*COS_DELTA/R3-alpha*((C+D)*SIN_DELTA/R3-Eta/R3-3.E0*C*Y*q/R5);
// 	double result=2.E0*coef_alpha_1*(D/R3-Y_0*SIN_DELTA)*SIN_DELTA-Y*COS_DELTA/R3-alpha*((C+D)*SIN_DELTA/R3-Eta/R3-C*Y*qR);
	double result=2.E0*coef_alpha_1*(D/R3-Y_0*SIN_DELTA)*SIN_DELTA-Y*COS_DELTA/R3-alpha*(CDR*SIN_DELTA-Eta/R3-C*Y*qR);
	return result;
}
double BE3D::fSY3C(double Xi,double Eta,double z)
{
/*	double result=-coef_alpha_1*q/R3+(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+alpha*((C+D)*COS_DELTA/R3+3.E0*C*D*q/R5-(Y_0*COS_DELTA+q*Z_0)*SIN_DELTA);*/
// 	double result=-coef_alpha_1*q/R3+(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+alpha*((C+D)*COS_DELTA/R3+C*D*qR-(Y_0*COS_DELTA+q*Z_0)*SIN_DELTA);
// 	double result=-coef_alpha_1*q/R3+(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+alpha*(CDR*COS_DELTA+C*D*qR-(Y_0*COS_DELTA+q*Z_0)*SIN_DELTA);
	double result=-coef_alpha_1*q/R3+YY0*SIN_DELTA+alpha*(CDR*COS_DELTA+C*D*qR-(Y_0*COS_DELTA+q*Z_0)*SIN_DELTA);
	return result;
}
double BE3D::fSZ1C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*Xi*Pp*COS_DELTA-alpha*Xi*Qp;
	return result;
}
double BE3D::fSZ2C(double Xi,double Eta,double z)
{
// 	double result=2.E0*coef_alpha_1*(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+D*COS_DELTA/R3-alpha*((C+D)*COS_DELTA/R3+3.E0*C*D*q/R5);
// 	double result=2.E0*coef_alpha_1*(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+D*COS_DELTA/R3-alpha*((C+D)*COS_DELTA/R3+C*D*qR);
// 	double result=2.E0*coef_alpha_1*(Y/R3-Y_0*COS_DELTA)*SIN_DELTA+D*COS_DELTA/R3-alpha*(CDR*COS_DELTA+C*D*qR);
	double result=2.E0*coef_alpha_1*YY0*SIN_DELTA+D*COS_DELTA/R3-alpha*(CDR*COS_DELTA+C*D*qR);
	return result;
}
double BE3D::fSZ3C(double Xi,double Eta,double z)
{
// 	double result=(Y/R3-Y_0*COS_DELTA)*COS_DELTA-alpha*((C+D)*SIN_DELTA/R3-3.E0*C*Y*q/R5-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
// 	double result=(Y/R3-Y_0*COS_DELTA)*COS_DELTA-alpha*((C+D)*SIN_DELTA/R3-C*Y*qR-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
// 	double result=(Y/R3-Y_0*COS_DELTA)*COS_DELTA-alpha*(CDR*SIN_DELTA-C*Y*qR-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
	double result=YY0*COS_DELTA-alpha*(CDR*SIN_DELTA-C*Y*qR-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
	return result;
}



// # Funzioni per una faglia di tipo dip-slip

double BE3D::fDX1A(double Xi,double Eta,double z)
{
	double result=-0.5*alpha*Xi*q/R3;
	return result;
}
double BE3D::fDX2A(double Xi,double Eta,double z)
{
	double result=-0.5*q*Y_11-0.5*alpha*Eta*q/R3;
	return result;
}
double BE3D::fDX3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1/R+0.5*alpha*q2/R3;
	return result;
}
double BE3D::fDY1A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*E;
	return result;
}
double BE3D::fDY2A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*D*X_11+0.5*Xi*Y_11*SIN_DELTA+0.5*alpha*Eta*G;
	return result;
}
double BE3D::fDY3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Y*X_11-0.5*alpha*q*G;
	return result;
}
double BE3D::fDZ1A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*Ep;
	return result;
}
double BE3D::fDZ2A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Y*X_11+0.5*Xi*Y_11*COS_DELTA+0.5*alpha*Eta*Gp;
	return result;
}
double BE3D::fDZ3A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*D*X_11-0.5*alpha*q*Gp;
	return result;
}


double BE3D::fDX1B(double Xi,double Eta,double z)
{
	double result=Xi*q/R3+coef_alpha_2*J_4*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDX2B(double Xi,double Eta,double z)
{
	double result=Eta*q/R3+q*Y_11+coef_alpha_2*J_5*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDX3B(double Xi,double Eta,double z)
{
	double result=-q2/R3+coef_alpha_2*J_6*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDY1B(double Xi,double Eta,double z)
{
	double result=-E+coef_alpha_2*J_1*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDY2B(double Xi,double Eta,double z)
{
	double result=-Eta*G-Xi*Y_11*SIN_DELTA+coef_alpha_2*J_2*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDY3B(double Xi,double Eta,double z)
{
	double result=q*G+coef_alpha_2*J_3*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDZ1B(double Xi,double Eta,double z)
{
	double result=-Ep-coef_alpha_2*K_3*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDZ2B(double Xi,double Eta,double z)
{
	double result=-Eta*Gp-Xi*Y_11*COS_DELTA-coef_alpha_2*Xi*D_11*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDZ3B(double Xi,double Eta,double z)
{
	double result=q*Gp-coef_alpha_2*K_4*SIN_DELTA*COS_DELTA;
	return result;
}


double BE3D::fDX1C(double Xi,double Eta,double z)
{
// 	double result=-coef_alpha_1*(Xi/R3)*COS_DELTA+Xi*q*Y_32*SIN_DELTA+alpha*3.E0*C*Xi*q/R5;
	double result=-coef_alpha_1*(Xi/R3)*COS_DELTA+Xi*q*Y_32*SIN_DELTA+alpha*C*Xi*qR;
	return result;
}
double BE3D::fDX2C(double Xi,double Eta,double z)
{
// 	double result=-coef_alpha_1*Y/R3+alpha*3.E0*C*Eta*q/R5;
	double result=-coef_alpha_1*Y/R3+alpha*C*Eta*qR;
	return result;
}
double BE3D::fDX3C(double Xi,double Eta,double z)
{
	double result=D/R3-Y_0*SIN_DELTA+alpha*(C/R3)*(1.E0-3.E0*q2/R2);
	return result;
}
double BE3D::fDY1C(double Xi,double Eta,double z)
{
// 	double result=-coef_alpha_1*Eta/R3+Y_0*SIN_DELTA_2-alpha*((C+D)*SIN_DELTA/R3-3.E0*C*Y*q/R5);
// 	double result=-coef_alpha_1*Eta/R3+Y_0*SIN_DELTA_2-alpha*((C+D)*SIN_DELTA/R3-C*Y*qR);
	double result=-coef_alpha_1*Eta/R3+Y_0*SIN_DELTA_2-alpha*(CDR*SIN_DELTA-C*Y*qR);
	return result;
}
double BE3D::fDY2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*(X_11-Y2*X_32)-alpha*C*((D+2.E0*q*COS_DELTA)*X_32-Y*Eta*q*X_53);
	return result;
}
double BE3D::fDY3C(double Xi,double Eta,double z)
{
	double result=Xi*P*SIN_DELTA+Y*D*X_32+alpha*C*((Y+2.E0*q*SIN_DELTA)*X_32-Y*q2*X_53);
	return result;
}
double BE3D::fDZ1C(double Xi,double Eta,double z)
{
// 	double result=-q/R3+Y_0*SIN_DELTA*COS_DELTA-alpha*(((C+D)/R3)*COS_DELTA+3.E0*C*D*q/R5);
// 	double result=-q/R3+Y_0*SIN_DELTA*COS_DELTA-alpha*(((C+D)/R3)*COS_DELTA+C*D*qR);
	double result=-q/R3+Y_0*SIN_DELTA*COS_DELTA-alpha*(CDR*COS_DELTA+C*D*qR);
	return result;
}
double BE3D::fDZ2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*Y*D*X_32-alpha*C*((Y-2.E0*q*SIN_DELTA)*X_32+D*Eta*q*X_53);
	return result;
}
double BE3D::fDZ3C(double Xi,double Eta,double z)
{
	double result=-Xi*Pp*SIN_DELTA+X_11-D2*X_32-alpha*C*((D-2.E0*q*COS_DELTA)*X_32-D*q2*X_53);
	return result;
}



// # Funzioni per una faglia di tipo tensile

double BE3D::fTX1A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*Xi*Y_11+0.5*alpha*Xi*q2*Y_32;
	return result;
}
double BE3D::fTX2A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1/R+0.5*alpha*q2/R3;
	return result;
}
double BE3D::fTX3A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*q*Y_11-0.5*alpha*q3*Y_32;
	return result;
}
double BE3D::fTY1A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*(COS_DELTA/R+q*Y_11*SIN_DELTA)-0.5*alpha*q*F;
	return result;
}
double BE3D::fTY2A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*Y*X_11-0.5*alpha*q*G;
	return result;
}
double BE3D::fTY3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*(D*X_11+Xi*Y_11*SIN_DELTA)+0.5*alpha*q*H;
	return result;
}
double BE3D::fTZ1A(double Xi,double Eta,double z)
{	
	double result=0.5*coef_alpha_1*(SIN_DELTA/R-q*Y_11*COS_DELTA)-0.5*alpha*q*Fp;
	return result;
}
double BE3D::fTZ2A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*D*X_11-0.5*alpha*q*Gp;
	return result;
}
double BE3D::fTZ3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*(Y*X_11+Xi*Y_11*COS_DELTA)+0.5*alpha*q*Hp;
	return result;
}


double BE3D::fTX1B(double Xi,double Eta,double z)
{
	double result=-Xi*q2*Y_32-coef_alpha_2*J_4*SIN_DELTA_2;
	return result;
}
double BE3D::fTX2B(double Xi,double Eta,double z)
{
	double result=-q2/R3-coef_alpha_2*J_5*SIN_DELTA_2;
	return result;
}
double BE3D::fTX3B(double Xi,double Eta,double z)
{
	double result=q3*Y_32-coef_alpha_2*J_6*SIN_DELTA_2;
	return result;
}

double BE3D::fTY1B(double Xi,double Eta,double z)
{
	double result=q*F-coef_alpha_2*J_1*SIN_DELTA_2;
	return result;
}
double BE3D::fTY2B(double Xi,double Eta,double z)
{
	double result=q*G-coef_alpha_2*J_2*SIN_DELTA_2;
	return result;
}
double BE3D::fTY3B(double Xi,double Eta,double z)
{
	double result=-q*H-coef_alpha_2*J_3*SIN_DELTA_2;
	return result;
}
double BE3D::fTZ1B(double Xi,double Eta,double z)
{
	double result=q*Fp+coef_alpha_2*K_3*SIN_DELTA_2;
	return result;
}
double BE3D::fTZ2B(double Xi,double Eta,double z)
{
	double result=q*Gp+coef_alpha_2*Xi*D_11*SIN_DELTA_2;
	return result;
}
double BE3D::fTZ3B(double Xi,double Eta,double z)
{
	double result=-q*Hp+coef_alpha_2*K_4*SIN_DELTA_2;
	return result;
}

double BE3D::fTX1C(double Xi,double Eta,double z)
{	
	double result=coef_alpha_1*(Xi/R3)*SIN_DELTA+Xi*q*Y_32*COS_DELTA+alpha*Xi*(3.E0*C*Eta/R5-2.E0*Z_32-Z_0);
	return result;
}
double BE3D::fTX2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*2.E0*Y_0*SIN_DELTA-D/R3+alpha*(C/R3)*(1.E0-3.E0*q2/R2);
	return result;
}
double BE3D::fTX3C(double Xi,double Eta,double z)
{
// 	double result=-coef_alpha_1*(Y/R3-Y_0*COS_DELTA)-alpha*(3.E0*C*Eta*q/R5-q*Z_0);
// 	double result=-coef_alpha_1*(Y/R3-Y_0*COS_DELTA)-alpha*(C*Eta*qR-q*Z_0);
	double result=-coef_alpha_1*YY0-alpha*(C*Eta*qR-q*Z_0);
	return result;
}
double BE3D::fTY1C(double Xi,double Eta,double z)
{	
// 	double result=coef_alpha_1*(q/R3+Y_0*SIN_DELTA*COS_DELTA)+alpha*(z*COS_DELTA/R3+3.E0*C*D*q/R5-q*Z_0*SIN_DELTA);
	double result=coef_alpha_1*(q/R3+Y_0*SIN_DELTA*COS_DELTA)+alpha*(z*COS_DELTA/R3+C*D*qR-q*Z_0*SIN_DELTA);
	return result;
}
double BE3D::fTY2C(double Xi,double Eta,double z)
{
	double result=-coef_alpha_1*2.E0*Xi*P*SIN_DELTA-Y*D*X_32+alpha*C*((Y+2*q*SIN_DELTA)*X_32-Y*q2*X_53);
	return result;
}
double BE3D::fTY3C(double Xi,double Eta,double z)
{
	double result=-coef_alpha_1*(Xi*P*COS_DELTA-X_11+Y2*X_32)+alpha*C*((D+2*q*COS_DELTA)*X_32-Y*Eta*q*X_53)+alpha*Xi*Q;
	return result;
}
double BE3D::fTZ1C(double Xi,double Eta,double z)
{
// 	double result=-Eta/R3+Y_0*COS_DELTA_2-alpha*(z*SIN_DELTA/R3-3.E0*C*Y*q/R5-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
	double result=-Eta/R3+Y_0*COS_DELTA_2-alpha*(z*SIN_DELTA/R3-C*Y*qR-Y_0*SIN_DELTA_2+q*Z_0*COS_DELTA);
	return result;
}
double BE3D::fTZ2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*2.E0*Xi*Pp*SIN_DELTA-X_11+D2*X_32-alpha*C*((D-2*q*COS_DELTA)*X_32-D*q2*X_53);
	return result;
}
double BE3D::fTZ3C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*(Xi*Pp*COS_DELTA+Y*D*X_32)+alpha*C*((Y-2*q*SIN_DELTA)*X_32+D*Eta*q*X_53)+alpha*Xi*Qp;
	return result;
}



// # Spostamento generato da una dislocazione di tipo Strike-slip

double BE3D::fSu1A(double Xi,double Eta,double z)
{
	double result=0.5*Teta+0.5*alpha*Xi*q*Y_11;
	return result;
}
double BE3D::fSu2A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*q/R;
	return result;
}
double BE3D::fSu3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Sing_iv-0.5*alpha*q2*Y_11;
	return result;
}


double BE3D::fSu1B(double Xi,double Eta,double z)
{
	double result=-Xi*q*Y_11-Teta-coef_alpha_2*I_1*SIN_DELTA;
	return result;
}
double BE3D::fSu2B(double Xi,double Eta,double z)
{
	double result=-q/R+coef_alpha_2*(Y/RpD)*SIN_DELTA;
	return result;
}
double BE3D::fSu3B(double Xi,double Eta,double z)
{
	double result=q2*Y_11-coef_alpha_2*I_2*SIN_DELTA;
	return result;
}


double BE3D::fSu1C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*Xi*Y_11*COS_DELTA-alpha*Xi*q*Z_32;
	return result;
}
double BE3D::fSu2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*(COS_DELTA/R+2.E0*q*Y_11*SIN_DELTA)-alpha*C*q/R3;
	return result;
}
double BE3D::fSu3C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*q*Y_11*COS_DELTA-alpha*(C*Eta/R3-z*Y_11+Xi2*Z_32);
	return result;
}



// # Spostamento generato da una dislocazione di tipo Dip-slip

double BE3D::fDu1A(double Xi,double Eta,double z)
{
	double result=0.5*alpha*q/R;
	return result;
}
double BE3D::fDu2A(double Xi,double Eta,double z)
{
	double result=0.5*Teta+0.5*alpha*Eta*q*X_11;
	return result;
}
double BE3D::fDu3A(double Xi,double Eta,double z)
{
	double result=0.5*coef_alpha_1*Sing_iii-0.5*alpha*q2*X_11;
	return result;
}

double BE3D::fDu1B(double Xi,double Eta,double z)
{
	double result=-q/R+coef_alpha_2*I_3*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDu2B(double Xi,double Eta,double z)
{
	double result=-Eta*q*X_11-Teta-coef_alpha_2*(Xi/RpD)*SIN_DELTA*COS_DELTA;
	return result;
}
double BE3D::fDu3B(double Xi,double Eta,double z)
{
	double result=q2*X_11+coef_alpha_2*I_4*SIN_DELTA*COS_DELTA;
	return result;
}


double BE3D::fDu1C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*COS_DELTA/R-q*Y_11*SIN_DELTA-alpha*C*q/R3;
	return result;
}
double BE3D::fDu2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*Y*X_11-alpha*C*Eta*q*X_32;
	return result;
}

double BE3D::fDu3C(double Xi,double Eta,double z)
{
	double result=-D*X_11-Xi*Y_11*SIN_DELTA-alpha*C*(X_11-q2*X_32);
	return result;
}



// # Spostamento generato da una dislocazione di tipo tensile

double BE3D::fTu1A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*Sing_iv-0.5*alpha*q2*Y_11;
	return result;
}
double BE3D::fTu2A(double Xi,double Eta,double z)
{
	double result=-0.5*coef_alpha_1*Sing_iii-0.5*alpha*q2*X_11;
	return result;
}
double BE3D::fTu3A(double Xi,double Eta,double z)
{
	double result=0.5*Teta-0.5*alpha*q*(Eta*X_11+Xi*Y_11);
	return result;
}


double BE3D::fTu1B(double Xi,double Eta,double z)
{
	double result=q2*Y_11-coef_alpha_2*I_3*SIN_DELTA_2;
	return result;
}
double BE3D::fTu2B(double Xi,double Eta,double z)
{
	double result=q2*X_11+coef_alpha_2*(Xi/RpD)*SIN_DELTA_2;
	return result;
}
double BE3D::fTu3B(double Xi,double Eta,double z)
{
	double result=q*(Eta*X_11+Xi*Y_11)-Teta-coef_alpha_2*I_4*SIN_DELTA_2;
	return result;
}


double BE3D::fTu1C(double Xi,double Eta,double z)
{
	double result=-coef_alpha_1*(SIN_DELTA/R+q*Y_11*COS_DELTA)-alpha*(z*Y_11-q2*Z_32);
	return result;
}
double BE3D::fTu2C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*2.E0*Xi*Y_11*SIN_DELTA+D*X_11-alpha*C*(X_11-q2*X_32);
	return result;
}
double BE3D::fTu3C(double Xi,double Eta,double z)
{
	double result=coef_alpha_1*(Y*X_11+Xi*Y_11*COS_DELTA)+alpha*q*(C*Eta*X_32+Xi*Z_32);
	return result;
}



void BE3D::Derivates(double burgers_vector[],struct_var FT[],double z)
{

	if(z < 0)
	{

		for(int ti=0; ti < 3; ti++)
		{

			if(burgers_vector[ti] != 0)
			{

				FT[ti].Ux_x=FT[ti].j1A-FT[ti].J1A+FT[ti].j1B+z*FT[ti].j1C;
				FT[ti].Uy_x=(FT[ti].j2A-FT[ti].J2A+FT[ti].j2B+z*FT[ti].j2C)*COS_DELTA-(FT[ti].j3A-FT[ti].J3A+FT[ti].j3B+z*FT[ti].j3C)*SIN_DELTA;
				FT[ti].Uz_x=(FT[ti].j2A-FT[ti].J2A+FT[ti].j2B-z*FT[ti].j2C)*SIN_DELTA+(FT[ti].j3A-FT[ti].J3A+FT[ti].j3B-z*FT[ti].j3C)*COS_DELTA;

				FT[ti].Ux_y=FT[ti].k1A-FT[ti].K1A+FT[ti].k1B+z*FT[ti].k1C;
				FT[ti].Uy_y=(FT[ti].k2A-FT[ti].K2A+FT[ti].k2B+z*FT[ti].k2C)*COS_DELTA-(FT[ti].k3A-FT[ti].K3A+FT[ti].k3B+z*FT[ti].k3C)*SIN_DELTA;
				FT[ti].Uz_y=(FT[ti].k2A-FT[ti].K2A+FT[ti].k2B-z*FT[ti].k2C)*SIN_DELTA+(FT[ti].k3A-FT[ti].K3A+FT[ti].k3B-z*FT[ti].k3C)*COS_DELTA;

				FT[ti].Ux_z=FT[ti].l1A+FT[ti].L1A+FT[ti].l1B+FT[ti].u1C+z*FT[ti].l1C;
				FT[ti].Uy_z=(FT[ti].l2A+FT[ti].L2A+FT[ti].l2B+FT[ti].u2C+z*FT[ti].l2C)*COS_DELTA-(FT[ti].l3A+FT[ti].L3A+FT[ti].l3B+FT[ti].u3C+z*FT[ti].l3C)*SIN_DELTA;
				FT[ti].Uz_z=(FT[ti].l2A+FT[ti].L2A+FT[ti].l2B-FT[ti].u2C-z*FT[ti].l2C)*SIN_DELTA+(FT[ti].l3A+FT[ti].L3A+FT[ti].l3B-FT[ti].u3C-z*FT[ti].l3C)*COS_DELTA;

			}

		}

	}
	else if(z == 0)
	{

		for(int ti=0; ti < 3; ti++)
		{

			if(burgers_vector[ti] != 0)
			{

				FT[ti].Ux_x=+FT[ti].j1B;
				FT[ti].Uy_x=(+FT[ti].j2B)*COS_DELTA-(+FT[ti].j3B)*SIN_DELTA;
				FT[ti].Uz_x=(+FT[ti].j2B)*SIN_DELTA+(+FT[ti].j3B)*COS_DELTA;

				FT[ti].Ux_y=+FT[ti].k1B;
				FT[ti].Uy_y=(+FT[ti].k2B)*COS_DELTA-(+FT[ti].k3B)*SIN_DELTA;
				FT[ti].Uz_y=(+FT[ti].k2B)*SIN_DELTA+(+FT[ti].k3B)*COS_DELTA;

				FT[ti].Ux_z=2*FT[ti].l1A+FT[ti].l1B+FT[ti].u1C;
				FT[ti].Uy_z=(2*FT[ti].l2A+FT[ti].l2B+FT[ti].u2C)*COS_DELTA-(2*FT[ti].l3A+FT[ti].l3B+FT[ti].u3C)*SIN_DELTA;
				FT[ti].Uz_z=(2*FT[ti].l2A+FT[ti].l2B-FT[ti].u2C)*SIN_DELTA+(2*FT[ti].l3A+FT[ti].l3B-FT[ti].u3C)*COS_DELTA;

			}

		}

	}
//	else if(CASE == 1)
//	{
//
//		for(int ti=0; ti < 3; ti++)
//		{
//
//			if(burgers_vector[ti] != 0)
//			{
//
//				FT[ti].Ux_x=-FT[ti].J1A;
//				FT[ti].Uy_x=(-FT[ti].J2A)*COS_DELTA-(-FT[ti].J3A)*SIN_DELTA;
//				FT[ti].Uz_x=(-FT[ti].J2A)*SIN_DELTA+(-FT[ti].J3A)*COS_DELTA;
//
//				FT[ti].Ux_y=-FT[ti].K1A;
//				FT[ti].Uy_y=(-FT[ti].K2A)*COS_DELTA-(-FT[ti].K3A)*SIN_DELTA;
//				FT[ti].Uz_y=(-FT[ti].K2A)*SIN_DELTA+(-FT[ti].K3A)*COS_DELTA;
//
//				FT[ti].Ux_z=+FT[ti].L1A;
//				FT[ti].Uy_z=(+FT[ti].L2A)*COS_DELTA-(+FT[ti].L3A)*SIN_DELTA;
//				FT[ti].Uz_z=(+FT[ti].L2A)*SIN_DELTA+(+FT[ti].L3A)*COS_DELTA;
//
//			}
//
//		}
//
//	}

}



void BE3D::CDisplacement(double burgers_vector[],struct_var FT[],double z)
{
  
    if(z < 0)
    {
      
      for(int ti=0; ti < 3; ti++)
      {      	      
      
    	  if(burgers_vector[ti] != 0)
    	  {
		      FT[ti].Ux=FT[ti].u1A-FT[ti].U1A+FT[ti].u1B+z*FT[ti].u1C;
    		  FT[ti].Uy=(FT[ti].u2A-FT[ti].U2A+FT[ti].u2B+z*FT[ti].u2C)*COS_DELTA-(FT[ti].u3A-FT[ti].U3A+FT[ti].u3B+z*FT[ti].u3C)*SIN_DELTA;
    		  FT[ti].Uz=(FT[ti].u2A-FT[ti].U2A+FT[ti].u2B-z*FT[ti].u2C)*SIN_DELTA+(FT[ti].u3A-FT[ti].U3A+FT[ti].u3B-z*FT[ti].u3C)*COS_DELTA;
    	  }

      }
     
   }
   else if(z == 0)
   {
     
     for(int ti=0; ti < 3; ti++)
     {     	      
     
    	 if(burgers_vector[ti] != 0)
    	 {
    		 FT[ti].Ux=FT[ti].u1B;
    		 FT[ti].Uy=FT[ti].u2B*COS_DELTA-FT[ti].u3B*SIN_DELTA;
    		 FT[ti].Uz=FT[ti].u2B*SIN_DELTA+FT[ti].u3B*COS_DELTA;
    	 }
	
     }	    
     
   }
//   else if(CASE == 1)
//   {
//
//     for(int ti=0; ti < 3; ti++)
//     {
//
//    	 if(burgers_vector[ti] != 0)
//    	 {
//    		 FT[ti].Ux=-FT[ti].U1A;
//    		 FT[ti].Uy=-FT[ti].U2A*COS_DELTA-(-FT[ti].U3A)*SIN_DELTA;
//    		 FT[ti].Uz=-FT[ti].U2A*SIN_DELTA+(-FT[ti].U3A)*COS_DELTA;
//    	 }
//
//     }
//
//   }

}



void BE3D::Strain_tensor(double Strain_components[],struct_var FF)
{

	Strain_components[0] = FF.Ux_x;
	Strain_components[1] = FF.Uy_y;
	Strain_components[2] = FF.Uz_z;
	Strain_components[3] = 0.5e0*(FF.Ux_y + FF.Uy_x);
	Strain_components[4] = 0.5e0*(FF.Ux_z + FF.Uz_x);
	Strain_components[5] = 0.5e0*(FF.Uy_z + FF.Uz_y);

}



void BE3D::Stress_tensor(double Stress_components[],struct_var FF,double mu,double lambda)
{

	double e_kk=FF.Ux_x + FF.Uy_y + FF.Uz_z;

	Stress_components[0]=coefpi*(lambda*e_kk + 2.E0*mu*FF.Ux_x);
	Stress_components[1]=coefpi*(lambda*e_kk + 2.E0*mu*FF.Uy_y);
	Stress_components[2]=coefpi*(lambda*e_kk + 2.E0*mu*FF.Uz_z);
	Stress_components[3]=coefpi*mu*(FF.Ux_y + FF.Uy_x);
	Stress_components[4]=coefpi*mu*(FF.Ux_z + FF.Uz_x);
	Stress_components[5]=coefpi*mu*(FF.Uy_z + FF.Uz_y);

}





void BE3D::SHIFT(double Dx,double Dy,double Dz)
{
	PARAM.GEOM.cc[0] = PARAM.GEOM.cc[0] + Dx;
	PARAM.GEOM.cc[1] = PARAM.GEOM.cc[1] + Dy;
	PARAM.GEOM.cc[2] = PARAM.GEOM.cc[2] + Dz;

	PARAM.GEOM.pc[0] = PARAM.GEOM.pc[0] + Dx;
	PARAM.GEOM.pc[1] = PARAM.GEOM.pc[1] + Dy;
	PARAM.GEOM.pc[2] = PARAM.GEOM.pc[2] + Dz;
}




