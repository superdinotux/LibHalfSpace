#ifndef FOO_H 
#define FOO_H 


const int max_dim=27000;


struct source_param
{
  
    double Delta_P; // overpressure
    double c1;	// source depth
    double h_1; // x dimension
    double h_2; // y dimension
    double h_3; // z dimension
    double max_side; // z dimension
    int type; // 1-TENSILE (T), 2-TENSILE, STRIKE-SLIP, DIP-SLIP (TSD)
    
};  



struct medium_param
{
  
    double mu;	// rigidity
    double lambda;
    double nu; // poisson ratio
    double K;
  
};  



struct discretization_param
{
  
    int NUM;
    double cc[max_dim][3];
    double pc[max_dim][3];
    double n_matrix[max_dim][3];
    double vettore_delta[max_dim];
    double LV[max_dim];
    double WV[max_dim];
    int flag[max_dim];
    int N_1;
    int N_2;
    int N_3;
    int TOTALE;
    int NUM_COND;    
    
};



struct obs_point_ti
{
    
    double xn;
    double yn;
    double z;  
    double n[3];
    int ti;    
    int I1;
    
};



struct subsource_param_tj
{
    
    double x_cc;
    double y_cc;
    double c_cc;
    double x_pc;
    double y_pc;
    double c_pc; 
    double delta_gradi;
    double L;
    double W;    
    int flag;
    int tj;      
    int I2;
    
};



struct struct_var
{

	double j1A,J1A,j1B,j1C;
	double j2A,J2A,j2B,j2C;
	double j3A,J3A,j3B,j3C;
	
	double k1A,K1A,k1B,k1C;	
	double k2A,K2A,k2B,k2C;
	double k3A,K3A,k3B,k3C;
	
	double l1A,L1A,l1B,u1C,l1C;
	double l2A,L2A,l2B,u2C,l2C;
	double l3A,L3A,l3B,u3C,l3C;
	
	double Ux_x,Uy_y,Uz_z;
	double Ux_y,Uy_x;
	double Ux_z,Uz_x;
	double Uy_z,Uz_y;  
	
	double u1A,U1A,u1B;
	double u2A,U2A,u2B;
	double u3A,U3A,u3B;
	
	double Ux,Uy,Uz;	
    
};


struct gravity_grid_param
{
  
	double vett_TOT1[4];  
	double vett_TOT2[4];  
	double vett_TOT3[4];  
	double vett_TOT4[4];  

};




#endif

