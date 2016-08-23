/*
 * Main.cpp
 *
 *  Created on: 04/lug/2014
 *      Author: dino
 */

#include <iostream>
#include <string>

using namespace std;

#include "Halfspace.h"


void example_1();
void example_2();
void example_3();



int main()
{
	example_2();

	return 0;
}





void example_1()
{
	CONSOLE CONSOLE_out;                                    // Creation of output objects:
	DATAFILE DATA;                     						// CONSOLE_out, for console output (CONSOLE type)
	DATAFILE DATA_sym(0e0,1000,10e3);                 		// DATA and DATA_sym, for output on datafile (DATAFILE type)


	double mu=4e9,nu=0.25;                                  // Creation of crust (MEDIUM object)
	MEDIUM crust(mu,nu);

	double x0=1e3,y0=2e3,z0=-3e3,phi=60,theta=0;
    double DP=690e6;

	double R1=5e2;
	MOGI my_sphere("MY SPHERE",x0,y0,z0,R1);						// Creation of my_sphere (MOGI object)

	double a=6e2,b=3e2;
	YANG my_spheroid("MY SPHEROID",x0,y0,z0,a,b,phi,theta);			// Creation of my_spheroid (YANG object)

	double R2=5e2;
	FIALKO my_sill("MY PENNY",x0,y0,z0,R2);							// Creation of my_sill (FIALKO object)

	int input=-1;
	while(input != 0)
	{
		cout << endl << endl << "Choose source type (1- MOGI,2- YANG, 3- FIALKO, 4- EXIT)" << endl;
		cin >> input;

		switch(input)												// MENU
		{															// Selection of the source type
		case (1):                                                	// and generation of the corresponding output
			my_sphere.SET(crust,DP);                         		//
			my_sphere.PRINT(CONSOLE_out);                    		// 1- MOGI object
			my_sphere.PRINT(DATA);                       		    //
			my_sphere.PRINT(DATA_sym);                       		//
		break;                                                  	//
		case (2):                                                	//
			my_spheroid.SET(crust,DP);                       		// 2- YANG object
			my_spheroid.PRINT(CONSOLE_out);                  		//
			my_spheroid.PRINT(DATA);                     		    //
			my_spheroid.PRINT(DATA_sym);                     		//
		break;                                                   	//
		case (3):													// 3- FIALKO object
			my_sill.SET(crust,DP);									//
			my_sill.PRINT(CONSOLE_out);								//
			my_sill.PRINT(DATA);								    //
			my_sill.PRINT(DATA_sym);								//
			break;
		case (4):
			exit(0);

		default:
			input=-1;
		}
	}
}






void example_2()
{
	CONSOLE CONSOLE_out;										// CONSOLE_out, for console output (CONSOLE type)
	DATAFILE DATA(0);											// DATA, for output on datafile (DATAFILE type)
	DATAFILE DATA_hires(0,15e3,15e3,200,200,"hires");			// DATA, for output on datafile (DATAFILE type)

	double mu=4e9,nu=0.25;
	MEDIUM crust(mu,nu);										// Creation of crust (MEDIUM object)

	double x0=0,y0=0,z0=-5e2,phi=0,theta=0,DP=10e6;

	int NDP;													//
	cout << "Insert discretization parameter" << endl;			// Assignment of NDP (discretization parameter)
	cin >> NDP;													//

	SOURCE_BE3D *p_source;

	double L=5e3,W=3e3;
	FAULT my_fault("my fault",x0,y0,z0,L,W,phi,theta,NDP);		// Creation of my_fault (FAULT object)
	my_fault.PRINT();											// PRINT() generates a graphic of my_fault model geometry

	double a=2.5e3,b=1.5e3;
	SILL my_sill("my sill",x0,y0,z0,a,b,phi,theta,NDP);			// Creation of my_sill (SILL object)
	my_sill.PRINT();											// PRINT() generates a graphic of my_sill model geometry

	double h1=1e2,h2=1e2,h3=1e2;
	PARAL my_paral("my paral",x0,y0,z0,h1,h2,h3,phi,NDP);		// Creation of my_paral (PARAL object)
	my_paral.PRINT();											// PRINT() generates a graphic of my_paral model geometry


	int input=-1;
	while(input != 0){
		cout << endl << endl << "Choose source type (1- FAULT,2- SILL, 3- PARAL, 4- EXIT)" << endl;
		cin >> input;

		switch(input){											// 				MENU
		case (1):												// Selection of the source type
				p_source = &my_fault;							// 1- FAULT object
			break;												//
		case (2):												//
				p_source = &my_sill;							// 2- SILL object
			break;												//
		case (3):												//
				p_source = &my_paral;							// 3- PARAL object
			break;												//
		case (4):												//
			exit(0);											// 4- Exit
		default:
			input=-1;
		}

		if(input != -1){
			p_source->SET(crust,2,0.,0.,DP);					//
			p_source->SOLVE();									//
			p_source->PRINT();									// Generation of the output for the selected source
			p_source->PRINT(CONSOLE_out);						//
			p_source->PRINT(DATA_hires);						//
		}														//
	}
}







#include "NOT_PLANAR.h"

void example_3()
{

	CONSOLE CONSOLE_out;												// CONSOLE_out, for console output (CONSOLE type)
	DATAFILE DATA;														// DATA, for output on datafile (DATAFILE type)

	double mu=4e9,nu=0.25;
	MEDIUM crust(mu,nu);												// Creation of crust (MEDIUM object)

	double x0=0,y0=0,z0=-5e3,phi=30;
	double DP=10e6;

	int NDP=10;
	double L=3e3,W_1=4e3,W_2=2e3,theta_1=70,theta_2=20;
	NOT_PLANAR my_source("my source",x0,y0,z0,L,W_1,W_2,phi,theta_1,theta_2,NDP);  // NOT_PLANAR object creation

	my_source.PRINT();													// PRINT() generates a graphic of my_fault model geometry
	my_source.SET(crust,2,0.,0.,DP);
	my_source.SOLVE();
	my_source.PRINT(CONSOLE_out);										// Information about my_source on standard output
	my_source.PRINT(DATA);												// Generation of DATAFILE output


	double y_range=10e3,z_range=10e3,y_0=y_range*0.5;					//
	int N_y=300,N_z=300;												// Parameter for the vertical grid of observation points
	double step_y=y_range/N_y,step_z=z_range/N_z;						//

	ofstream out_file;
	out_file.open("Vertical_plane.dat");								// Output file opening

	double U[3],Strain[6],e_kk;
	double y,x_obs,y_obs,z,phi_rad=acos(-1e0)*phi/180e0;

	for(int ti=0; ti < N_y; ti++)															//
	{																						//
		y=step_y*ti-y_0;																	//
																							//
		x_obs =  sin(phi_rad)*y;															//
		y_obs =  cos(phi_rad)*y;															//
																							//
		for(int tj=0; tj < N_z; tj++)														//
		{																					//
			z=-step_z*tj;																	//
																							//
			my_source.DISPLACEMENT(x_obs,y_obs,z,U);										// Grid generation
			my_source.STRAIN(x_obs,y_obs,z,Strain);											//
																							//
			e_kk=Strain[0]+Strain[1]+Strain[2];												//
																							//
			out_file << y << " " << z << " " << U[2] << " " << e_kk << endl;				//
		}																					//
		out_file << endl;																	//
	}																						//
	out_file.close();																		//
}












