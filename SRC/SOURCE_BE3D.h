/*
 * SOURCE_BE3D.h
 *
 *  Created on: 18/lug/2014
 *      Author: dino
 */

#ifndef SOURCE_BE3D_H_
#define SOURCE_BE3D_H_


#include "BE3D.h"
#include "GRID.h"
#include "SOURCE.h"


class SOURCE_BE3D: public SOURCE {

protected:	// PROTECTED MEMBERS of the CLASS

	double phi_degree;								//
													//
	int NBE;										//
	BE3D *BE;										//	PROTECTED DATA MEMBERS
													//
	int NDP,FLAG_COMP;								//
	GRID GRID_PAR;									//


	void PRINT_GEOM_PROP(ostream *out_stream);		//
	void PRINT_STRESS_PROP(ostream *stream);		//	PROTECTED METHOD MEMBERS
	void PRINT_BV_COMP(ostream *stream);			//

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print = 0);


public:		// All PUBLIC MEMBERS of the CLASS are METHODS

	SOURCE_BE3D();

	SOURCE_BE3D(string s_label,double x0,double y0,double z0,double phi_angle_i,int N)
				:SOURCE(s_label,x0,y0,z0)
				{phi_degree=phi_angle_i; NDP=N;};

	virtual ~SOURCE_BE3D();

	/*!Association of the medium param with the source buried in the Half-space
	 * This method uses the integer variable option_flag to choose the type of input expected:
	 * - option_flag = 1, we consider the direct problem in fracture theory; a uniform slip is assigned over the source boundaries, with strike-slip component equal to INPUT_1,  with dip-slip component equal to INPUT_2 and with tensile component  with strike-slip component equal to INPUT_3;
	 * - option_flag = 2, we consider the inverse problem in fracture theory, uniform boundary conditions are assigned over source boundaries.
	 * The optional variable FLAG_COMP_i can be specified in order to select the Okada implementation adopted:
	 * - FLAG_COMP_i = 1, the implementation given by Ferrari et Al. (default option);
	 * - FLAG_COMP_i = 2, the original implementation by Okada (dc3d.f code)*/

	virtual void SET(MEDIUM MEDIUM_PAR_i,int option_flag,double INPUT_1,double INPUT_2,double INPUT_3,int FLAG_COMP_i = 1)
	{
		MEDIUM_PAR = MEDIUM_PAR_i;

		switch (option_flag)
		{
			case (1):

					for(int i=0; i < NBE; i++)
					{
						BE[i].PUT_BV(INPUT_1,INPUT_2,INPUT_3);
					}

					break;

			case (2):

					for(int i=0; i < NBE; i++)
					{
						BE[i].PUT_BC(INPUT_1,INPUT_2,INPUT_3);
					}

					break;

			default:

				cout << "Option " << option_flag << " is not available" << endl;
		}

		FLAG_COMP = FLAG_COMP_i;
	};

/*! Computation of the displacement components (Ux=U[0], Uy=U[1], Uz=U[2]) at the observation point (x,y,z) */

	void DISPLACEMENT(double x,double y,double z,double U[]);

/*! Computation of the strain tensor components
 *  (Sxx=Strain[0], Syy=Strain[1], Szz=Strain[2],Sxy=Strain[3],Sxz=Strain[4], Syz=Strain[5])
 *  at the observation point (x,y,z) */

	void STRAIN(double x,double y,double z,double Strain[]);

/*! Computation of the stress tensor components
 *  (Sxx=Stress0], Syy=Stress[1], Szz=Stress[2],Sxy=Stress[3],Sxz=Stress[4], Syz=Stress[5])
 *  at the observation point (x,y,z) */

	void STRESS(double x,double y,double z,double Stress[]);

	/**< SOLVE is the method used to calculate the solution relative to the inverse problem in fracture theory.	SOLVE implements the discontinuity element method (BEM) */

    virtual void SOLVE();

/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

	virtual void PRINT();

/*! Print to standard output useful information about the source */

	virtual void PRINT(CONSOLE &out);

/*! Generation of a datafile and a gnuplot script.
 *  According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface. */

	virtual void PRINT(DATAFILE &out);

/*! With this method, we query the object in order to retrieve the phi angle associated to the source. */

	void GET_phi(double phi_o){phi_o = phi_degree;};

/*! With this method, we query the object in order to retrieve the NBE parameter associated to the source. */

	void GET_NBE(int NBE_o){NBE_o = NBE;};

/*! With this method, we query the object in order to retrieve the NBE parameter associated to the source. */

	void GET_NDP(int NDP_o){NDP_o = NDP;};
	void GET_FLAG_COMP(int FLAG_COMP_o){FLAG_COMP_o = FLAG_COMP;};

/*! We change the coordinates of the source center (x0,y0,z0) */

	void PUT_x0_y0_z0(double x0_i, double y0_i, double z0_i);

};


#endif






