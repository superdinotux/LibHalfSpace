/*
 * OKADA.h
 *
 *  Created on: 24/ott/2014
 *      Author: dino
 */

#ifndef OKADA_H_
#define OKADA_H_


///*! \brief Brief description.
// *         Brief description continued.
//
// */


#include "string.h"


#include "BE3D.h"
#include "SOURCE_3D.h"


class OKADA: public SOURCE_3D{

private:	// PRIVATE MEMBERS of the CLASS

	double L,W;									//
	double theta_deg,phi_deg;					//
												//	PRIVATE DATA MEMBERS
	BE3D BE;									//
	int FLAG_COMP;								//


	void NEW(double L_i,double W_i,double theta_deg_i,double phi_deg_i);

	void PRINT_GEOM_PROP(ostream *out_stream);	//  PRIVATE METHOD MEMBER

	void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print=0);

public:		// All PUBLIC MEMBERS of the CLASS are METHODS

	OKADA();

	OKADA(string label,double x0_i,double y0_i,double z0_i,double L,double W,double theta_deg,double phi_deg)
		:SOURCE_3D(label,x0_i,y0_i,z0_i){NEW(L,W,theta_deg,phi_deg);};

	virtual ~OKADA();

	/*! Association of the medium param with the source buried in the Half-space.*/
	/*!
	 *  We also fix:
	 * - the strike slip component (Us)
	 * - the dip slip component (Ud)
	 * - the tensile component (Ut)
	 * The value assigned to FLAG_COMP assigns the Okada implementation adopted:
	 * - 1, Ferrari et Al. implementation
	 * - 2, Okada implementation (dc3d.f).
	 * If the value is not assigned, the default value is 2.
	 */

	void SET(MEDIUM MEDIUM_PAR_i,double U_s,double U_d,double U_t,int FLAG_COMP = 2);

	/*! Computation of the displacement components (Ux=U[0], Uy=U[1], Uz=U[2]) at the observation point (x,y,z) */

	void DISPLACEMENT(double x,double y,double z,double U[]);

	/*! Computation of the strain tensor components
	 * (Sxx=Strain[0], Syy=Strain[1], Szz=Strain[2],Sxy=Strain[3],Sxz=Strain[4], Syz=Strain[5])
	 *  at the observation point (x,y,z) */

	void STRAIN(double x,double y,double z,double Strain[]);

	/*! Computation of the stress tensor components
	 * (Sxx=Stress0], Syy=Stress[1], Szz=Stress[2],Sxy=Stress[3],Sxz=Stress[4], Syz=Stress[5])
	 *  at the observation point (x,y,z) */

	void STRESS(double x,double y,double z,double Stress[]);

	/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

	void PRINT();

	/*! Print to standard output useful information about the source */

	void PRINT(CONSOLE &out);

	/*! Generation of a datafile and a gnuplot script.
	 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

	void PRINT(DATAFILE &out);


};

#endif /* OKADA_H_ */





