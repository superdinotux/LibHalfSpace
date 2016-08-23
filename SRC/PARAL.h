/*
 * PARAL.h
 *
 *  Created on: 07/lug/2014
 *      Author: dino
 */

#ifndef PARAL_H_
#define PARAL_H_


#include "cmath"


#include "BE3D.h"
#include "CONSOLE.h"
#include "DATAFILE.h"
#include "SOURCE_BE3D.h"


extern "C" {
void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
}



class PARAL: public SOURCE_BE3D {

	double h_1;
	double h_2;
	double h_3;

	void NEW(double h1,double h2,double h3);

	void PRINT_GEOM_PROP(ostream *out_stream);
    void PRINT_STRESS_PROP(ostream *stream);
    void PRINT_BV_COMP(ostream *stream);

    void PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print = 0);


public:

	PARAL(string s_label,double x0,double y0,double z0,double h1,double h2,double h3,double fi,int NUM)
			:SOURCE_BE3D(s_label,x0,y0,z0,fi,NUM)
			{NEW(h1,h2,h3);};

	virtual ~PARAL(){};

	/*! Association of the medium param with the source buried in the Half-space
	 * This method uses the integer variable option_flag to choose the type of input expected:
	 * - option_flag = 1, we consider the direct problem in fracture theory; a uniform slip is assigned over the source boundaries, with strike-slip component equal to INPUT_1,  with dip-slip component equal to INPUT_2 and with tensile component  with strike-slip component equal to INPUT_3;
	 * - option_flag = 2, we consider the inverse problem in fracture theory, uniform boundary conditions are assigned over source boundaries.
	 * The optional variable FLAG_COMP_i can be specified in order to select the Okada implementation adopted:
	 * - FLAG_COMP_i = 1, the implementation given by Ferrari et Al. (default option);
	 * - FLAG_COMP_i = 2, the original implementation by Okada (dc3d.f code)*/

	void SET(MEDIUM MEDIUM_PAR_i,int option_flag,double INPUT_1,double INPUT_2,double INPUT_3,int FLAG_COMP_i = 1);

	/*! SOLVE implements the discontinuity element method (BEM). The implementation of the method for the PARAL model source is designed to solve an issue related to closed sources. */

	void SOLVE();

	/*! Generation of a datafile and a gnuplot script suitable to obtain a graphic representation of the source model geometry */

    void PRINT(){SOURCE_BE3D::PRINT();};

	/*! Print to standard output useful information about the source */

    void PRINT(CONSOLE &out);

	/*! Generation of a datafile and a gnuplot script.
	 * According to the DATAFILE object passed to PRINT, we obtain displacement and/or strain and/or stress maps at the free surface */

    void PRINT(DATAFILE &out);

};




#endif /* PARAL_H_ */











//    PARAL operator=(PARAL &obj)
//    {
//
//    	cout << "PARAL - Overload dell'operatore =" << endl;
//
//    	if(obj.NBE != NBE)
//    	{
//    		delete BE;
//
//    		NBE=obj.NBE;
//    		BE = new BE3D [NBE];
//
//    		if(!BE)
//    		{
//    			cout << "Impossibile allocare." << endl;
//    //			exit(1);
//    		}
//    	}
//
//    	for(int i=0; i < NBE; i++)
//    	{
//    		BE[i] = obj.BE[i];
//    	}
//
//    	return *this;
//
//    };





