/*
 * DATAFILE_class.h
 *
 *  Created on: 13/mag/2014
 *      Author: dino
 */

#ifndef DATAFILE_CLASS_H_
#define DATAFILE_CLASS_H_

#include <fstream>
#include <string>
#include <vector>


using namespace std;


#include "OUTPUT.h"


class DATAFILE: public OUTPUT
{
	int flag_datatype;

	string filename;

	int N_1,N_2;
	double DX1,DX2;
	double X01,X02;

	string SUFFIX;

	int Nr;
	double ri;
	double rf;
	vector<double> a_r;

public:

/*! This method generates a DATAFILE object usable only by symmetric sources. When PRINT is invoked with this DATAFILE object as argument, then the radial and the vertical components of the displacement field at the free surface are computed for the selected source. */
/*! If PRINT is invoked by a asymmetric source, then a warning message is printed by the library */

	DATAFILE(double ri_i,int Nr_i,double rf_i,string extra = "");

/*! This method generates a DATAFILE object which builds a grid of observation points with these features: */
/*! - equal width on both x and y axes   (DX1 = 10 km, DX2 = 10 km) */
/*! - equal density on both x and y axes (N_1 = 100, N_2 = 100) */
/*! When print is invoked with this DATAFILE object as argument, then displacement, strain and stress fields produced by the source are printed.	*/
/*! The optional string extra can be used to modify the default name of the files created by the PRINT method. */

	DATAFILE(string extra = "");

/* This method has the same behavior of DATAFILE(string extra = ""), but with int_flag we select the field that must be printed: */
/*! 0  - all */
/*! 1  - displacement */
/*! 2  - strain */
/*! 3  - stress */

	DATAFILE(int int_flag_i,string extra = "");

/*! This method is used to build DATAFILE object which defines a non standard grid of observation points.*/
/*!	Displacement, strain and stress fields are selected to be printed and in this case the string must be specified in order to distinguish the non standard grid. */

	DATAFILE(double DX1_i,double DX2_i,int N_1_i,int N_2_i,string extra);

/*! Same behavior of the previous method but int_flag allows to select the desired output: */
/*! 0  - all */
/*! 1  - displacement */
/*! 2  - strain */
/*! 3  - stress */

	DATAFILE(int int_flag_i,double DX1_i,double DX2_i,int N_1_i,int N_2_i,string extra);



	void GET_flag_datatype(int &flag_datatype_o);

	void GET_file_parameters(string &filename_o,int &Nr_o,double &ri_o,double &rf_o);

	void GET_file_parameters(string &filename_o,int &N_1_o,int &N_2_o,double &DX1_o,double &DX2_o);

};



#endif /* DATAFILE_CLASS_H_ */
