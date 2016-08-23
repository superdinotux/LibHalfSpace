/*
 * DATAFILE_class.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: dino
 */

#include <iostream>
#include <fstream>
#include <string>


using namespace std;


#include "DATAFILE.h"



DATAFILE::DATAFILE(double ri_i,int Nr_i,double rf_i,string extra)
{

	int_flag = 1;

	flag_datatype = 0;

	if(extra != "")
	{
		extra = "-" + extra;
	}

	filename = "-DISPL" + extra;

	Nr = Nr_i;

	ri = ri_i;
	rf = rf_i;
}



DATAFILE::DATAFILE(string extra)
{
	flag_datatype = 1;

	if(extra != "")
	{
		extra = "-" + extra;
	}

	int_flag = 0;

	filename = "_MAPS" + extra;

	N_1=100;
	N_2=100;

	DX1=10.e3;
	DX2=10.e3;
}






DATAFILE::DATAFILE(int int_flag_i,string extra)
{
	flag_datatype = 1;

	if(extra != "")
	{
		extra = "-" + extra;
	}

	int_flag = int_flag_i;

	switch(int_flag)
	{
	    case (0):
		filename = "_MAPS" + extra;
		break;
	    case (1):
		filename = "_MAP" + extra;
		break;
	    case (2):
		filename = "_MAP" + extra;
		break;
	    case (3):
		filename = "_MAP" + extra;
		break;

	    default:
		cout << "Error" << endl;
//		exit(1);
	}

	N_1=100;
	N_2=100;

	DX1=10.e3;
	DX2=10.e3;
}





DATAFILE::DATAFILE(double DX1_i,double DX2_i,int N_1_i,int N_2_i,string extra)
{
	flag_datatype = 1;

	string suffix = "-" + extra;


	int_flag = 0;

	filename = "_MAPS" + suffix;

	N_1=N_1_i;
	N_2=N_2_i;

	DX1=DX1_i;
	DX2=DX2_i;
}






DATAFILE::DATAFILE(int int_flag_i,double DX1_i,double DX2_i,int N_1_i,int N_2_i,string extra)
{
	flag_datatype = 1;

	string suffix = "-" + extra;


	int_flag = int_flag_i;

	switch(int_flag)
	{
	    case (0):
		filename = "_MAPS" + suffix;
		break;
	    case (1):
		filename = "_MAP" + suffix;
		break;
	    case (2):
		filename = "_MAP" + suffix;
		break;
	    case (3):
		filename = "_MAP" + suffix;
		break;
		
	    default:
		cout << "Error" << endl;
//		exit(1);
	}

	N_1=N_1_i;
	N_2=N_2_i;

	DX1=DX1_i;
	DX2=DX2_i;
}









void DATAFILE::GET_flag_datatype(int &flag_datatype_o)
{
	flag_datatype_o = flag_datatype;
}



void DATAFILE::GET_file_parameters(string &filename_o,int &Nr_o,double &ri_o,double &rf_o)
{
	filename_o = filename;

	Nr_o = Nr;
	ri_o = ri;
	rf_o = rf;
}



void DATAFILE::GET_file_parameters(string &filename_o,int &N_1_o,int &N_2_o,double &DX1_o,double &DX2_o)
{
	filename_o = filename;

	N_1_o=N_1;
	N_2_o=N_2;

	DX1_o=DX1;
	DX2_o=DX2;
}






