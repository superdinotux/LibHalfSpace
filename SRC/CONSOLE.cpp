/*
 * OUTPUT_class.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: dino
 */

#include <iostream>
#include <fstream>
#include <string>


using namespace std;


#include "CONSOLE.h"


CONSOLE::CONSOLE()
{
	int_flag = 0;

	out_stream = &cout;

	flag_print = 0;

}



CONSOLE::CONSOLE(ostream &out_stream_i)
{
	int_flag = 0;

	flag_print = 0;

	out_stream = &out_stream_i;
}



CONSOLE::CONSOLE(int int_flag_i)
{
	int_flag = int_flag_i;

	flag_print = 0;

	out_stream =&cout;
}



CONSOLE::CONSOLE(int int_flag_i,ostream &out_stream_i)
{
	int_flag = int_flag_i;

	flag_print = 0;

	out_stream = &out_stream_i;
}




CONSOLE::CONSOLE(int int_flag_i,ostream *out_stream_i,int flag_print_i)
{
	int_flag = int_flag_i;

	out_stream = out_stream_i;

	flag_print = flag_print_i;
}




