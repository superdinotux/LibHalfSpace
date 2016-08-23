/*
 * CONSOLE_class.h
 *
 *  Created on: 13/mag/2014
 *      Author: dino
 */

#ifndef CONSOLE_CLASS_H_
#define CONSOLE_CLASS_H_


/*! This class is designed to configure the output for the model source objects.
 */



#include <fstream>
#include <string>


using namespace std;


#include "OUTPUT.h"

class CONSOLE: public OUTPUT
{
	int flag_print;

public:

	ostream *out_stream;

/*! This method generates a CONSOLE object used to print to the standard output all the available information about the source  */

	CONSOLE();

/*! This method has the same output of CONSOLE() but this output is redirected to the text file associated to out_stream.   */

	CONSOLE(ostream &out_stream_i);

/*! This method generates a CONSOLE object used to print to the standard output some selected information about the source: */
/*! 0  - same output of CONSOLE() method */
/*! 1  - geometric properties */
/*! 2  - stress properties */
/*! 3  - Burgers vector components */
/*!	10 - main features */
/*! According to the type of the source, some of the option are not available and in this case a warning message will be printed by the library. */

	CONSOLE(int int_flag_i);

/*! This method has the same behavior of CONSOLE(int int_flag_i) but this output is redirected to the text file associated to out_stream.  */

	CONSOLE(int int_flag_i,ostream &out_stream_i);


	CONSOLE(int int_flag_i,ostream *out_stream_i,int flag_print = 0);

	void GET_flag_print(int &flag_print_o){flag_print_o = flag_print;};


};



#endif /* CONSOLE_CLASS_H_ */
