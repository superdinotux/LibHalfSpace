// PRINT_TEMPLATE.h

#ifndef PRINT_TEMPLATE_H_
#define PRINT_TEMPLATE_H_


#include <string>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <typeinfo>
#include <ctime>

#include <typeinfo>


using namespace std;


#include "Utilities.h"
#include "CONSOLE.h"


template <class TSOURCE>
void print_map_displ(TSOURCE source,string filename,int N_1,int N_2,double DX1,double DX2)
{
	int type=0;
	print_displ_script(source,filename,DX1,DX2,type);

	string label;
	source.GET_label(label);

	filename = label + filename + "-DISPL.dat";

	double step_X1,step_X2;
	double x,y,z=0;
	double U[3];


	//        !-------------------------------------------------------------------!
	//        ! OUTPUT FILE OPENING

	ofstream fid_out;
	fid_out.open(filename.c_str());

	ostream *stream;
	stream = &fid_out;

	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->setf(ios::showpos);
	stream->precision(5);


	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	char* dt = ctime(&now);


	*stream << "# Datafile creation: " << dt << endl;

	*stream << endl;

	*stream << "# Parameters of the datafile: " << endl;
	*stream << "# " << "N_1 = " << N_1 << " ,N_2 = " << N_2 << endl;
	*stream << "# " << "DX_1 = " << DX2 << " m, DX_2 = " << DX2 << " m" << endl;

	*stream << endl;


	int input=10;
	int flag_print=1;

	CONSOLE CONSOLE_out(input,stream,flag_print);
	source.PRINT(CONSOLE_out);


	//        !-------------------------------------------------------------------!
	//        ! SET GRID OF OBSERVATION POINTS

	step_X1=DX1/N_1;
	step_X2=DX2/N_2;

	double X_01=DX1*0.5;
	double X_02=DX2*0.5;

	int cont=0;
	int np=N_1*N_2;
	int c[10]={0,0,0,0,0,0,0,0,0,0};

	cout << endl << "Generation of the grid about surface displacement" << endl;

	*stream << endl;
	*stream << "#   " <<  "x (m)" << "         " << "y (m)" << "            " << "Ux (m)" << "          " <<  "Uy (m)" << "           " <<  "Uz (m)" << endl;
	*stream << endl;

	for(int i=0; i <= N_1; i++)
	{

		x=step_X1*i-X_01;

		for(int j=0; j <= N_2; j++)
		{

			y=step_X2*j-X_02;

			cont = i * N_1 + j + 1;


			source.DISPLACEMENT(x,y,z,U);

			//     	!-------------------------------------------------------------------
			//     	! WRITE OUTPUT

			stream->fill(' ');
			*stream << right << setw(8) << x;

			stream->fill(' ');
			*stream << setw(2) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << y ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << U[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << U[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << U[2] ;
			*stream << endl;

		}

		fid_out << endl;

		Cycle_counter(cont,np,c);

	}

	fid_out.close();

}




template <class TSOURCE>
void print_displ_script(TSOURCE source,string filename,double DX1,double DX2,int type)
{
	string source_label;
	source.GET_label(source_label);

	string PLOT_DISPL = source_label + filename + "-DISPL.gnuplot";

	string IMAGE_NAME_Ux = source_label + filename + "-Ux.eps";

	string IMAGE_NAME_Uy = source_label + filename + "-Uy.eps";

	string IMAGE_NAME_Uz = source_label + filename + "-Uz.eps";

	if(type == 0)
	{
		filename = source_label + filename + "-DISPL.dat";
	}
	else if(type != 0)
	{
		filename = source_label + filename + ".dat";
	}


	int i_pDX1 = (int) DX1*0.5*1e-3;
	int i_mDX1 = -i_pDX1;

	int i_pDX2 = (int) DX2*0.5*1e-3;
	int i_mDX2 = -i_pDX2;


	ofstream PLOT_SOURCE_DISPL;
	PLOT_SOURCE_DISPL.open(PLOT_DISPL.c_str());

	PLOT_SOURCE_DISPL << "set term postscript eps enhanced														" << endl;

	PLOT_SOURCE_DISPL << "set size 0.7,1																		" << endl;

	PLOT_SOURCE_DISPL << "set xlabel 'x (km)'																	" << endl;
	PLOT_SOURCE_DISPL << "set ylabel 'y (km)'																	" << endl;

	PLOT_SOURCE_DISPL << "set xtics 1																			" << endl;
	PLOT_SOURCE_DISPL << "set ytics 1																			" << endl;

	PLOT_SOURCE_DISPL << "set view map																			" << endl;
	PLOT_SOURCE_DISPL << "set pm3d																				" << endl;
	PLOT_SOURCE_DISPL << "set cbrange [-2:2]																	" << endl;
	PLOT_SOURCE_DISPL << "set cblabel '(m)' rotate by 0 offset -2.5,-12.5										" << endl;
	PLOT_SOURCE_DISPL << "set palette defined" <<
						 "( 0 '#0d0375',  4 '#1e158c',  8 '#3228a3', 12 '#4a40b7', 16 '#655cc9', 20 '#7C71F6',  " <<
						 " 24 '#7787ff', 28 '#77a7ff', 32 '#77c6ff', 36 '#77e6ff', 38 '#9bf5ff', 40 '#b2fdff',  " <<
						 " 42 '#dcffff', 44 '#ffffff', 46 '#ffffdc', 48 '#fffdb2', 50 '#FFF59B', 52 '#ffe677',  " <<
						 " 56 '#ffc677', 60 '#ffa777', 64 '#ff8777', 68 '#ef614f', 72 '#ed3017', 76 '#dd0404',  " <<
						 " 80 '#c40d19', 84 '#b50e19', 88 '#ad202a')											" << endl;

	PLOT_SOURCE_DISPL << "set nokey																				" << endl;
	PLOT_SOURCE_DISPL << "set style increment userstyles 														" << endl;
	PLOT_SOURCE_DISPL << "set contour																			" << endl;
	PLOT_SOURCE_DISPL << "set cntrparam levels auto																" << endl;
	PLOT_SOURCE_DISPL << "set nosurface																			" << endl;

	PLOT_SOURCE_DISPL << "set xrange [" << i_mDX1 << ":" << i_pDX1 << "]										" << endl;
	PLOT_SOURCE_DISPL << "set yrange [" << i_mDX2 << ":" << i_pDX2 << "]										" << endl;

	PLOT_SOURCE_DISPL 																							  << endl;
	PLOT_SOURCE_DISPL << "set output '" << IMAGE_NAME_Ux <<"'													" << endl;
	PLOT_SOURCE_DISPL << "set title  '" << source_label << " - U_{x}'											" << endl;
	PLOT_SOURCE_DISPL << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):3 with lines						" << endl;
	PLOT_SOURCE_DISPL 																							  << endl;
	PLOT_SOURCE_DISPL << "set output '" << IMAGE_NAME_Uy <<"'													" << endl;
	PLOT_SOURCE_DISPL << "set title  '" << source_label << " - U_{y}'											" << endl;
	PLOT_SOURCE_DISPL << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):4 with lines						" << endl;
	PLOT_SOURCE_DISPL 																							  << endl;
	PLOT_SOURCE_DISPL << "set output '" << IMAGE_NAME_Uz <<"'													" << endl;
	PLOT_SOURCE_DISPL << "set title  '" << source_label << " - U_{z}'											" << endl;
	PLOT_SOURCE_DISPL << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):5 with lines   					" << endl;

}




template <class TSOURCE>
void print_map_strain(TSOURCE source,string filename,int N_1,int N_2,double DX1,double DX2)
{
	int type=0;
	print_strain_script(source,filename,DX1,DX2,type);

	string label;
	source.GET_label(label);

	filename = label + filename + "-STRAIN.dat";

	double step_X1,step_X2;
	double x,y,z=0;
	double S[6];

	//        !-------------------------------------------------------------------!
	//        ! OUTPUT FILE OPENING

	ofstream fid_out;
	fid_out.open(filename.c_str());

	ostream *stream;
	stream = &fid_out;


	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->setf(ios::showpos);
	stream->precision(5);


	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	char* dt = ctime(&now);


	*stream << "# Datafile creation: " << dt << endl;

	*stream << endl;

	*stream << "# Parameters of the datafile: " << endl;
	*stream << "# " << "N_1 = " << N_1 << " ,N_2 = " << N_2 << endl;
	*stream << "# " << "DX_1 = " << DX2 << " m, DX_2 = " << DX2 << " m" << endl;

	*stream << endl;


	int input=1;
	int flag_print=10;

	CONSOLE CONSOLE_out(input,stream,flag_print);
	source.PRINT(CONSOLE_out);


	//        !-------------------------------------------------------------------!
	//        ! SET GRID OF OBSERVATION POINTS

	step_X1=DX1/N_1;
	step_X2=DX2/N_2;

	double X_01=DX1*0.5;
	double X_02=DX2*0.5;

	int cont=0;
	int np=N_1*N_2;
	int c[10]={0,0,0,0,0,0,0,0,0,0};

	cout << endl << "Generation of the grid about surface strain" << endl;


	*stream << endl;
	*stream << "#   " << "x (m)"     << "         "    <<  "y (m)"     << "        "
					  << "Exx      " << "         "    <<  "Eyy      " << "        " <<  "Ezz      " << "        "
	  	  	  	  	  << "Exy      " << "         "    <<  "Exz      " << "        " <<  "Eyz      " << endl;
	*stream << endl;

	for(int i=0; i <= N_1; i++)
	{

		x=step_X1*i-X_01;

		for(int j=0; j <= N_2; j++)
		{

			y=step_X2*j-X_02;

			cont = i * N_1 + j + 1;


			source.STRAIN(x,y,z,S);

			//     	!-------------------------------------------------------------------
			//     	! WRITE OUTPUT

			stream->fill(' ');
			*stream << right << setw(8) << x;

			stream->fill(' ');
			*stream << setw(2) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << y ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[2] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[3] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[4] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[5] ;

			*stream << endl;

		}

		fid_out << endl;

		Cycle_counter(cont,np,c);

	}

	fid_out.close();

}




template <class TSOURCE>
void print_strain_script(TSOURCE source,string filename,double DX1,double DX2,int type,int int_shift = 0)
{

	string source_label;
	source.GET_label(source_label);

	string PLOT_Strain = source_label + filename + "-STRAIN.gnuplot";


	string IMAGE_NAME_Sxx = source_label + filename + "-e_xx.eps";

	string IMAGE_NAME_Syy = source_label + filename + "-e_yy.eps";

	string IMAGE_NAME_Szz = source_label + filename + "-e_zz.eps";

	string IMAGE_NAME_Sxy = source_label + filename + "-e_xy.eps";

	string IMAGE_NAME_Sxz = source_label + filename + "-e_xz.eps";

	string IMAGE_NAME_Syz = source_label + filename + "-e_yz.eps";


	if(type == 0)
	{
		filename = source_label + filename + "-STRAIN.dat";
	}
	else if(type != 0)
	{
		filename = source_label + filename + ".dat";
	}


	int i_Sxx = 3 + int_shift;
	int i_Syy = 4 + int_shift;
	int i_Szz = 5 + int_shift;
	int i_Sxy = 6 + int_shift;
	int i_Sxz = 7 + int_shift;
	int i_Syz = 8 + int_shift;

	int i_pDX1 = (int) DX1*0.5*1e-3;
	int i_mDX1 = -i_pDX1;

	int i_pDX2 = (int) DX2*0.5*1e-3;
	int i_mDX2 = -i_pDX2;


	ofstream PLOT_STRAIN;
	PLOT_STRAIN.open(PLOT_Strain.c_str());

	PLOT_STRAIN << "set term postscript eps enhanced															" << endl;

	PLOT_STRAIN << "set size 0.7,1																				" << endl;

	PLOT_STRAIN << "set xlabel 'x (km)'																			" << endl;
	PLOT_STRAIN << "set ylabel 'y (km)'																			" << endl;

	PLOT_STRAIN << "set xtics 1																					" << endl;
	PLOT_STRAIN << "set ytics 1																					" << endl;

	PLOT_STRAIN << "set view map																				" << endl;
	PLOT_STRAIN << "set pm3d																					" << endl;
	PLOT_STRAIN << "set cbrange [-1:1]																			" << endl;
	PLOT_STRAIN << "set cblabel '(x 10^{-4})' rotate by 0 offset -2.5,-12.5										" << endl;
	PLOT_STRAIN << "set palette defined" <<
						  "( 0 '#0d0375',  4 '#1e158c',  8 '#3228a3', 12 '#4a40b7', 16 '#655cc9', 20 '#7C71F6', " <<
			              " 24 '#7787ff', 28 '#77a7ff', 32 '#77c6ff', 36 '#77e6ff', 38 '#9bf5ff', 40 '#b2fdff', " <<
			              " 42 '#dcffff', 44 '#ffffff', 46 '#ffffdc', 48 '#fffdb2', 50 '#FFF59B', 52 '#ffe677', " <<
			              " 56 '#ffc677', 60 '#ffa777', 64 '#ff8777', 68 '#ef614f', 72 '#ed3017', 76 '#dd0404', " <<
			              " 80 '#c40d19', 84 '#b50e19', 88 '#ad202a')											" << endl;

	PLOT_STRAIN << "set nokey																					" << endl;
	PLOT_STRAIN << "set style increment userstyles 																" << endl;
	PLOT_STRAIN << "set contour																					" << endl;
	PLOT_STRAIN << "set cntrparam levels auto																	" << endl;
	PLOT_STRAIN << "set nosurface																				" << endl;
	PLOT_STRAIN 																								  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Sxx <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{xx}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"		      <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Sxx << "*1e4) with lines			" << endl;
	PLOT_STRAIN 																							 	  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Syy <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{yy}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"			  <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Syy << "*1e4) with lines			" << endl;
	PLOT_STRAIN 																								  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Szz <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{zz}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"			  <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Szz << "*1e4) with lines			" << endl;
	PLOT_STRAIN 																								  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Sxy <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{xy}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"		  	  <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Sxy << "*1e4) with lines			" << endl;
	PLOT_STRAIN 																							 	  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Sxz <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{xz}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"			  <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Sxz << "*1e4) with lines			" << endl;
	PLOT_STRAIN 																								  << endl;
	PLOT_STRAIN << "set output '" << IMAGE_NAME_Syz <<"'														" << endl;
	PLOT_STRAIN << "set title '" << source_label << " - e_{yz}'													" << endl;
	PLOT_STRAIN << "splot [" << i_mDX1 << ":" << i_pDX1 << "] [" << i_mDX2 << ":" << i_pDX2 << "] '"			  <<
							filename <<	"' using ($1*1e-3):($2*1e-3):($" << i_Syz << "*1e4) with lines			" << endl;

}




template <class TSOURCE>
void print_map_stress(TSOURCE source,string filename,int N_1,int N_2,double DX1,double DX2)
{
	int type=0;
	print_stress_script(source,filename,DX1,DX2,type);

	string label;
	source.GET_label(label);

	filename = label + filename + "-STRESS.dat";

	double step_X1,step_X2;
	double x,y,z=0;
	double S[6];

	//        !-------------------------------------------------------------------!
	//        ! OUTPUT FILE OPENING

	ofstream fid_out;
	fid_out.open(filename.c_str());

	ostream *stream;
	stream = &fid_out;

	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->setf(ios::showpos);
	stream->precision(5);


	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	char* dt = ctime(&now);


	*stream << "# Datafile creation: " << dt << endl;

	*stream << endl;

	*stream << "# Parameters of the datafile: " << endl;
	*stream << "# " << "N_1 = " << N_1 << " ,N_2 = " << N_2 << endl;
	*stream << "# " << "DX_1 = " << DX2 << " m, DX_2 = " << DX2 << " m" << endl;

	*stream << endl;


	int input=1;
	int flag_print=10;

	CONSOLE CONSOLE_out(input,stream,flag_print);
	source.PRINT(CONSOLE_out);


	//        !-------------------------------------------------------------------!
	//        ! SET GRID OF OBSERVATION POINTS

	step_X1=DX1/N_1;
	step_X2=DX2/N_2;

	double X_01=DX1*0.5;
	double X_02=DX2*0.5;

	int cont=0;
	int np=N_1*N_2;
	int c[10]={0,0,0,0,0,0,0,0,0,0};

	cout << endl << "Generation of the grid about surface stress" << endl;

	for(int i=0; i <= N_1; i++)
	{

		x=step_X1*i-X_01;

		for(int j=0; j <= N_2; j++)
		{

			y=step_X2*j-X_02;

			cont = i * N_1 + j + 1;


			source.STRESS(x,y,z,S);

			//     	!-------------------------------------------------------------------
			//     	! WRITE OUTPUT

			stream->fill(' ');
			*stream << right << setw(6) << x;

			stream->fill(' ');
			*stream << setw(2) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << y ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[2] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[3] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[4] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

			stream->fill(' ');
			*stream << right << setw(6) << S[5] ;

			*stream << endl;

		}

		fid_out << endl;

		Cycle_counter(cont,np,c);

	}

	fid_out.close();

}




template <class TSOURCE>
void print_stress_script(TSOURCE source,string filename,double DX1,double DX2,int type,int int_shift = 0)
{

	string source_label;
	source.GET_label(source_label);

	string PLOT_Stress = source_label + filename + "-STRESS.gnuplot";


	string IMAGE_NAME_Sxx = source_label + filename + "-S_xx.eps";

	string IMAGE_NAME_Syy = source_label + filename + "-S_yy.eps";

	string IMAGE_NAME_Szz = source_label + filename + "-S_zz.eps";

	string IMAGE_NAME_Sxy = source_label + filename + "-S_xy.eps";

	string IMAGE_NAME_Sxz = source_label + filename + "-S_xz.eps";

	string IMAGE_NAME_Syz = source_label + filename + "-S_yz.eps";


	if(type == 0)
	{
		filename = source_label + filename + "-STRESS.dat";
	}
	else if(type != 0)
	{
		filename = source_label + filename + ".dat";
	}


	int i_Sxx = 3 + int_shift;
	int i_Syy = 4 + int_shift;
	int i_Szz = 5 + int_shift;
	int i_Sxy = 6 + int_shift;
	int i_Sxz = 7 + int_shift;
	int i_Syz = 8 + int_shift;


	int i_pDX1 = (int) DX1*0.5*1e-3;
	int i_mDX1 = -i_pDX1;

	int i_pDX2 = (int) DX2*0.5*1e-3;
	int i_mDX2 = -i_pDX2;


	ofstream PLOT_STRESS;
	PLOT_STRESS.open(PLOT_Stress.c_str());

	PLOT_STRESS << "set term postscript eps enhanced															" << endl;

	PLOT_STRESS << "set size 0.7,1																				" << endl;

	PLOT_STRESS << "set xlabel 'x (km)'																			" << endl;
	PLOT_STRESS << "set ylabel 'y (km)'																			" << endl;

	PLOT_STRESS << "set xtics 1																					" << endl;
	PLOT_STRESS << "set ytics 1																					" << endl;

	PLOT_STRESS << "set view map																				" << endl;
	PLOT_STRESS << "set pm3d																					" << endl;
	PLOT_STRESS << "set cbrange [-2:2]																			" << endl;
	PLOT_STRESS << "set cblabel '(MPa)' rotate by 0 offset -2.5,-12.5											" << endl;
	PLOT_STRESS << "set palette defined" <<
						  "( 0 '#0d0375',  4 '#1e158c',  8 '#3228a3', 12 '#4a40b7', 16 '#655cc9', 20 '#7C71F6', " <<
			              " 24 '#7787ff', 28 '#77a7ff', 32 '#77c6ff', 36 '#77e6ff', 38 '#9bf5ff', 40 '#b2fdff', " <<
			              " 42 '#dcffff', 44 '#ffffff', 46 '#ffffdc', 48 '#fffdb2', 50 '#FFF59B', 52 '#ffe677', " <<
			              " 56 '#ffc677', 60 '#ffa777', 64 '#ff8777', 68 '#ef614f', 72 '#ed3017', 76 '#dd0404', " <<
			              " 80 '#c40d19', 84 '#b50e19', 88 '#ad202a')											" << endl;

	PLOT_STRESS << "set nokey																					" << endl;
	PLOT_STRESS << "set style increment userstyles 																" << endl;
	PLOT_STRESS << "set contour																					" << endl;
	PLOT_STRESS << "set cntrparam levels auto																	" << endl;
	PLOT_STRESS << "set nosurface																				" << endl;

	PLOT_STRESS << "set xrange [" << i_mDX1 << ":" << i_pDX1 << "]												" << endl;
	PLOT_STRESS << "set yrange [" << i_mDX2 << ":" << i_pDX2 << "]												" << endl;

	PLOT_STRESS 																							  	  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Sxx <<"'														" << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{xx}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Sxx << "*1e-6) with lines		" << endl;
	PLOT_STRESS 																							 	  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Syy <<"'														" << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{yy}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Syy << "*1e-6) with lines		" << endl;
	PLOT_STRESS 																							  	  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Szz <<"'													"     << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{zz}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Szz << "*1e-6) with lines		" << endl;
	PLOT_STRESS 																								  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Sxy <<"'														" << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{xy}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Sxy << "*1e-6) with lines		" << endl;
	PLOT_STRESS 																							  	  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Sxz <<"'														" << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{xz}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Sxz << "*1e-6) with lines		" << endl;
	PLOT_STRESS 																							  	  << endl;
	PLOT_STRESS << "set output '" << IMAGE_NAME_Syz <<"'														" << endl;
	PLOT_STRESS << "set title '" << source_label << " - {/Symbol \163}_{yz}'									" << endl;
	PLOT_STRESS << "splot '" << filename << "' using ($1*1e-3):($2*1e-3):($" << i_Syz << "*1e-6) with lines		" << endl;

}





template <class TSOURCE>
void print_maps(TSOURCE source,string filename,int N_1,int N_2,double DX1,double DX2)
{

	print_displ_script(source,filename,DX1,DX2,1);
	print_strain_script(source,filename,DX1,DX2,1,3);
	print_stress_script(source,filename,DX1,DX2,1,9);


	string label;
	source.GET_label(label);

	filename = label + filename + ".dat";


	double step_X1,step_X2;
	double x,y,z=0;
	double U[3];
	double Strain[6],Stress[6];

	//        !-------------------------------------------------------------------!
	//        ! OUTPUT FILE OPENING

	ofstream fid_out;
	fid_out.open(filename.c_str());

	ostream *stream;
	stream = &fid_out;


	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->setf(ios::showpos);
	stream->precision(5);


	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	char* dt = ctime(&now);


	*stream << "# Datafile creation: " << dt << endl;

	*stream << endl;

	*stream << "# Parameters of the datafile: " << endl;
	*stream << "# " << "N_1 = " << N_1 << " ,N_2 = " << N_2 << endl;
	*stream << "# " << "DX_1 = " << DX2 << " m, DX_2 = " << DX2 << " m" << endl;

	*stream << endl;


	int input=10;
	int flag_print=1;

	CONSOLE CONSOLE_out(input,stream,flag_print);
	source.PRINT(CONSOLE_out);


	//        !-------------------------------------------------------------------!
	//        ! SET GRID OF OBSERVATION POINTS

	step_X1=DX1/N_1;
	step_X2=DX2/N_2;

	double X_01=DX1*0.5;
	double X_02=DX2*0.5;


	MEDIUM ELASTIC_PARAM;
	source.GET_MEDIUM_PAR(ELASTIC_PARAM);

	double mu = ELASTIC_PARAM.GET_mu();
	double lambda = ELASTIC_PARAM.GET_lambda();

	double Skk_lambda;

	int cont=0;
	int np=N_1*N_2;
	int c[10]={0,0,0,0,0,0,0,0,0,0};

	cout << endl << "Generation of the grid about surface displacement, strain and stress" << endl;

	*stream << endl;
	*stream << "#   " << "x (m)"     << "            " <<  "y (m)"     << "            "
					  << "Ux (m)   " << "         "    <<  "Uy (m)   " << "        " <<  "Uz (m)   " << "   	 "
					  << "Exx      " << "         "    <<  "Eyy      " << "        " <<  "Ezz      " << "        "
	  	  	  	  	  << "Exy      " << "         "    <<  "Exz      " << "        " <<  "Eyz      " << "        "
					  << "Sxx (MPa)" << "         "    <<  "Syy (MPa)" << "        " <<  "Szz (MPa)" << "        "
					  << "Sxy (MPa)" << "         "    <<  "Sxz (MPa)" << "        " <<  "Syz (MPa)" << endl;
	*stream << endl;

	for(int i=0; i <= N_1; i++)
	{

		x=step_X1*i-X_01;

		for(int j=0; j <= N_2; j++)
		{

			y=step_X2*j-X_02;

			cont = i * N_1 + j + 1;


			source.DISPLACEMENT(x,y,z,U);
			source.STRAIN(x,y,z,Strain);

			Skk_lambda = (Strain[0]+Strain[1]+Strain[2]) * lambda;

			Stress[0] = Skk_lambda + 2 * mu * Strain[0];
			Stress[1] = Skk_lambda + 2 * mu * Strain[1];
			Stress[2] = Skk_lambda + 2 * mu * Strain[2];
			Stress[3] = 2 * mu * Strain[3];
			Stress[4] = 2 * mu * Strain[4];
			Stress[5] = 2 * mu * Strain[5];


			//     	!-------------------------------------------------------------------
			//     	! WRITE OUTPUT

     		stream->width(12);

			stream->fill(' ');
			*stream << right << x;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << y ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << U[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << U[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << U[2] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";
			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[2] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[3] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[4] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Strain[5] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[0] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[1] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[2] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[3] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[4] ;

			stream->fill(' ');
			*stream << setw(5) << " ";

     		stream->width(12);

			stream->fill(' ');
			*stream << right << Stress[5] ;

			stream->fill(' ');
			*stream << setw(5) << " ";


			*stream << endl;

		}

		fid_out << endl;

		Cycle_counter(cont,np,c);

	}

	fid_out.close();

}










template <class TSOURCE>
void print_map_radial_displ(TSOURCE source,string filename,int N_r,double r_i,double r_f)
{

	print_radial_displ_script(source,filename,r_i,r_f);

	double z=0;

	double Uz,Ur;

	//        !-------------------------------------------------------------------!
	//        ! OUTPUT FILE OPENING

	string label;
	source.GET_label(label);


	ofstream fid_out;
	fid_out.open(filename.c_str());

	ostream *stream;
	stream = &fid_out;


	stream->fill(' ');
	stream->setf(ios::scientific);
	stream->setf(ios::right);
	stream->setf(ios::showpoint);
	stream->setf(ios::showpos);
	stream->precision(5);


	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	char* dt = ctime(&now);


	*stream << "# Datafile creation: " << dt << endl;

	*stream << endl;

	*stream << "# Parameters of the datafile: " << endl;
	*stream << "# " << "N_r = " << N_r << endl;
	*stream << "# " << "r_i = " << r_i << " m, r_fD = " << r_f << " m" << endl;

	*stream << endl;


	int input=1;
	int flag_print=10;

	CONSOLE CONSOLE_out(input,stream,flag_print);
	source.PRINT(CONSOLE_out);

	double step_r=(r_f-r_i)/N_r;

	double r;

	int dim_r = N_r + 1;

	int c[10]={0,0,0,0,0,0,0,0,0,0};

	cout << endl << "Plot generation of surface displacement (available for symmetric sources)" << endl;


	for(int i=0; i < dim_r; i++)
	{
		r = step_r*i;

		source.DISPLACEMENT(r,z,Uz,Ur);

		//     	!-------------------------------------------------------------------
		//     	! WRITE OUTPUT

 		stream->width(12);

		stream->fill(' ');
		*stream << right << r;

		stream->fill(' ');
		*stream << setw(5) << " ";

 		stream->width(12);

		stream->fill(' ');
		*stream << right << Uz ;

		stream->fill(' ');
		*stream << setw(5) << " ";

 		stream->width(12);

		stream->fill(' ');
		*stream << right << Ur ;

		stream->fill(' ');
		*stream << setw(5) << " ";

		*stream << endl;

		Cycle_counter(i,dim_r,c);

	}

	fid_out.close();

}






template <class TSOURCE>
void print_radial_displ_script(TSOURCE source,string filename,double r_i,double r_f)
{

	string SUFFIX = "_DISPL.gnuplot";

	string source_label;
	source.GET_label(source_label);

	string PLOT_DISPL = source_label + SUFFIX;


	string SUFFIX_Uz = "_Uz.eps";
	string IMAGE_NAME_Uz = source_label + SUFFIX_Uz;

	string SUFFIX_Ur = "_Ur.eps";
	string IMAGE_NAME_Ur = source_label + SUFFIX_Ur;


	int i_r_i = (int) r_i*1e-3;
	int i_r_f = (int) r_f*1e-3;


	ofstream PLOT_SOURCE_DISPL;
	PLOT_SOURCE_DISPL.open(PLOT_DISPL.c_str());

	PLOT_SOURCE_DISPL << "set term postscript eps enhanced														" << endl;

	PLOT_SOURCE_DISPL << "set size 0.7,1																		" << endl;

	PLOT_SOURCE_DISPL << "set nokey																				" << endl;

	PLOT_SOURCE_DISPL << "set xrange [" << i_r_i << ":" << i_r_f << "]											" << endl;
	PLOT_SOURCE_DISPL << "set yrange [0:5]																		" << endl;

	PLOT_SOURCE_DISPL << "set xlabel 'r (km)'																	" << endl;
	PLOT_SOURCE_DISPL << "set ylabel 'U_{z} (m)'																" << endl;

	PLOT_SOURCE_DISPL << "set output '" << IMAGE_NAME_Uz <<"'													" << endl;
	PLOT_SOURCE_DISPL << "set title '" << source_label << " - U_{z}'											" << endl;
	PLOT_SOURCE_DISPL << "plot '" << filename << "' using ($1*1e-3):2 with lines								" << endl;

	PLOT_SOURCE_DISPL << "set xlabel 'r (km)'																	" << endl;
	PLOT_SOURCE_DISPL << "set ylabel 'U_{r} (m)'																" << endl;

	PLOT_SOURCE_DISPL << "set output '" << IMAGE_NAME_Ur <<"'													" << endl;
	PLOT_SOURCE_DISPL << "set title '" << source_label << " - U_{r}'											" << endl;
	PLOT_SOURCE_DISPL << "plot '" << filename << "' using ($1*1e-3):3 with lines								" << endl;

}




#endif  // PRINT_TEMPLATE_H_














