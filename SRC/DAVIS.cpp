/*
 * DAVIS.cpp
 *
 *  Created on: 04 nov 2015
 *      Author: dino
 */


#include "DAVIS.h"

#include <iomanip>
#include "Utilities.h"

double pi=acos(-1);



DAVIS::DAVIS() {
	// TODO Auto-generated constructor stub

}

DAVIS::~DAVIS() {
	// TODO Auto-generated destructor stub
}



void DAVIS::SET(MEDIUM medium_par_i,double MXX,double MYY,double MZZ,double MXY,double MXZ,double MYZ)
{
	MEDIUM_PAR = medium_par_i;

	M[0]=MXX;
	M[1]=MYY;
	M[2]=MZZ;
	M[3]=MXY;
	M[4]=MXZ;
	M[5]=MYZ;

	int dimension=3;
	double Matrix[9];
	double eigenvalues[3];

	Matrix[0] = MXX;
	Matrix[1] = MXY;
	Matrix[2] = MXZ;
	Matrix[3] = MXY;
	Matrix[4] = MYY;
	Matrix[5] = MYZ;
	Matrix[6] = MXZ;
	Matrix[7] = MYZ;
	Matrix[8] = MZZ;

	gsl_eigenvalueproblem(dimension,eigenvalues,Matrix);

	M_1 = eigenvalues[0];
	M_2 = eigenvalues[1];
	M_3 = eigenvalues[2];


	double R21=M_2/M_1;
	double R31=M_3/M_1;

	for(int ti=0; ti < 250; ti++)
	{

		double x_left  = TABLE_range[ti+1][1];
		double x_right = TABLE_range[ti][1];

		if((R21 >= x_left) && (R21 < x_right))
		{

			double y_left  = TABLE_range[ti+1][0];
			double y_right = TABLE_range[ti][0];

			double m1=(y_right-y_left)/(x_right-x_left);

			double q1=-m1*x_left+y_left;

			double y_inf,y_sup;

			if(R21 > 0.428572)
			{
				y_inf = m1*R21+q1;
				y_sup = R21;
			}
			else
			{
				y_sup = m1*R21+q1;
				y_inf = R21;
			}

			if((R31 >= y_inf) && (R31 <= y_sup))
			{
				FLAG = 1;
			}

		}
		else
		{
			FLAG = -1;
		}
	}


};




void A2M(float nu, float a, float b, float c, float &M1, float &M2, float &M3)
{

	M1=0;
	M2=0;
	M3=0;


	float R = (1e0 / (8e0 * pi)) * (1e0 - 2e0 * nu) / (1e0 - nu);

	float Q = (3e0 / (8e0 * pi)) * (1e0 / (1e0 - nu));


	if ((a == b) && (a == c))
	{
		M1 = 1;
		M2 = 1;
		M3 = 1;
	}
	else if ((a != b) && (a != c) && (b != c))
	{

		float a2 = pow(a, 2);
		float b2 = pow(b, 2);
		float c2 = pow(c, 2);

		float coef_abc = (4e0 / 3e0) * pi * a * b * c;

		float Ia = coef_abc * rd(b2, c2, a2);

		float Ic = coef_abc * rd(a2, b2, c2);

		float Ib = 4 * pi - Ia - Ic;


		float Iab = (Ib - Ia) / (3 * (a2 - b2));
		float Iba = Iab;

		// Provo a interpretare
		float Iac = (Ic - Ia) / (3 * (a2 - c2));
		float Ica = Iac;
		//
		float Iaa = ((4 * pi) / (3 * a2)) - Iab - Iac;

		// Provo a interpretare
		float Ibc = (Ic - Ib) / (3 * (b2 - c2));
		float Icb = Ibc;

		float Ibb = ((4 * pi) / (3 * b2)) - Iba - Ibc;

		// Provo a interpretare
		float Icc = ((4 * pi) / (3 * c2)) - Ica - Icb;


		float S11 = Q * a2 * Iaa + R * Ia;
		float S12 = Q * b2 * Iab - R * Ia;
		float S13 = Q * c2 * Iac - R * Ia;

		float S21 = Q * a2 * Iba - R * Ib;
		float S22 = Q * b2 * Ibb + R * Ib;
		float S23 = Q * c2 * Ibc - R * Ib;

		float S31 = Q * a2 * Ica - R * Ic;
		float S32 = Q * b2 * Icb - R * Ic;
		float S33 = Q * c2 * Icc + R * Ic;


		float A[9] = {0,0,0,0,0,0,0,0,0};

		A[0] = (S11 - 1) - nu * S12 - nu * S13;
		A[3] = -nu * (S11 - 1) + S12 - nu * S13;
		A[6] = -nu * (S11 - 1) - nu * S12 + S13;

		A[1] = S21 - nu * (S22 - 1) - nu * S23;
		A[4] = -nu * S21 + (S22 - 1) - nu * S23;
		A[7] = -nu * S21 - nu * (S22 - 1) + S23;

		A[2] = S31 - nu * S32 - nu * (S33 - 1);
		A[5] = -nu * S31 + S32 - nu * (S33 - 1);
		A[8] = -nu * S31 - nu * S32 + (S33 - 1);


		float B[3] = {0,0,0};

		B[0] = 1 - 2 * nu;
		B[1] = B[0];
		B[2] = B[0];

		int TOT = 3;

		int ncB = 1;
		int IPIV[TOT];
		int INFO;

//		cout << endl << "Try to find a solution (using SGESV routine)" << endl;

		sgesv_(&TOT, &ncB, A, &TOT, IPIV, B, &TOT, &INFO);

		if(INFO == 0)
		{

//			cout << endl << "Solution identified (INFO = " << INFO << ")" << endl;

			M1 = B[0];
			M2 = B[1];
			M3 = B[2];

		}
		else
		{

//			cout << endl << "INFO =" << INFO << endl;

			if (INFO < 0)
			{
				cout << "If INFO = -i, the i-th argument had an illegal value." << endl;
			}

			if (INFO > 0)
			{
				cout << "The factor U(i,i) is exactly singular (A = P * L * U) and so the solution could not be computed." << endl;
				//			exit(0);
			}

		}

	}
	else if((a == b) || (a == c) || (b == c))
	{

		float a1=a;
		float b1=b;
		float c1=c;

		float scambio;

		if(a == c)
		{
			scambio = b;
			b = c;
			c = scambio;
		}
		else if(b == c)
		{
			scambio = a;
			a = b;
			c = scambio;
		}

		float a2 = pow(a, 2);
		float b2 = pow(b, 2);
		float c2 = pow(c, 2);


		float e = c / a;

		float Ia;

		Ia = 0;

		if(e < 1)
		{
			Ia = (2 * pi * e / sqrt(pow(1 - pow(e, 2), 3)))	* (acos(e) - e * sqrt(1 - pow(e, 2)));
		}
		else if(e > 1)
		{
			Ia = (2 * pi * e / sqrt(pow(pow(e, 2) - 1, 3)))	* (e * sqrt(pow(e, 2) - 1) - acosh(e));
		}

		float Ib = Ia;

		float Ic = 4e0 * pi - 2e0 * Ia;


		float Iac, Ica, Icb;

		Iac = (1e0 / (pow(a, 2) * (pow(e, 2) - 1))) * (Ia - (4e0 / 3e0) * pi);
		Ica = Iac;
		Icb = Iac;

		float Ibc = Icb;

		float Iab;

		Iab = (1e0 / (4e0 * pow(a, 2) * (pow(e, 2) - 1))) * (-Ia + (4e0 / 3e0) * pi * pow(e, 2));

		float Iba = Iab;


		float Iaa;

		Iaa = 3e0 * Iab;

		float Ibb = Iaa;


		float Icc;

		Icc = (2e0 / (pow(a, 2) * (pow(e, 2) - 1)))	* (-Ia + 2e0 * pi * (1 - (1e0 / (3e0 * pow(e, 2)))));


		float S21 = Q * a2 * Iba - R * Ib;
		float S22 = Q * b2 * Ibb + R * Ib;
		float S23 = Q * c2 * Ibc - R * Ib;

		float S31 = Q * a2 * Ica - R * Ic;
		float S32 = Q * b2 * Icb - R * Ic;
		float S33 = Q * c2 * Icc + R * Ic;


		float A[4]={0,0,0,0};


		A[0] = (1e0 - nu) * (S21 + S22 - 1e0) - 2e0 * nu * S23;
		A[2] = -nu * (S21 + S22 - 1e0) + S23;

		A[1] = (1e0 - nu) * (S31 + S32) - 2e0 * nu * (S33 - 1e0);
		A[3] = -nu * (S31 + S32) + (S33 - 1e0);


		float B[2]={0,0};

		B[0] = 1e0 - 2e0 * nu;
		B[1] = B[0];


		int TOT = 2;

		int ncB = 1;
		int IPIV[TOT];
		int INFO;

//		cout << endl << "Try to find a solution (using SGESV routine)" << endl;

		sgesv_(&TOT, &ncB, A, &TOT, IPIV, B, &TOT, &INFO);

		if (INFO == 0)
		{
//			cout << endl << "Solution identified (INFO = " << INFO << ")" << endl;

			if(a1 == b1)
			{
				M1 = B[0];
				M2 = B[0];
				M3 = B[1];
			}
			if(a1 == c1)
			{
				M1 = B[0];
				M2 = B[1];
				M3 = B[0];
			}
			else if(b1 == c1)
			{
				M1 = B[1];
				M2 = B[0];
				M3 = B[0];
			}
		}
		else
		{

			if (INFO < 0)
			{
				cout << "If INFO = -i, the i-th argument had an illegal value."	<< endl;
			}

			if (INFO > 0)
			{
				cout << "The factor U(i,i) is exactly singular (A = P * L * U) and so the solution could not be computed." << endl;
//				exit(0);
			}

		}

	}

}















