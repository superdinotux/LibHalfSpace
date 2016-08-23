#include <stdio.h>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "cmath"


#include "Utilities.h"


using namespace std;



void Cycle_counter(int jf,int nsg,int c[])
{

        int perc_cont;
        double rjf,cont;

        rjf=jf;
        cont = rjf/nsg;
        perc_cont = cont*1e2;

        switch (perc_cont)
		{
        case (10):
            if(c[0] == 0)
            {
            	cout << endl << "Execution: 10% completed" << endl;
                c[0]=1;
            }
        	break;

        case (20):
            if(c[1] == 0)
            {
                cout << "Execution: 20% completed" << endl;
                c[1]=1;
            }
    	break;

        case (30):
            if(c[2] == 0)
            {
                cout << "Execution: 30% completed" << endl;
                c[2]=1;
            }
    	break;

        case (40):
            if(c[3] == 0)
            {
                cout << "Execution: 40% completed" << endl;
                c[3]=1;
            }
    	break;

        case (50):
            if(c[4] == 0)
            {
                cout << "Execution: 50% completed" << endl;
                c[4]=1;
            }
    	break;

        case (60):
            if(c[5] == 0)
            {
                cout << "Execution: 60% completed" << endl;
                c[5]=1;
            }
    	break;

        case (70):
            if(c[6] == 0)
            {
                cout << "Execution: 70% completed" << endl;
                c[6]=1;
            }
        break;

        case (80):
            if(c[7] == 0)
            {
                cout << "Execution: 80% completed" << endl;
                c[7]=1;
            }
    	break;

        case (90):
            if(c[8] == 0)
            {
                cout << "Execution: 90% completed" << endl;
                c[8]=1;
            }
    	break;

        case (100):
            if(c[9] == 0)
            {
                cout << "Execution completed" << endl;
                c[9]=1;
            }

		}


}





void gsl_eigenvalueproblem(int dimension,double eigenvalues[],double Matrix[])
{

	  gsl_matrix_view m = gsl_matrix_view_array (Matrix, dimension, dimension);

	  gsl_vector *eval = gsl_vector_alloc (dimension);
	  gsl_matrix *evec = gsl_matrix_alloc (dimension, dimension);

	  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dimension);

	  gsl_eigen_symmv (&m.matrix, eval, evec, w);

	  gsl_eigen_symmv_free (w);

	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);


      for(int i = 0; i < dimension; i++)
      {
           double eval_i = gsl_vector_get (eval, i);
           gsl_vector_view evec_i = gsl_matrix_column (evec, i);

           eigenvalues[i]=eval_i;

//           printf ("\neigenvalue = %g\n", eval_i);
//           printf ("eigenvector = \n");
//           gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
      }


     gsl_vector_free (eval);
     gsl_matrix_free (evec);

     return;


};





void ref_func(double &x,double &y,double &z,double X1,double X2,double X3_plane,int flag_plane)
{

	switch (flag_plane)
	{
		case (1):
			// Horizontal plane (parallel to xy plane)

			x = X1;
			y = X2;
			z = X3_plane;

			break;

		case (2):
			// Vertical plane (parallel to xz plane)

			x = X1;
			y = X3_plane;
			z = X2;

			break;

		case (3):
			// Vertical plane (parallel to yz plane)

			x = X3_plane;
			y = X1;
			z = X2;

			break;

		break;

	}

}





float rd(float x, float y, float z)
//Computes Carlson’s elliptic integral of the second kind, rd (x, y, z). x and y must be non-
//negative, and at most one can be zero. z must be positive. TINY must be at least twice the
//negative 2/3 power of the machine overflow limit. BIG must be at most 0.1 × ERRTOL times
//the negative 2/3 power of the machine underflow limit.
{
	float alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty,
	sqrtz, sum, xt, yt, zt;

	//	if (FMIN(x,y) < 0.0 || FMIN(x+y,z) < TINY || FMAX(FMAX(x,y),z) > BIG)
	//		cout << "invalid arguments in rd") << endl;

	xt = x;
	yt = y;
	zt = z;

	sum = 0.0;
	fac = 1.0;

	do{
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);

		alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
		sum += fac / (sqrtz * (zt + alamb));
		fac = 0.25 * fac;

		xt = 0.25 * (xt + alamb);
		yt = 0.25 * (yt + alamb);
		zt = 0.25 * (zt + alamb);

		ave = 0.2 * (xt + yt + 3.0 * zt);

		delx = (ave - xt) / ave;
		dely = (ave - yt) / ave;
		delz = (ave - zt) / ave;

		//} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	} while (FMAX(FMAX(abs(delx),abs(dely)),abs(delz)) > ERRTOL);

	ea = delx * dely;
	eb = delz * delz;
	ec = ea - eb;
	ed = ea - 6.0 * eb;
	ee = ed + ec + ec;

	return 3.0 * sum
			+ fac
			* (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee)
					+ delz
					* (C2 * ee
							+ delz * (-C3 * ec + delz * C4 * ea)))
							/ (ave * sqrt(ave));

}






extern "C" {
void dc3d_(float *alpha,float *X,float *Y,float *Z,float *DEPTH,float *DIP,
			float *AL1,float *AL2,float *AW1,float *AW2,float *DISL1,float *DISL2,float *DISL3,
			float *UX,float *UY,float *UZ,
			float *UXX,float *UYX,float *UZX,float *UXY,float *UYY,float *UZY,float *UXZ,float *UYZ,float *UZZ,int *IRET);
}



void DC3D(float alpha,float X,float Y,float Z,float DEPTH,float DIP,
		float AL1,float AL2,float AW1,float AW2,float DISL1,float DISL2,float DISL3,
		float &UX,float &UY,float &UZ,
		float &UXX,float &UYX,float &UZX,float &UXY,float &UYY,float &UZY,float &UXZ,float &UYZ,float &UZZ,int &IRET)
{
	  dc3d_(&alpha,&X,&Y,&Z,&DEPTH,&DIP,&AL1,&AL2,&AW1,&AW2,&DISL1,&DISL2,&DISL3,&UX,&UY,&UZ,&UXX,&UYX,&UZX,&UXY,&UYY,&UZY,&UXZ,&UYZ,&UZZ,&IRET);
}






