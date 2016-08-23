/*
 * FIALKO.cpp
 *
 *  Created on: 20/ott/2014
 *      Author: dino
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>


using namespace std;


#include "FIALKO.h"
#include "MEDIUM.h"
#include "PRINT_TEMPLATE.h"


FIALKO::FIALKO() {
	// TODO Auto-generated constructor stub

}

FIALKO::~FIALKO() {
	// TODO Auto-generated destructor stub
}




vector<double> fpkernel(double h,double ti,vector<double> r,int n);

//void fredholm(vector<double> &fi,vector<double> &psi,vector<double> &t,vector<double> &Wt,double h,int m,double er)
void FIALKO::SET(MEDIUM MEDIUM_PAR_i,double DP_i)
{

	int m = 2;
	double eps = 1e-5;

	MEDIUM_PAR = MEDIUM_PAR_i;

	DP = DP_i;

// function [fi,psi,t,Wt]=fredholm(h,m,er)
// fi,psi: basis functions
// t: interval of integration
// m: size(t)

	double pi=acos(-1);
	double lamda=2/pi;

	vector<double> Root(16);

	Root =
	{
		-0.989400934991649932596,
		-0.944575023073232576078,
		-0.865631202387831743880,
		-0.755404408355003033895,
		-0.617876244402643748447,
		-0.458016777657227386342,
		-0.281603550779258913230,
		-0.095012509837637440185,
		 0.095012509837637440185,
		 0.281603550779258913230,
		 0.458016777657227386342,
		 0.617876244402643748447,
		 0.755404408355003033895,
		 0.865631202387831743880,
		 0.944575023073232576078,
		 0.989400934991649932596
	};


	vector<double> Weight(16);

	Weight =
	{
		0.027152459411754094852,
		0.062253523938647892863,
		0.095158511682492784810,
		0.124628971255533872052,
		0.149595988816576732081,
		0.169156519395002538189,
		0.182603415044923588867,
		0.189450610455068496285,
		0.189450610455068496285,
		0.182603415044923588867,
		0.169156519395002538189,
		0.149595988816576732081,
		0.124628971255533872052,
		0.095158511682492784810,
		0.062253523938647892863,
		0.027152459411754094852
	};


	int NumLegendreTerms=Root.size();

	int dim=m*NumLegendreTerms;

	t.resize(dim,0);

	Wt=t;

	for(int k=0; k < m; k++)
	{
		for(int i=0; i < NumLegendreTerms; i++)
		{
			double d1=1e0/m;
			double t1=d1*k;
			double r1=d1*(k+1);

			int j = NumLegendreTerms*k+i;

			t[j] = Root[i]*(r1-t1)*0.5+(r1+t1)*0.5;
			Wt[j] = 0.5/m*Weight[i];
		}
	}

	vector<double> fi1(m*NumLegendreTerms);

	for(int i=0; i < dim; i++)
	{
		fi1[i] = -lamda*t[i];
	}

	vector<double> psi1(m*NumLegendreTerms,0);

	fi.resize(dim);
	psi.resize(dim);

	fi.assign(dim,0);
	psi.assign(dim,0);

	double res=1e9;

	vector<double> K1(dim),K2(dim),K3(dim),K4(dim);
	vector<double> prod_fi1(dim),prod_fi4(dim);

	while(res > eps)
	{
		double sum_fi=0,sum_psi=0;

		for(int i=0; i < m*NumLegendreTerms; i++)
		{
			sum_fi=0;
			sum_psi=0;

			K1 = fpkernel(h,t[i],t,1);
			K3 = fpkernel(h,t[i],t,3);

			K2 = fpkernel(h,t[i],t,2);
			K4 = fpkernel(h,t[i],t,4);

			for(int j=0; j < m*NumLegendreTerms; j++)
			{
				sum_fi += Wt[j]*(fi1[j]*K1[j] + psi1[j]*K3[j]);
				sum_psi += Wt[j]*(psi1[j]*K2[j] + fi1[j]*K4[j]);
			}

			fi[i]=(-t[i]+sum_fi)*lamda;
			psi[i]=sum_psi*lamda;
		}

		vector<double> diff(dim,0);

		for(int i=0; i < m*NumLegendreTerms; i++)
		{
			diff[i]=abs(fi1[i]-fi[i]);
		}

		double Max_fi = *max_element(diff.begin(),diff.end());
		int Pos_Max_fi = distance(diff.begin(), max_element(diff.begin(), diff.end()));

		double fim = Max_fi / abs(fi[Pos_Max_fi]);

		for(int i=0; i < m*NumLegendreTerms; i++)
		{
			diff[i]=abs(psi1[i]-psi[i]);
		}

		double Max_psi = *max_element(diff.begin(),diff.end());
		int Pos_Max_psi = distance(diff.begin(), max_element(diff.begin(), diff.end()));

		double psim = Max_psi / abs(psi[Pos_Max_psi]);

		res=max(fim,psim);

		fi1=fi;
		psi1=psi;
	}

}




vector<double> KG(vector<double> s,double p);
vector<double> KERN(vector<double> w,double p);

vector<double> fpkernel(double h,double t,vector<double> r,int n)
{

	double p=4*pow(h,2);

	double Dlt=1e-6;

	int dim = r.size();


	vector<double> K(dim);

	vector<double> R1(dim),R2(dim),R3(dim),R4(dim);

	vector<double> a(dim),b(dim);
	vector<double> y(dim),z(dim);

	vector<double> g(dim),s(dim);
	vector<double> trbl(dim),rs(dim);

	vector<double> c(dim),d(dim);


	switch (n)
	{
	case 1:	//KN

		for(int i=0; i < dim; i++)
		{
			R1[i] = t-r[i];
			R2[i] = t+r[i];
		}

		for(int i=0; i < dim; i++)
		{
			R3 = KG(R1,p);
			R4 = KG(R2,p);
		}

		for(int i=0; i < dim; i++)
		{
			K[i]=p*h*(R3[i] - R4[i]);
		}

		break;

	case 2:	//KN1

		Dlt=1e-6;

		for(int i=0; i < dim; i++)
		{
			a[i] = t+r[i];
			b[i] = t-r[i];
		}

		for(int i=0; i < dim; i++)
		{
			y[i] = pow(a[i],2);
			z[i] = pow(b[i],2);
		}

		for(int i=0; i < dim; i++)
		{
			g[i] = 2*p*h*(pow(p,2)+6*p*(pow(t,2)+pow(r[i],2))+5*pow(a[i]*b[i],2));
		}

		for(int i=0; i < dim; i++)
		{
			s[i] = pow((p+z[i])*(p+y[i]),2);
		}

		for(int i=0; i < dim; i++)
		{
			s[i] = g[i]/s[i];
		}


		for(int i=0; i < dim; i++)
		{
			if(t < Dlt)
			{
				trbl[i] = -4*h/(p+pow(r[i],2));
			}
			else
			{
				trbl[i] = h/t/r[i]*log((p+z[i])/(p+y[i]));
			}
		}

		for(int i=0; i < dim; i++)
		{
			R1 = KERN(b,p);
			R2 = KERN(a,p);
		}

		for(int i=0; i < dim; i++)
		{
			K[i] = trbl[i]+s[i]+h*(R1[i]+R2[i]);
		}

		break;

	case 3:	//KM

		for(int i=0; i < dim; i++)
		{
			y[i] = pow(t+r[i],2);
			z[i] = pow(t-r[i],2);
		}

		for(int i=0; i < dim; i++)
		{
			a[i] = pow((p+y[i])*(p+z[i]),2);
			c[i] = t+r[i];
			d[i] = t-r[i];
		}

		for(int i=0; i < dim; i++)
		{
			b[i] = p*t*((3*pow(p,2)-pow(c[i]*d[i],2)+2*p*(pow(t,2)+pow(r[i],2)))/a[i]);
		}


		R1=KG(c,p);
		R2=KG(d,p);

		for(int i=0; i < dim; i++)
		{
			a[i] = p/2*(c[i]*R1[i]+d[i]*R2[i]);
		}

		for(int i=0; i < dim; i++)
		{
			K[i] = b[i]-a[i];
		}

		break;

	case 4:	////KM1(t,r)=KM(r,t);

		for(int i=0; i < dim; i++)
		{
			y[i] = pow(t+r[i],2);
			z[i] = pow(t-r[i],2);
		}

		for(int i=0; i < dim; i++)
		{
			a[i] = pow((p+y[i])*(p+z[i]),2);
			c[i] =  t+r[i];
			d[i] = -t+r[i];
		}


		for(int i=0; i < dim; i++)
		{
			b[i] = p*r[i]*((3*pow(p,2)-pow(c[i]*d[i],2)+2*p*(pow(t,2)+pow(r[i],2)))/a[i]);
		}

		R1=KG(c,p);
		R2=KG(d,p);

		for(int i=0; i < dim; i++)
		{
			a[i] = p/2*(c[i]*R1[i]+d[i]*R2[i]);
		}

		for(int i=0; i < dim; i++)
		{
			K[i] = b[i]-a[i];
		}

		break;

	} //switch n

	return K;


}




vector<double> KG(vector<double> s,double p)
{

	int dim = s.size();

	vector<double> KGo(dim);
	vector<double> y(dim),z(dim);

	for(int i=0; i < dim; i++)
	{
		z[i] = pow(s[i],2);
	}

	for(int i=0; i < dim; i++)
	{
		y[i] = p+z[i];
	}

	for(int i=0; i < dim; i++)
	{
		KGo[i] = (3*p-z[i])/pow(y[i],3);
	}

	return KGo;
}




vector<double> KERN(vector<double> w,double p)
{
	int dim = w.size();

	vector<double> KERNo(dim);
	vector<double> u(dim);

	for(int i=0; i < dim; i++)
	{
		u[i] = pow((p+pow(w[i],2)),3);
	}

	for(int i=0; i < dim; i++)
	{
		KERNo[i] = 2*(pow(p,2)/2+pow(w[i],4)-p*pow(w[i],2)/2)/u[i];
	}

	return KERNo;
}




vector<double> Q(double h,vector<double> t,double r,int n);



void FIALKO::DISPLACEMENT(double r,double z,double &Ur,double &Uz)
{
	if(z == 0)
	{

		// function [Uz,Ur]=intgr(r,fi,psi,h,Wt,t)
		// Uz(r),Ur(r) - vertical and radial displacements
		// fi,psi: basis functions
		// t: interval of integration

		int dim_t=t.size();

		vector<double> Q1(dim_t),Q2(dim_t),Q3(dim_t);
		vector<double> Q4(dim_t),Q5(dim_t),Q6(dim_t),Q7(dim_t),Q8(dim_t);

		double sum_uz=0,sum_ur=0;

		double rp=r/R;

		Q1=Q(h,t,rp,1);
		Q2=Q(h,t,rp,2);
		Q3=Q(h,t,rp,3);

		Q4=Q(h,t,rp,4);
		Q5=Q(h,t,rp,5);
		Q6=Q(h,t,rp,6);
		Q7=Q(h,t,rp,7);
		Q8=Q(h,t,rp,8);


		sum_uz=0;
		sum_ur=0;

		for(int j=0; j < dim_t; j++)
		{
			sum_ur += Wt[j]*(psi[j]*((Q4[j]-h*Q5[j])/t[j] - Q6[j]+h*Q7[j])- h*fi[j]*Q8[j]);
			sum_uz += Wt[j]*(fi[j]*(Q1[j]+h*Q2[j]) + psi[j]*(Q1[j]/t[j]-Q3[j]));
		}


		double coef=2e0*(1-MEDIUM_PAR.GET_nu())*R*DP/MEDIUM_PAR.GET_mu();

		Ur =  sum_ur*coef;
		Uz = -sum_uz*coef;

	}
	else
	{
		cout << "FIALKO: displacement computation allowed only for z=0." << endl;
	}
}




void FIALKO::DISPLACEMENT(double x,double y,double z,double U[])
{
	double xp=x-x0;
	double yp=y-y0;

	double r = sqrt(pow(xp,2)+pow(yp,2));

	double Uz,Ur;

	DISPLACEMENT(r,z,Uz,Ur);

	double teta=atan2(yp,xp);

	U[0]=Ur*cos(teta);
	U[1]=Ur*sin(teta);
	U[2]=Uz;
}




vector<double> Q(double h,vector<double> t,double r,int n)
{

	int dim = t.size();

	vector<double> K(dim);

	// Kernels calculation

	vector<double> E(dim);

	double term=pow(h,2)+pow(r,2);

	for(int i=0; i < dim; i++)
	{
		E[i] = term-pow(t[i],2);
	}

	vector<double> D(dim);

	for(int i=0; i < dim; i++)
	{
		D[i] = sqrt(pow(E[i],2)+4*pow(h,2)*pow(t[i],2));
	}

	vector<double> D3(dim);

	for(int i=0; i < dim; i++)
	{
		D3[i] = pow(D[i],3);
	}


	switch (n)
	{
	case 1:	//Q1

		for(int i=0; i < dim; i++)
		{
			K[i] = sqrt(2)*h*t[i]/(D[i]*sqrt(D[i]+E[i]));
		}

		break;

	case 2:	//Q2

		for(int i=0; i < dim; i++)
		{
			K[i] = 1/sqrt(2)/D3[i]*(h*sqrt(D[i]-E[i])*(2*E[i]+D[i])-t[i]*sqrt(D[i]+E[i])*(2*E[i]-D[i]));
		}

		break;

	case 3:	//Q3

		for(int i=0; i < dim; i++)
		{
			K[i] = 1/sqrt(2)/D3[i]*(h*sqrt(D[i]+E[i])*(2*E[i]-D[i])+t[i]*sqrt(D[i]-E[i])*(2*E[i]+D[i]));
		}

		break;

	case 4:	//Q4

		for(int i=0; i < dim; i++)
		{
			K[i] = t[i]/r-sqrt(D[i]-E[i])/r/sqrt(2);
		}

		break;

	case 5:	//Q5

		for(int i=0; i < dim; i++)
		{
			K[i] = -(h*sqrt(D[i]-E[i])-t[i]*sqrt(D[i]+E[i]))/D[i]/r/sqrt(2);
		}

		break;

	case 6:	//Q6

		for(int i=0; i < dim; i++)
		{
			K[i] = 1/r-(h*sqrt(D[i]+E[i])+t[i]*sqrt(D[i]-E[i]))/D[i]/r/sqrt(2);
		}
		break;

	case 7:	//Q7

		for(int i=0; i < dim; i++)
		{
			K[i] = r*sqrt(D[i]+E[i])*(2*E[i]-D[i])/D3[i]/sqrt(2);
		}

		break;

	case 8:	//Q8

		for(int i=0; i < dim; i++)
		{
			K[i] = r*sqrt(D[i]-E[i])*(2*E[i]+D[i])/D3[i]/sqrt(2);
		}

		break;

	case 41:	//Q4*r

		for(int i=0; i < dim; i++)
		{
			K[i] = t[i]-sqrt(D[i]-E[i])/sqrt(2);
		}

		break;

	case 51:	//Q5*r

		for(int i=0; i < dim; i++)
		{
			K[i] = -(h*sqrt(D[i]-E[i])-t[i]*sqrt(D[i]+E[i]))/D[i]/sqrt(2);
		}

		break;

	case 61:	//Q6*r

		for(int i=0; i < dim; i++)
		{
			K[i] = 1-(h*sqrt(D[i]+E[i])+t[i]*sqrt(D[i]-E[i]))/D[i]/sqrt(2);
		}

		break;

	}

	return K;

}






void FIALKO::PRINT(CONSOLE &out)
{

	int int_flag,flag_print;

	string filename;

	out.GET_int_flag(int_flag);
	out.GET_flag_print(flag_print);

	ostream *stream;
	stream = out.out_stream;

	switch(int_flag)
	{
		case (0):
				PRINT_GEOM_PROP(stream);
				PRINT_STRESS_PROP(stream);
		break;

		case (1):
				PRINT_GEOM_PROP(stream);
		break;

		case (2):
				PRINT_STRESS_PROP(stream);
		break;

		case (10):
			    PRINT_MAIN_FEATURES(stream,flag_print);
				break;

		default:
				*stream << "Not implemented for this source type" << endl;
	}

}




void FIALKO::PRINT_GEOM_PROP(ostream *out_stream)
{

	*out_stream << label << " - Geometrical properties" << endl;

	*out_stream << endl;

	*out_stream << "Position of source center at:" << endl;
	*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
	*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
	*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

	*out_stream << endl;

	*out_stream << "Source radius = " << R << " m" << endl;

	*out_stream << endl;

}




void FIALKO::PRINT_STRESS_PROP(ostream *out_stream)
{
	*out_stream << label << " - Stress properties" << endl;

	*out_stream << endl;

	*out_stream << "Source overpressure:" << endl;
	*out_stream << "DP = " << DP*1e-6 << " MPa" << endl;

	*out_stream << endl;
}




void FIALKO::PRINT(DATAFILE &out)
{

	int flag_datatype;
	out.GET_flag_datatype(flag_datatype);

	if(flag_datatype == 0)
	{
		int int_flag;
		out.GET_int_flag(int_flag);

		string filename;

		int Nr;
		double ri,rf;

		out.GET_file_parameters(filename,Nr,ri,rf);

		filename = label + filename + ".dat";

		switch(int_flag)
		{
			case (1):
				print_map_radial_displ(*this,filename,Nr,ri,rf);

			break;

			default:
				cout << "Opzione " << int_flag << " non disponibile";
				break;
		}

	}
	else if(flag_datatype == 1)
	{
		int int_flag;
		out.GET_int_flag(int_flag);

		string filename;

		int N_1,N_2;
		double DX1,DX2;

		out.GET_file_parameters(filename,N_1,N_2,DX1,DX2);

		switch(int_flag)
		{
			case (1):
            	print_map_displ(*this,filename,N_1,N_2,DX1,DX2);
				break;

			default:
				cout << "Opzione " << int_flag << " non disponibile";
				break;
		}
	}

}




void FIALKO::PRINT_MAIN_FEATURES(ostream *out_stream,int flag_print)
{

    out_stream->fill(' ');
    out_stream->setf(ios::right);
	out_stream->setf(ios::showpoint);
	out_stream->precision(5);

	double mu=MEDIUM_PAR.GET_mu();
	double nu=MEDIUM_PAR.GET_nu();

	if(flag_print == 0)
	{
		*out_stream << label << " - Main features" << endl;

		*out_stream << endl;

		*out_stream << "Position of source center at:" << endl;
		*out_stream << "x0 = " << x0*1e-3 << " km"     << endl;
		*out_stream << "y0 = " << y0*1e-3 << " km"     << endl;
		*out_stream << "z0 = " << z0*1e-3 << " km"     << endl;

		*out_stream << endl;

		*out_stream << "Geometric properties:" << endl;
		*out_stream << "R = " << R << " m" << endl;

		*out_stream << endl;

		*out_stream << "Overpressure:" << endl;
		*out_stream << "DP = " << DP << " Pa" << endl;

		*out_stream << endl;

		*out_stream << "The elastic parameter of the medium are:" << endl;
		*out_stream << "mu = " << mu << endl;
		*out_stream << "nu = " << nu << endl;

		*out_stream << endl;
	}
	else
	{
		*out_stream << "# " << label << " - Main features" << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "Position of source center at:" << endl;
		*out_stream << "# " << "x0 = " << x0*1e-3 << " km"     << endl;
		*out_stream << "# " << "y0 = " << y0*1e-3 << " km"     << endl;
		*out_stream << "# " << "z0 = " << z0*1e-3 << " km"     << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "Overpressure:" << endl;
		*out_stream << "# " << "DP = " << DP << " Pa" << endl;

		*out_stream << "# " << endl;

		*out_stream << "# " << "The elastic parameter of the medium are:" << endl;
		*out_stream << "# " << "mu = " << mu << endl;
		*out_stream << "# " << "nu = " << nu << endl;

		*out_stream << "# " << endl;
	}

}
