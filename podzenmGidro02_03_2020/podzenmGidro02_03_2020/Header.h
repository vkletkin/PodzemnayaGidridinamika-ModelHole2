#pragma once
#define spaCe1 "                              "
#define spaCe "     "
#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <fstream>

using namespace std;

extern double L;
extern double Left;
extern double Right;

extern double Rbegin;
extern double Rend;


int progonka_3d(
	const vector<double>& c,
	const vector<double>& b,
	const vector<double>& a,
	const vector<double>& d,
	vector <double>& X,
	int n
);

void vvodN(
	int& N
);

double f(
	double x
);

double W(
	double x
);

void vvodNach_bezrazmDan(
	double& P0,
	double& PN,
	double& Rc,
	double& Rk
);

void Get_Setka(
	vector<double>& X,
	int Get
);

double P_An(
	double x
);

double P_An_poluuzel(
	const vector<double>& R, 
	int i
);

double U_An_poluuzel(
	const vector<double>& R, 
	int i
);

double Pop(
	const vector<double>& R,
	int i 
);

double U_chis_poluuzel(
	const vector<double>& P, 
	const vector<double>& R,
	int i
);

double U_chis_poluuzel_pop(
	const vector<double>& P,
	const vector<double>& R,
	int i
);

double Q_chis(
	const vector<double>& P, 
	const vector<double>& R,
	int i
);
double Q_An(
);
void vvod_razmernDan(
	double& L,
	double& k0,
	double& m,
	double& p,
	double& mu
);

void VvodANDVivod_razmernDan(
	double& L,
	double& k0,
	double& m, 
	double& delta_p,
	double& mu, 
	const vector<double>& P,
	const vector<double>& X
);


void resaults(
	const vector<double>& x,
	const vector<double>& P,
	int N
);


double viviod_newyazka(
	const vector<double>& x,
	const vector<double>& P,
	int N
);

void Solve(
	vector<double>& B,
	vector<double>& C,
	vector<double>& A,
	vector<double>& D,
	vector <double>& P,
	vector<double>& x,
	int N,
	double l,
	double volleft,
	double volright
);