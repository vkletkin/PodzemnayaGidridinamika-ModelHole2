#pragma once

#include "Header.h"
double L;
double Left;
double Right;

double Rbegin;
double Rend;

double Pbegin;
double Pend;

void vvodN(int& N)
{
	cout << "N(точки) = ";
	cin >> N;
}

double f(double x)
{
	if (x < L)return Left;
	else return Right;
}

double W(double x)
{
	return x;
}

double P_An(double r)
{
	return Pbegin + (Pend - Pbegin) * (log( r / Rbegin ) / log(Rend / Rbegin));
}

double P_An_poluuzel(const vector<double>& R, int i = 0)
{
	double R_poluzel = (R[i] + R[i + 1]) * 0.5;
	return P_An(R[i]) + (P_An(R[i + 1]) - P_An(R[i])) * (log(R_poluzel / R[i]) / log(R[i + 1] / R[i]));
}

double U_An_poluuzel(const vector<double>& R, int i = 0)
{
	double R_poluzel = (R[i] + R[i + 1]) * 0.5;
	return (-f(R_poluzel) / R_poluzel) * (P_An(R[i + 1]) - P_An(R[i])) / log(R[i + 1] / R[i]);
}

double Pop(const vector<double>& R, int i = 0 )
{
	if (i == 0)
	{
		double R_poluzel = (R[i] + R[i + 1]) * 0.5;
		return (R[i + 1] - R[i]) / (log(R[i + 1] / R[i]) * R_poluzel);
	}
	else
		return 1;
}

double U_chis_poluuzel(const vector<double>& P, const vector<double>& R, int i)
{
	double R_poluzel = (R[i] + R[i + 1]) * 0.5;
	return -f(R_poluzel) * (P[i + 1] - P[i]) / (R[i + 1] - R[i]);
}

double U_chis_poluuzel_pop(const vector<double>& P, const vector<double>& R, int i)
{
	double R_poluzel = (R[i] + R[i + 1]) * 0.5;
	return -f(R_poluzel) * (P[i + 1] - P[i]) / (R[i + 1] - R[i]) * Pop(R);
}

double Q_chis(const vector<double>& P, const vector<double>& R, int i = 0)
{
	double R_poluzel = (R[i] + R[i + 1]) / 2;

	double U_poluzel = U_chis_poluuzel(P, R, i);

	return -2*3.14159*(R_poluzel* U_poluzel);
}

double Q_An()
{
	return 2 * 3.14159 * (Pend - Pbegin) / log(Rend / Rbegin);
}


void vvodNach_bezrazmDan(double& P0, double& PN, double& Rc, double& Rk)
{
	Pbegin = P0 = 0;
	Pend = PN = 1;
	Rbegin = Rc = 1e-3;
	Rend = Rk = 1;
}

void Get_Setka(vector<double>& X, int Get)
{
	int N = X.size();

	if(Get==1)
	{
		for (int i = 0; i < N; i++)
		{
			X[i] = Rbegin + i * (1 - Rbegin) / (N - 1.);
		}
	}

	if(Get == 0)
	{
		for (int i = 0; i < N; i++)
		{
			X[i] = Rbegin * pow(Rend/Rbegin, i * 1. / (N - 1.));
		}
	}

	
}

void Solve(vector<double>& B, vector<double>& C, vector<double>& A, vector<double>& D, vector<double>& P, vector<double>& R, int N, double l, double volleft, double volright)
{
	L = l;
	Right = volright;
	Left = volleft;

	double P0, PN, Rc, Rk;
	vvodNach_bezrazmDan(P0, PN, Rc, Rk);
	Get_Setka(R, 0); // 1 - это сетка прямая, 0 - exp

	double R_plus, R_minus;

	for (int i = 1; i < N - 1 ; i++)
	{
		R_plus = 1. / 2 * (R[i + 1] + R[i]);
		R_minus = 1. / 2 * (R[i - 1] + R[i]);

		A[i] = W(R_minus) * f(R_minus) * Pop(R, i - 1)/ (R[i] - R[i - 1]);
		B[i] = W(R_plus) * f(R_plus) * Pop(R, i) / (R[i + 1] - R[i]);
		C[i] = -(A[i] + B[i]);
		D[i] = 0;
	}

	// первая строчка матрицы
	A[0] = 0;		B[0] = 0;		C[0] = 1;		D[0] = P0;		//P[0] = P0;
	// последняя строчка матрицы
	A[N - 1] = 0;	B[N - 1] = 0;	C[N - 1] = 1;	D[N - 1] = PN;	//P[N - 1] = PN;

	progonka_3d(
		B,											// выше диагонали
		C,											// cама диагональ
		A,											// ниже диагонали
		D,											// свободные значение
		P,											// искомый вектор
		N);
}


void resaults(const vector<double>& X, const vector<double>& P, int N)
{
	cout << "l = " << L << "   f(r=Rc)= " << Left << "   f(r=Rg)= " << Right << endl;

	cout << "i" << spaCe << "x[i]    " << spaCe << "f(x)  "<< spaCe << spaCe << spaCe << "P " << spaCe << spaCe << spaCe <<"P an" << endl << endl;
	cout << spaCe1 << spaCe << "U" << spaCe << spaCe << "U_pop    " << spaCe << spaCe << "U an" << endl << endl;
	cout << fixed;

	for (int i = 0; i < N; i++)
	{
		cout << i << spaCe << X[i] << spaCe << f(X[i]) << spaCe << spaCe << P[i] << spaCe << spaCe << P_An(X[i]) << endl << endl;;
		if (i != N - 1)cout << spaCe1 << U_chis_poluuzel(P, X, i) << spaCe << U_chis_poluuzel_pop(P, X, i) << spaCe << spaCe << U_An_poluuzel(X, i) << endl << endl;;
	}
}

double viviod_newyazka(const vector<double>& x, const vector<double>& P, int N)
{
	double S = 0;
	for (int i = 0; i < N; i++)
	{
		S += pow(P_An(x[i]) - P[i], 2);
	}
	S = S / N;
	S = sqrt(S);
	cout << "Невязка(" << N << ") = " << S << endl;
	return S;
}




void vvod_razmernDan(double& L, double& k0, double& m, double& delta_p, double& mu)
{
	L = 100;
	k0 = 1e-12;
	m = 0.2;
	delta_p = 1e6;
	mu = 1e-3;
}

void VvodANDVivod_razmernDan(double& L, double& k0, double& m, double& delta_p, double& mu, const vector<double>& P, const vector<double>& X)
{
	vvod_razmernDan(L, k0, m, delta_p, mu);

	vector<double> U_chisl_bezraz(X.size() - 1);

	double U0 = k0 / mu * delta_p / L;
	double U = 0;
	double V = 0;
	double t = 0;

	cout << fixed;
	for (int i = 0; i < U_chisl_bezraz.size(); i++)
	{
		U_chisl_bezraz[i] = U_chis_poluuzel_pop(P, X, i);
		U = U0 * U_chisl_bezraz[i];
		V = U / m;
		t += (X[i]-X[i + 1]) / V;
		cout << i << " Скор   U0[м/с] = " << U0 << "     Скор фильтр U[м/с] = " << U << "      Скор истин   V[м/с] = " << V << endl;
	}

	cout << "время               t[сек] = " << t*10 << endl
		 << "время               t[сут] = " << (t*10 / (60 * 60 * 24)) << endl;
}


int progonka_3d(const vector<double>& c, const vector<double>& b, const vector<double>& a, const vector<double>& d, vector <double>& X, int n)
{
	//																			расчёт производится по принципуAX = D
	//========================================================================================================================================================================================================================   
	vector <double> Alpha(n, 0);
	vector <double> Betta(n, 0);
	//==========================================================================ПРЯМАЯ ПРОГОНКА
	Alpha[0] = -c[0] * 1. / b[0];
	Betta[0] = d[0] * 1. / b[0];
	for (int i = 1; i < n; i++)
	{
		double y = b[i] + a[i] * Alpha[i - 1];
		double t = d[i] - a[i] * Betta[i - 1];
		Alpha[i] = -c[i] * 1. / y;
		Betta[i] = t / y;
	}
	//==========================================================================ОБРАТНАЯ ПРОГОНКА				
	X[n - 1] = Betta[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		X[i] = Alpha[i] * X[i + 1] + Betta[i];
	}
	return 0;
}
