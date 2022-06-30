
#include "Header.h"


int main()
{
	setlocale(LC_ALL, "Russian");

	double l = 0.5, volleft = 1, volright = 1;
	int N;

	vvodN(N);

	vector<double> B(N, 0);											// выше диагонали
	vector<double> C(N, 0);											// cама диагональ
	vector<double> A(N, 0);											// ниже диагонали
	vector<double> D(N, 0);											// свободные значение
	vector<double> P(N, 0);											// искомый вектор
	vector<double> X(N, 0);


	Solve(
		B,											// выше диагонали
		C,											// cама диагональ
		A,											// ниже диагонали
		D,											// свободные значение
		P,											// искомый вектор
		X,											// сетка
		N,											// точки
		l,
		volleft,
		volright
	);
	resaults(X, P, N);

	//cout << Q_chis(P, X, 0)<<endl << Q_chis(P, X, N - 2)<<endl << Q_An();
	/*
	ofstream File("1.txt");
	for (int i = 0; i < X.size(); i++)
	{
		File << X[i] << " " << P[i] << endl;
	}
	File.close();
	*/
	/*
	ofstream File("1.txt");
	for (int N = 2; N <= 1000; N++)
	{
		vector<double> B(N, 0);											// выше диагонали
		vector<double> C(N, 0);											// cама диагональ
		vector<double> A(N, 0);											// ниже диагонали
		vector<double> D(N, 0);											// свободные значение
		vector<double> P(N, 0);											// искомый вектор
		vector<double> X(N, 0);

		Solve(
			B,											// выше диагонали
			C,											// cама диагональ
			A,											// ниже диагонали
			D,											// свободные значение
			P,											// искомый вектор
			X,											// сетка
			N,											// точки
			l, 
			volleft, 
			volright
		);

		File << N << " " << viviod_newyazka(X,P,N) << endl;
	}
	File.close();
	*/
	
	double L, k0, m, delta_p, mu;
	VvodANDVivod_razmernDan(L, k0, m, delta_p, mu, P, X);
	
}

