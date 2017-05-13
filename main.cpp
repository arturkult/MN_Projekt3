#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
//#define test

using namespace std;

#define N 100 //liczba wezlow wyjsciowych 
#define PRC_NA_KRANCACH 3 // 1/PRC_NA_KRANCACH wszystkich wezlow wyjsciowych na krancach przedzialow

int main() {
	double x[N];
	double Wm[N];
#ifndef test
	for (int K = 20;K < 100;K += 10)
	{
#endif
		vector<double> xData;
		vector<double> yData;
		vector<double> xSource;
		vector<double> ySource;
		ifstream plik;
		plik.open("dane.csv");
#pragma region PobieranieDanych
		int ileLinii = 0;
		if (plik.is_open())
		{
			double tmp1, tmp2;
			char buff[100];
			char delimiter;
			plik.getline(buff, 100);
			while (!plik.eof())
			{
				plik >> tmp1;
				plik.get(delimiter);
				plik >> tmp2;
#ifdef test
				if (ileLinii % 10 == 0)
					//if(ileLinii<ILOSC_NA_KRANCACH || ileLinii>N-ILOSC_NA_KRANCACH  || ileLinii % (N - 2 * ILOSC_NA_KRANCACH) == 0)

#else
				if (ileLinii % K == 0)
#endif
				{
					xData.push_back(tmp1);
					yData.push_back(tmp2);
				}
				ileLinii++;

			}
			plik.close();
		}
		long size = xData.size();
#pragma endregion
		//pobieranie danych z pliku dane.csv
//#pragma region PobieranieDanych
//		int ileLinii = 0;
//		if (plik.is_open())
//		{
//			double tmp1, tmp2;
//			char buff[100];
//			char delimiter;
//			plik.getline(buff, 100);
//			while (!plik.eof())
//			{
//				plik >> tmp1;
//				plik.get(delimiter);
//				plik >> tmp2;
//				xSource.push_back(tmp1);
//				ySource.push_back(tmp2);
//			}
//			plik.close();
//		}
//#ifndef test
//		int iloscWezlow = (xSource.size() / K);
//#else
//		int iloscWezlow = (xSource.size() / N);
//#endif
//		int krance = iloscWezlow / PRC_NA_KRANCACH;
//		for (int i = 0;i < krance;i++)
//		{
//			xData.push_back(xSource[i]);
//			yData.push_back(ySource[i]);
//		}
//		for (int i = krance;i < xSource.size() - krance;i++) {
//#ifndef test
//			if (i%K == 0)
//#else
//			if (i%N == 0)
//#endif
//			{
//				xData.push_back(xSource[i]);
//				yData.push_back(ySource[i]);
//			}
//
//		}
//		for (int i = 0;i < krance;i++)
//		{
//			xData.push_back(xSource[xSource.size() - i - 1]);
//			yData.push_back(ySource[ySource.size() - i - 1]);
//		}
//
//
//		long size = xData.size();
//#pragma endregion
		/*
#pragma region AproksymacjaLagrange

		x[0] = 0;
		for (int k = 0;k < N;k++)
		{
			if (k > 0)
				x[k] = x[k - 1] + xData.at(size - 1) / N;
			double suma = 0;
			for (int i = 0;i < size;i++)
			{
				double iloczyn = 1;
				for (int j = 0;j < size;j++)
					if (i != j)
						iloczyn *= (x[k] - xData.at(j)) / (xData.at(i) - xData.at(j));
				suma += iloczyn*yData.at(i);
			}
			Wm[k] = suma;

		}
		ofstream outplik;
		outplik.precision(16);
		outplik.open("wyniki.txt",ios::app);
		outplik << "K\t" << size << endl;
		for (int i = 0;i < N;i++) {
			outplik << x[i] << "\t" << Wm[i] << endl;
		}
		outplik.close();
#pragma endregion
	*/

#pragma region AproksymacjaSplajn
		int n = size - 1;
		double* h = new double[n + 1];
		double* mi = new double[n + 1];
		double* del = new double[n + 1];
		double* lam = new double[n + 1];
		double* ci = new double[n + 1];
		double* di = new double[n + 1];
		double* M = new double[n + 1];

		double* a = mi;
		double* b = new double[n + 1];
		double* c = lam;
		double* d = del;
		double* cp = new double[n + 1];
		double* dp = new double[n + 1];

		lam[0] = del[0] = mi[n] = del[n] = 0;

		for (int i = 0;i < n;i++)
			h[i + 1] = xData[i + 1] - xData[i];

		for (int i = 1;i < n;i++)
		{
			mi[i] = h[i] / (h[i] + h[i + 1]);
			lam[i] = h[i + 1] / (h[i] + h[i + 1]);
			del[i] = (6 / (h[i] + h[i + 1]))*(((yData[i + 1] - yData[i]) / h[i + 1]) - ((yData[i] - yData[i - 1]) / h[i]));
		}
		for (int i = 0;i < n+1;i++)
			b[i] = 2;
		cp[0] = c[0] / b[0];
		dp[0] = d[0] / b[0];

		for (int i = 1;i < n;i++)
		{
			cp[i] = c[i] / (b[i] - a[i] * cp[i - 1]);
		}

		for(int i=1;i<n+1;i++)
			dp[i] = (d[i] - a[i] * dp[i - 1]) / (b[i] - a[i] * cp[i - 1]);

		M[n] = dp[n];

		for (int i = n - 1;i >= 0;i--)
		{
			M[i] = dp[i] - cp[i] * M[i + 1];
		}

		for (int i = 0;i < n;i++)
		{
			a[i] = yData[i];
			b[i] = (yData[i + 1] - yData[i]) / h[i + 1] - (2 * M[i] + M[i + 1])*h[i + 1] / 6;
			c[i] = M[i] / 2;
			d[i] = (M[i + 1] - M[i]) / (6 * h[i + 1]);
		}

		int j = 0;
		int ileWPrzedziale = ceil(double(N) / n);
		x[j] = 0;
		for (int i = 0;i < n;i++)
		{
			for (int k = 0;k < ileWPrzedziale && j<N;j++, k++)
			{
				if (j > 0)
					x[j] = x[j - 1] + xData.at(size - 1) / N;
				Wm[j] = a[i] + b[i] * (x[j] - xData[i]) + c[i] * pow(x[j] - xData[i], 2) + d[i] * pow(x[j] - xData[i], 3);
			}
		}

		ofstream outplik;
		outplik.precision(16);
#ifdef test		
		outplik.open("wynikiSplajn2.txt");
#else
		outplik.open("wynikiSplajn2.txt", ios::app);
#endif
		outplik << "K\t" << size << endl;
		for (int i = 0;i < N;i++) {
			outplik << x[i] << "\t" << Wm[i] << endl;
		}
		outplik.close();
#pragma endregion
#ifndef test
}
#endif


	return 0;
}