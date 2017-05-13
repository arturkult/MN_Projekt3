#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#define test

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
		//pobieranie danych z pliku dane.csv
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
				xSource.push_back(tmp1);
				ySource.push_back(tmp2);
			}
			plik.close();
		}
#ifndef test
		int iloscWezlow = (xSource.size() / K);
#else
		int iloscWezlow = (xSource.size() / N);
#endif
		int krance = iloscWezlow / PRC_NA_KRANCACH;
		for (int i = 0;i < krance;i++)
		{
			xData.push_back(xSource[i]);
			yData.push_back(ySource[i]);
		}
		for (int i = krance;i < xSource.size() - krance;i++) {
#ifndef test
			if (i%K == 0)
#else
			if (i%N == 0)
#endif
			{
				xData.push_back(xSource[i]);
				yData.push_back(ySource[i]);
			}

		}
		for (int i = 0;i < krance;i++)
		{
			xData.push_back(xSource[xSource.size() - i - 1]);
			yData.push_back(ySource[ySource.size() - i - 1]);
		}

		/*for (int i = 0;i < xSource.size();i++)
		{
			if (i<krance ||
				i>(iloscWezlow - krance) ||
#ifndef test
			i % K/(iloscWezlow-2*krance) == 0
#else
			i%  N == 0
#endif
				)
			{
				xData.push_back(xSource[i]);
				yData.push_back(ySource[i]);
			}
		}*/


		long size = xData.size();
#pragma endregion
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
		double* h = new double[n];
		double* a = new double[n];
		double* al = new double[n];
		double* l = new double[n + 1];
		double* mi = new double[n];
		double* z = new double[n + 1];
		double* c = new double[n + 1];
		double* b = new double[n + 1];
		double* d = new double[n + 1];

		for (int i = 0;i < n;i++)
			h[i] = xData[i + 1] - xData[i];
		for (int i = 0;i < n + 1;i++)
			a[i] = yData[i];
		for (int i = 1;i < n;i++)
			al[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];
		l[0] = 1;
		mi[0] = z[0] = 0;
		for (int i = 1;i < n;i++)
		{
			l[i] = 2 * (xData[i + 1] - xData[i - 1]) - h[i - 1] * mi[i - 1];
			mi[i] = h[i] / l[i];
			z[i] = (al[i] - h[i - 1] * z[i - 1]) / l[i];
		}
		l[n] = 1;
		z[n] = 0;
		c[n] = 0;
		for (int j = n - 1;j >= 0;j--)
		{
			c[j] = z[j] - mi[j] * c[j + 1];
			b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
			d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
		}
		int j = 0;
		x[j] = 0;
		for (int i = 0;i < n;i++)
		{
			for (int k = 0;k < ceil(N / n);j++, k++)
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