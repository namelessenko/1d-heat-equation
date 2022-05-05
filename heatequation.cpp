#include <math.h>
#include <iostream>
#include <stdio.h>
#include <cmath>

using namespace std;

#define N 30

const float T0(20), Tl(0), Tr(0);
const float L(0.1);
const float lamda = 46, ro = 7800, c = 460;
const float alpha = lamda / (ro * c);
float T[N], TT[N];
float Ta[N], TTa[N];


double h = L / (N - 1);
float tao = 0.25 * pow(h, 2) / alpha;
float timer = 0.0;
float t_end = 60;

int main(int argc, char **argv) {
	setlocale(LC_ALL, "rus");

	cout << "Данные:" << endl;
	cout << "Длина = " << L << endl;
	cout << "Коэф. теплопроводности = " << alpha << endl;
	cout << "Начальное распределение темпереатуры = " << T0 << endl;
	cout << "Температура на левой границе = " << Tl << endl;
	cout << "Температура на правой границе = " << Tr << endl;
	cout << "Шаг по времени = " << tao << endl;
	cout << "Шаг по сетке = " << h << endl;
	cout << "Конечное время = " << t_end << endl;


	cout << "ЧИСЛЕННОЕ РЕШЕНИЕ" << endl;

	for (int i = 1; i < N; i++) {
		T[i] = T0;
	}

	T[0] = Tl;
	T[N - 1] = Tr;

	while (timer < t_end) {
		timer += tao;
		for (int i = 0; i < N; i++)
			TT[i] = T[i];
		for (int i = 1; i < N - 1; i++)
			T[i] = TT[i] + ((lamda * pow(tao, 2)) / (ro * c)) * ((TT[i + 1] - 2.0 * TT[i] + TT[i - 1]) / pow(h, 2));
	}	
	
	cout << endl;
	for (int i = 1; i < N-1; i++) {
		cout <<"В узле " << i << ": " << T[i] << endl;
	}
	cout << endl;

	cout << "РЕШЕНИЕ В ПРИГРАНИЧНЫХ ОБАСТЯХ: " << endl;
	double left = (2 * (lamda / 2) * T[14] + h * 20) / (2 * (lamda / 2) + h);
	cout << "x = L: " << left << endl;
	double right = (2 * (lamda / 2) * T[0] + h * 20) / (2 * (lamda / 2) + h);
	cout << "x = 0: " << right << endl;
	cout << endl;

	cout << "АНАЛИТИЧЕСКОЕ РЕШЕНИЕ" << endl;
    cout << "ЧЕТНАЯ" << endl;
	float x = 0;
	double sum = 0;
	double R = 0.05;
	for (int i = 0; i <= 5; i++) {
		for (int m = 1; m < 10000; m++) {
			double mju = (2 * m - 1)*(3.14 / 2);
			double Dm = (2 * T0 / mju)*pow(-1, m + 1);
			sum = sum + (Dm * cos(mju*(x / R)) * exp(((-1)*(mju*mju)) * ((alpha*t_end) / (R*R))));
		}
		cout << "В точке " << x << ": " <<  sum << endl;
		sum = 0;
		x += 0.01;
    }
    cout << endl;
    
    cout << "НЕЧЕТНАЯ" << endl;
    x = 0;
    sum = 0;
    R = 0.05;
    for (int i = 0; i <= 5; i++) {
        for (int m = 1; m < 10000; m++) {
            double mju = (2 * m - 1)*(3.14 / 2);
            double Dm = (2 * T0 / mju)*pow(-1, m + 1);
            sum = sum + (Dm * sin(mju*(x / R)) * exp(((-1)*(mju*mju)) * ((alpha*t_end) / (R*R))));
        }
        cout << "В точке " << x << ": " <<  sum << endl;
        sum = 0;
        x += 0.01;
    }
}
