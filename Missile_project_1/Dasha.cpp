// Dasha.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "dig.h"
#include <conio.h>
#include <cstdlib>
#include <iostream>
#define _USE_MATH_DEFINES
#include <iomanip>
#include <math.h>
#include "Matrix+Vector.h"
#include "Atmosphere.h"
#include "Aerodynamics.h"
#include <vector>
#include <cmath>

using namespace std;

#define	pr1Speed_gX Left[0]
#define	pr1Speed_gY Left[1]
#define	pr1Speed_gZ Left[2]
#define	pr1Point_gX Left[3]
#define	pr1Point_gY Left[4]
#define	pr1Point_gZ Left[5]
#define	pr1omegaX   Left[6]
#define	pr1omegaY   Left[7]
#define	pr1omegaZ   Left[8]
#define	pr1ro_rg    Left[9]
#define	pr1lamda_rg Left[10]
#define	pr1mu_rg    Left[11]
#define	pr1nu_rg    Left[12]

#define	Speed_gX el[0]
#define	Speed_gY el[1]
#define	Speed_gZ el[2]
#define	Point_gX el[3]
#define	Point_gY el[4]
#define	Point_gZ el[5]
#define	omegaX   el[6]
#define	omegaY   el[7]
#define	omegaZ   el[8]
#define	ro_rg    el[9]
#define	lamda_rg el[10]
#define	mu_rg    el[11]
#define	nu_rg    el[12]

double* el;
double* Left;
//--------------------------------------------------------------
double const
R = 8.31, // газовая постоянная 
c = 5.2E3,// теплоемкость гелия 
g = 9.81,// усокрение свободного падения
ro_He = 0.173, // плотность гелия 
ro_CH4 = 0.717, // плотность метана
P0 = 9.5E4, // давление у Земли
T0 = 273, // начальная температура
ro0 = 1.275, //плотность н.у.
k_atm = 6.5E-3,
n_atm = 5.255;
//--------------------------------------------------------------			
double
/* Требокания к ЛА (Величины представлены в файле начальных условий) */
m_pg, // Масса полезного груза 
L_dal, // Дальность полета
V_kr, // Крейсерская скорость
h_opt, // Высота полета на крейсерском режиме
// Проектные параметры (Величины представлены в файле начальных условий)
ro_ferma, // Условная плотность несущей фермы
k_ferma, //Коэффициент заполнения объема фермы
ro_generator, //Удельная масса генераторов
ro_shell, //Поверхностная площадь оболочки
ro_engine, //Удельная масса электродвигателей
i_du, // Полное число двигателей
ro_frame, // Плотность каркаса оболочки
q1, //Энергоемкость метана
ro_fb, //Плотность топливного бака
q_acc, //Энергоемкость аккумуляторов
nu_ad, // КПД воздушных винтов
nu_generator, //КПД генераторов электроэнергии
delta_T, // Величина нагрева несущего газа
P_el, // Стоимость электроэнергии
P_metan; // Стоимость метана
//--------------------------------------------------------------
double
n_He = -0.369, //коэф плавучести без пг, т.е. подъемная сила на n_нг*100% больше, чем сухой вес аппарата
n_DU = 0.1, //коэф запаса мощности двигательной установки в режиме взлета
n_m = 1.05; //коэф неучтенных факторов по массе 
//--------------------------------------------------------------
/* Параметры, определяемые в ходе решения краевой задачи */
double 
G_He, // вес несущего газа
G_fuel, // вес топлива
G_suh, // вес конструкции
V_b, // объем баллонетов
V_He, // объем несущего газа
V_fuel, // объем топлива на борту
V, // объем дирижабля
Cy, // коэффициент подъемной силы
N; // мощность двигателя
//--------------------------------------------------------------
double V0x = 0, V0y = 0, V0z = 0, Alfa = 0, Betta = 0, Pitch = 0, Yaw = 0, Roll = 0, Mass, ro_m,
Ix = 0, Iy = 0, Iz = 0, Length = 0, DiametrMiddle = 0, L_har,
SquareMiddle = 0, TETA, PSI, Step, Speed, a; bool part;
double
Cx = 0.075, // Коэффициент продольной силы (в СВ СК)
Cp = 1.499; // Коэффициент нормальной силы при вертикальном подъеме

double mz(double alfa)
{
	const int q = 21;
	double Alfa[q] = { -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	double mz[q] = { 0.33088, 0.303, 0.2757, 0.2482, 0.2215, 0.1948, 0.168, 0.14185, 0.11651, 0.09169, 0.0677, 0.04212, 0.1996, -0.00144, -0.02093, -0.04164, -0.06045, -0.07795, -0.09408, -0.10977, -0.12377 };

	int j = 0;
	for (j; j <= 21; j += 1) {
		if (alfa >= Alfa[j] && alfa < Alfa[j + 1]) {
			return (alfa * (mz[j + 1] - mz[j]) / (Alfa[j + 1] - Alfa[j])) + mz[j] -
				Alfa[j] * ((mz[j + 1] - mz[j]) / (Alfa[j + 1] - Alfa[j]));
		}

		if ((alfa >= Alfa[20]) || (alfa <= Alfa[0])) {
			return  -0.019147 * alfa + 0.0677;
		}
	}
}

double Time;
//--------------------------------------------------------------
double ToDegree() //Перевод из радиан в градусы
{
	return 180 / M_PI;
}
//--------------------------------------------------------------
double Volume(double P, double T) // Функция вычисления объема в зависимости от различных внешних условий
{
	return (T * P0 / (P * T0));
}
//--------------------------------------------------------------
double T(double h) //СТАНДАРТНАЯ АТМОСФЕРА (температура)
{
	return T0 - k_atm * h;
}
//--------------------------------------------------------------
double P(double h) //СТАНДАРТНАЯ АТМОСФЕРА (давление)
{
	return P0 * pow(1 - k_atm * h / T0, n_atm);
}
//--------------------------------------------------------------
double ro(double h) //СТАНДАРТНАЯ АТМОСФЕРА (плотность)
{
	return ro0 * pow(1 - k_atm * h / T0, n_atm);
}
//--------------------------------------------------------------
/* Функции правых частей НЛАУ */
void f0(double* x, double*f)
{
	 f[0] = x[0] - ro_He * g * x[4];
	 f[1] = x[1] - ro_CH4 * g * x[5];
	 f[2] = x[6] - x[3] - x[4] - x[5];
	 f[3] = x[3] - (1 - P(h_opt) * T0 / (T(h_opt) * P0)) * x[6];
	 f[4] = ((ro0 / ro_He) + n_He - 1) * x[0] + x[2] * (n_He - 1);
	 f[5] = -(1 + n_DU)* m_pg*g -(1 - ro0/ro_He)*(1+n_DU)*x[0]-(1 + n_DU)*x[2] + x[4]*ro0*g*delta_T * (1 + n_DU) /T0 + 1.84 * i_du * pow(x[8], 2.0 / 3.0);
	 f[6] = m_pg * g + x[0] + x[1] + x[2] - ((ro(h_opt) / ro(0)) * (ro0 / ro_He * x[0] + ro0 / ro_CH4 * x[1] + ro0 * g * x[3]) + x[7] * ro(h_opt) * pow(V_kr, 2.0) * 0.5 * pow(x[6], 2.0 / 3.0));
	 f[7] = x[2] - g * n_m * (Cx * ro(h_opt) * pow(V_kr, 3.0) * pow(x[6], 2.0 / 3.0) * ro_generator / (2 * nu_ad) + i_du * ro_engine * x[8] + ro_ferma * k_ferma * x[6] + (ro_frame + ro_shell) * 6.46 * pow(x[6], 2.0 / 3.0) + 10 * ro_fb * pow(x[5], 2.0 / 3.0));
	 f[8] = x[1] - g / q1 * (2 * x[0] * c * delta_T / g + Cx * ro(h_opt) * pow(V_kr, 2.0) * pow(x[6], 2.0 / 3.0) * L_dal / (2 * nu_generator) + 2 * 1.84 * h_opt * i_du * pow(x[8], 2.0 / 3.0) / nu_generator);
}
//--------------------------------------------------------------
/* Частные производные для матрицы Якоби */
double df6_x8(double x8)
{
	double d = pow(x8, -1.0 / 3.0);
	return 1.84 * i_du * 2.0 * d / 3.0;
}
//--------------------------------------------------------------
double df7_dx6(double x6, double x7)
{
	return -x7 * ro(h_opt) * V_kr * V_kr * 0.5 * 2 / 3 * pow(x6, -1.0 / 3.0);
}
//--------------------------------------------------------------
double df7_dx7(double x6, double x7)
{
	return -ro(h_opt) * V_kr * V_kr * 0.5 * pow(x6, 2.0 / 3.0);
}
//--------------------------------------------------------------
double df8_dx5(double x5)
{
	return -g * n_m * ro_fb * 10 * pow(x5, -1.0 / 3.0) * 2.0 / 3.0;
}
//--------------------------------------------------------------
double df8_dx6(double x6)
{
	double
		a1 = g * n_m * Cx * ro(h_opt) * V_kr * V_kr * V_kr * ro_generator / (2 * nu_ad),
		a2 = g * n_m * ro_ferma * k_ferma,
		a3 = g * n_m * (ro_frame + ro_shell) * 6.46;
	return (-(a1 + a3) * 2 / 3 * pow(x6, -1.0 / 3.0) - a2);
}
//--------------------------------------------------------------
double df9_dx6(double x6)
{
	return g * Cx * ro(h_opt) * L_dal * 2 * V_kr * V_kr * pow(x6, -1.0 / 3.0) / (3 * q1 * 2 * nu_generator);
}
//--------------------------------------------------------------
double df9_dx8(double x8)
{
	double d = pow(x8, -1.0 / 3.0);
	return  - g * 2 * 1.84 * h_opt * i_du * d * 2.0 / (3.0 * q1 * nu_generator);
}
//--------------------------------------------------------------
/* Функция заполнения матрицы Якоби */
void _Jakoby(double** Jakoby, double* x, int n)
{
	Jakoby[0][0] = 1;
	Jakoby[0][4] = -ro_He * g;
	for (int i = 0; i < n; i++)
	{
		if ((i == 0) || (i == 4)) continue;
		Jakoby[0][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[1][1] = 1;
	Jakoby[1][5] = -ro_CH4 * g;
	for (int i = 0; i < n; i++)
	{
		if ((i == 1) || (i == 5)) continue;
		Jakoby[1][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[2][3] = -1;
	Jakoby[2][4] = -1;
	Jakoby[2][5] = -1;
	Jakoby[2][6] = 1;
	for (int i = 0; i < n; i++)
	{
		if ((i == 3) || (i == 4) || (i == 5) || (i == 6)) continue;
		Jakoby[2][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[3][3] = 1;
	Jakoby[3][6] = -(1 - P(h_opt) * T0 / (T(h_opt) * P0));
	for (int i = 0; i < n; i++)
	{
		if ((i == 3) || (i == 6)) continue;
		Jakoby[3][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[4][0] = (ro0 / ro_He) + n_He - 1;
	Jakoby[4][2] = n_He - 1;
	for (int i = 0; i < n; i++)
	{
		if ((i == 0) || (i == 2)) continue;
		Jakoby[4][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[5][0] = -(1 - ro0 / ro_He) * (1 + n_DU);
	Jakoby[5][2] = -(1 + n_DU);
	Jakoby[5][4] = ro0 * g * delta_T * (1 + n_DU) / T0;
	Jakoby[5][8] = df6_x8(x[8]);
	for (int i = 0; i < n; i++)
	{
		if ((i == 0) || (i == 2) || (i == 4) || (i == 8)) continue;
		Jakoby[5][i] = 0;
	}
	//--------------------------------------------------------------
	Jakoby[6][0] = (1 - ro(h_opt) * ro0 / (ro(0) * ro_He));
	Jakoby[6][1] = (1 - ro(h_opt) * ro0 / (ro(0) * ro_CH4));
	Jakoby[6][2] = 1;
	Jakoby[6][3] = -ro(h_opt) * g * ro0 / (ro(0));
	Jakoby[6][4] = 0;
	Jakoby[6][5] = 0;
	Jakoby[6][6] = df7_dx6(x[6], x[7]);
	Jakoby[6][7] = df7_dx7(x[6], x[7]);
	Jakoby[6][8] = 0;
	//--------------------------------------------------------------
	Jakoby[7][0] = 0;
	Jakoby[7][1] = 0;
	Jakoby[7][2] = 1;
	Jakoby[7][3] = 0;
	Jakoby[7][4] = 0;
	Jakoby[7][5] = df8_dx5(x[5]);
	Jakoby[7][6] = df8_dx6(x[6]);
	Jakoby[7][7] = 0;
	Jakoby[7][8] = -g * n_m * i_du * ro_engine;
	//--------------------------------------------------------------
	Jakoby[8][0] = -2 * c * delta_T / q1;
	Jakoby[8][1] = 1;
	Jakoby[8][2] = 0;
	Jakoby[8][3] = 0;
	Jakoby[8][4] = 0;
	Jakoby[8][5] = 0;
	Jakoby[8][6] = df9_dx6(x[6]);
	Jakoby[8][7] = 0;
	Jakoby[8][8] = df9_dx8(x[8]);
}
//--------------------------------------------------------------
/* Функции печати матрицы NxN */
void PrintMatr(double** mas, int m) {
	int i, j;
	cout << "------------------------------------------------------------------------------------------------------" << endl;
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++)
			cout << setw(10) << setprecision(4) << mas[i][j] << " ";
		cout << endl;
	}
	cout << "------------------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------
void PrintVector(double* vector, int n)
{
	cout << "---------------------------" << endl;
	for (int j = 0; j < n; j++) cout << setw(10) << setprecision(4) << vector[j] << endl;
	cout << "---------------------------" << endl;
};
//--------------------------------------------------------------
/* Получение матрицы без i-й строки и j-го столбца */
void GetMatr(double** mas, double** p, int i, int j, int m) {
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki < m - 1; ki++) { // проверка индекса строки
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) { // проверка индекса столбца
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}
//--------------------------------------------------------------
/* Рекурсивное вычисление определителя */
double Determinant(double** mas, int m) {
	int i, j, k, n;
	double d;
	double** p;
	p = new double* [m];
	for (i = 0; i < m; i++)
		p[i] = new double[m];
	j = 0; d = 0;
	k = 1; //(-1) в степени i
	n = m - 1;
	if (m < 0.00001) cout << "Определитель вычислить невозможно!";
	if (m == 1) {
		d = mas[0][0];
		return(d);
	}
	if (m == 2) {
		d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
		return(d);
	}
	if (m > 2) {
		for (i = 0; i < m; i++) {
			GetMatr(mas, p, i, 0, m);
			//cout << mas[i][j] << endl;
			//PrintMatr(p, n);
			d = d + k * mas[i][0] * Determinant(p, n);
			k = -k;
		}
	}
	return(d);
}
//--------------------------------------------------------------
/* Функция вычеркивания строки и столбца*/
void Get_matr(double** matr, int n, double** temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}
//--------------------------------------------------------------
/* Функция освобождения буфера */
void Free(double** matr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] matr[i];
	delete[] matr;
}
//--------------------------------------------------------------
/* Функция нахождения обратной матрицы */
void obr_matr(double** matr, double** obr_matr, int n)
{
	double Det = Determinant(matr, n);
	//--------------------------------------------------------------
	/* Нахождение матрицы алгебраических доолнений */
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int m = n - 1;
			double** temp_matr = new double* [m];
			for (int k = 0; k < m; k++)
				temp_matr[k] = new double[m];
			Get_matr(matr, n, temp_matr, i, j);
			obr_matr[i][j] = pow(-1.0, i + j + 2) * Determinant(temp_matr, m) / Det;
			Free(temp_matr, m);
		}
	}
	/* Транспонирование матрицы алгебраических дополнений */
	double memory;
	for (int i = 0; i < n - 1; i++)
		for (int j = i + 1; j < n; j++)
		{
			memory = obr_matr[i][j];
			obr_matr[i][j] = obr_matr[j][i];
			obr_matr[j][i] = memory;
		}
}
//--------------------------------------------------------------
/* Нулевая инициализация матрицы */
void InitNull_matr(double** matr, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) matr[i][j] = 0.0;
};
//--------------------------------------------------------------
/* Нулевая инициализация вектора */
void InitNull_vector(double* vector, int n)
{
		for (int j = 0; j < n; j++) vector[j] = 0.0;
};
//--------------------------------------------------------------
/* Перемножение матрицы на вектор */
void Multiplication(double** matr, double* x, int n, double* temp_vector)
{
	InitNull_vector(temp_vector, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			temp_vector[i] += matr[i][j] * x[j];
		}
	}

}
//--------------------------------------------------------------
/* Суммирование векторов */
void sum_vector(double* x, double* y, double* sum, int n)
{
	for (int i = 0; i < n; i++)
		sum[i] = x[i] + y[i];
}
//--------------------------------------------------------------
/* Вычитание векторов */
void difference_vector(double* x, double* y, double* sum, int n)
{
	for (int i = 0; i < n; i++)
		sum[i] = x[i] - y[i];
}
//--------------------------------------------------------------
/* Вычисление углов Эйлера через парамерты Р-Г */
void SetAngle(double r, double l, double m, double n)
{
	FILE* fError;

	if (fabs(2 * r * n + 2 * l * m) <= 1.0)
		Pitch = asin(2 * r * n + 2 * l * m);
	else
	{
		Pitch = asin((2 * r * n + 2 * l * m) / fabs(2 * r * n + 2 * l * m));
	}

	Roll = atan2(2 * r * l - 2 * n * m, r * r + m * m - n * n - l * l);
	Yaw = atan2(2 * r * m - 2 * l * n, r * r + l * l - m * m - n * n);
};
//--------------------------------------------------------------
/* Функция знака */
double sign(double a)
{
	if ((a > -0.00000001) && (a < 0.00000001))
		return 0.0;
	return a / abs(a);
}
/* Блок расчета векторов правых частей */
void RightPart(double Time)
{
	FILE* fError;
	double h_geopotencial, Temprature, DavlenieAir = 0, DensityAir = 0, SpeedSound = 0, M_t, q;



	SetAngle(ro_rg, lamda_rg, mu_rg, nu_rg);
	Atm objAtm(Point_gY);
	objAtm.Calculation(Point_gY, DavlenieAir, DensityAir, SpeedSound);
	Vector vg(Speed_gX, Speed_gY, Speed_gZ);
	double velocity = vg.GetLength();
	//q = DensityAir * velocity * velocity / 2;
	q = ro(Point_gY) * velocity * velocity / 2;
		double mahh = velocity / SpeedSound;
	Matrix A_;//Матрица перехода от норм. земн. с.к. к связанной
	A_.GetMatrix(ro_rg, lamda_rg, mu_rg, nu_rg);

	Matrix A = A_;//Матрица перехода от связанной с.к. к норм. земн.
	A = ~A; // М-ца перехода от связанной к норм - это транспонированная от норм к связанной
	Vector v = A_ * vg;//Скорость в связанной системе коррдинат 

	Alfa = -atan2(v.Y, v.X);

	if (fabs(v.Z / velocity) <= 1.0)
		Betta = asin(v.Z / velocity);
	else
	{
		Betta = asin((v.Z / velocity) / fabs(v.Z / velocity));
	}

	double AlfaSpace = sqrt(pow((Alfa), 2) + pow((Betta), 2));
	Aerodynamic objAerodynamic(mahh, AlfaSpace, Alfa, Betta);

	Vector F, M, Gg, P, B, Pr_m;

	if (part)
	{
		F.Y = - sign(Speed_gY) * Cp * q * SquareMiddle;
		F.X = 0; 

		M.X = 0;
		M.Y = 0;
		M.Z = 0;

		P.X = 0;
		P.Y = 1.84 * i_du * pow(N, 2.0 / 3.0);
	}
	else
	{
		F.X =  - Cx * ro(Point_gY) * pow(Speed_gX, 2.0) * 0.5 * SquareMiddle;
		F.Y = Cy * ro(Point_gY) * pow(Speed_gX, 2.0) * 0.5 * SquareMiddle - sign(Speed_gY) * Cp * ro(Point_gY) * pow(Speed_gY, 2.0) * 0.5 * SquareMiddle;;
		
		M.X = 0;
		M.Y = 0;
		M.Z = 0;
		//M.Z = mz(Alfa) * q * SquareMiddle * Length;

		P.X = Cx * ro(Point_gY) * pow(V_kr, 2.0) * 0.5 * SquareMiddle;
		P.Y = 0;
	}
	

	Gg.X = 0;
	Gg.Y = - (G_suh + G_He + G_fuel + m_pg * g);
	Gg.Z = 0;

	B.X = 0;
	B.Y = ro(Point_gY) / ro(0) * (ro0/ro_He * G_He + ro0 / ro_CH4 * G_fuel + ro0 * g * V_b);
	B.Z = 0;

	Vector fg;
	//fg = A * F;
	fg = A * (F + P);
	fg = fg + Gg + B; //Записали все в нормальной земной
	//  TempPower=fg;
	
	if (fabs(Speed_gY / velocity) <= 1.0)
		 TETA = asin(Speed_gY / velocity);
	else
	{
	 TETA = asin((Speed_gY / velocity) / fabs(Speed_gY / velocity));
	}

	 PSI = atan2(-Speed_gZ, Speed_gX);

	pr1Speed_gX = fg.X / Mass;
	pr1Speed_gY = fg.Y / Mass;
	pr1Speed_gZ = fg.Z / Mass;
	pr1Point_gX = Speed_gX;
	pr1Point_gY = Speed_gY;
	pr1Point_gZ = Speed_gZ;
	pr1omegaX = M.X / Ix - (Iz - Iy) * omegaY * omegaZ / Ix;
	pr1omegaY = M.Y / Iy - (Ix - Iz) * omegaX * omegaZ / Iy;
	pr1omegaZ = M.Z / Iz - (Iy - Ix) * omegaX * omegaY / Iz;
	pr1ro_rg = -(omegaX * lamda_rg + omegaY * mu_rg + omegaZ * nu_rg) / 2;
	pr1lamda_rg = (omegaX * ro_rg - omegaY * nu_rg + omegaZ * mu_rg) / 2;
	pr1mu_rg = (omegaX * nu_rg + omegaY * ro_rg - omegaZ * lamda_rg) / 2;
	pr1nu_rg = (-omegaX * mu_rg + omegaY * lamda_rg + omegaZ * ro_rg) / 2;
};

void Rks4(int Size, double Time)
{

	double a[5];
	double* elTemp;
	double* elh;

	elTemp = new double[Size];
	elh = new double[Size];
	a[0] = a[1] = a[4] = Step / 2; a[2] = a[3] = Step;

	for (int i = 0; i < Size; i++)
		elTemp[i] = elh[i] = el[i];
	for (int i = 0; i < 4; i++)
	{
		RightPart(Time);
		for (int i = 0; i < Size; i++)
			el[i] = elTemp[i];

		for (int j = 0; j < Size; j++)
		{
			elh[j] += a[i + 1] * Left[j] / 3.0;
			el[j] += a[i] * Left[j];
		}

	}

		for (int i = 0; i < Size; i++)
		{
			el[i] = elh[i];
		};
	
	delete[] elTemp;
	delete[] elh;

};

/* Основной блок программы */
void main(void)
{
	setlocale(LC_ALL, "rus");
	el = new double[14];
	Left = new double[14];
	FILE *fError;
    fError=fopen("fError.err","w");

/*Начальные условия фигурирующие в интерфейсе, но напрямую неиспользующиеся в программном обеспечении*/
	FILE *FDATE;
	if ((FDATE=fopen("dasha.txt","r"))==NULL) return;
	char cc[256];
	int num=sizeof(cc)/sizeof(char);
	while (fgets(cc,num,FDATE)!=NULL)
	{
/*Считывание требований к ЛА*/
	if (strstr(cc, "m_pg") != NULL)m_pg = StrFileToDouble(cc, fError);
	if (strstr(cc, "L_dal") != NULL)L_dal = StrFileToDouble(cc, fError);
	if (strstr(cc, "V_kr") != NULL)V_kr = StrFileToDouble(cc, fError);
	if (strstr(cc, "h_opt") != NULL)h_opt = StrFileToDouble(cc, fError);
/*Считывание проектных параметров*/
	if (strstr(cc, "ro_ferma") != NULL)ro_ferma = StrFileToDouble(cc, fError);
	if (strstr(cc, "k_ferma") != NULL)k_ferma = StrFileToDouble(cc, fError);
	if (strstr(cc, "ro_generator") != NULL)ro_generator = StrFileToDouble(cc, fError);
	if (strstr(cc, "ro_shell") != NULL)ro_shell = StrFileToDouble(cc, fError);
	if (strstr(cc, "ro_engine") != NULL)ro_engine = StrFileToDouble(cc, fError);
	if (strstr(cc, "i_du") != NULL)i_du = StrFileToDouble(cc, fError);
	if (strstr(cc, "ro_frame") != NULL)ro_frame = StrFileToDouble(cc, fError);
	if (strstr(cc, "q1") != NULL)q1 = StrFileToDouble(cc, fError);
	if (strstr(cc, "ro_fb") != NULL)ro_fb = StrFileToDouble(cc, fError);
	if (strstr(cc, "q_acc") != NULL)q_acc = StrFileToDouble(cc, fError);
	if (strstr(cc, "nu_ad") != NULL)nu_ad = StrFileToDouble(cc, fError);
	if (strstr(cc, "nu_generator") != NULL)nu_generator = StrFileToDouble(cc, fError);
	if (strstr(cc, "delta_T") != NULL)delta_T = StrFileToDouble(cc, fError);
	if (strstr(cc, "P_el") != NULL)P_el = StrFileToDouble(cc, fError);
	if (strstr(cc, "P_metan") != NULL)P_metan = StrFileToDouble(cc, fError);
/*Считывание начальных данных для расчета траектории ЛА*/
	if (strstr(cc, "X0") != NULL)Point_gX = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Y0") != NULL)Point_gY = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Z0") != NULL)Point_gZ = StrFileToDouble(cc, fError);//el
	if (strstr(cc,"V0")!=NULL)Speed = StrFileToDouble(cc,fError);//el
	if (strstr(cc, "Alfa0") != NULL) Alfa = StrFileToDouble(cc, fError) / ToDegree();
	if (strstr(cc, "Betta0") != NULL) Betta = StrFileToDouble(cc, fError) / ToDegree();
	if (strstr(cc, "Pitch0") != NULL) Pitch = StrFileToDouble(cc, fError) / ToDegree();
	if (strstr(cc, "Yaw0") != NULL) Yaw = StrFileToDouble(cc, fError) / ToDegree();
	if (strstr(cc, "Roll0") != NULL) Roll = StrFileToDouble(cc, fError) / ToDegree();
	if (strstr(cc, "Wx0") != NULL) omegaX = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Wy0") != NULL) omegaY = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Wz0") != NULL) omegaZ = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Ix") != NULL) Ix = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Iy") != NULL) Iy = StrFileToDouble(cc, fError);//el
	if (strstr(cc, "Iz") != NULL) Iz = StrFileToDouble(cc, fError);//el
	} /*конец блока считывания данных*/
//--------------------------------------------------------------
/* Блок решения системы уравнений для определения проектных параметров */
	double** a, * y, *x, *f;
	int n = 9;
	double** Jakoby = new double* [n];
	double** obr_Jakoby = new double* [n];

	for (int i = 0; i < n; i++) 
	{
		Jakoby[i] = new double[n];
		obr_Jakoby[i] = new double[n];
	}
	x = new double[n];
	a = new double* [n];
	y = new double[n];
	f = new double[n];
//--------------------------------------------------------------
/* Начальные приближения */
	x[0] = 1E5;
	x[1] = 1E4;
	x[2] = 1E6;
	x[3] = 1E5;
	x[4] = 1E5;
	x[5] = 1E5;
	x[6] = 1E5;
	x[7] = 1;
	x[8] = 1E5;
//--------------------------------------------------------------
	double* x1, *y1;
	x1 = new double[n];
	y1 = new double[n];
	int stop = 0;
//--------------------------------------------------------------
	
/* Алгоритм решения краевой задачи */
	do {
		f0(x, f);// Вектор функций от НУ
		_Jakoby(Jakoby, x, n);// Заполнение матрицы Якоби 9х9 
		cout << "-----------------------------------------Матрица Якоби------------------------------------------------" << endl;
		PrintMatr(Jakoby, n);
		obr_matr(Jakoby, obr_Jakoby, n);
		cout << "-------------------------------------Обратная матрица Якоби-------------------------------------------" << endl;
		PrintMatr(obr_Jakoby, n);
		Multiplication(obr_Jakoby, f, n, y1); // Умножение обратной матрицы Якоби на F(x0), где х0 - вектор НУ
		//PrintVector(y1, n);
		difference_vector(x, y1, x1, n); // х - J(-1)*F(x0)
		PrintVector(x1, n);
		difference_vector(x, x1, y, n); // Вектор поправок
		PrintVector(y, n);
		for (int i = 0; i < n; i++)
		{
			x[i] = x1[i];
		}
		stop++;
	} while (stop < 3);

		/*G_He = 1.87271E5,
		G_fuel = 1.25536E5,
		G_suh = 8.20893E5,
		V_b = 5.14091E4,
		V_He = 1.10346E5,
		V_fuel = 1.784762E4,
		V = 1.79602E5,
		Cy = 0.36527,
		N = 6.79495E5;*/

		G_He = x[0],
		G_fuel = x[1],
		G_suh = x[2],
		V_b = x[3],
		V_He = x[4],
		V_fuel = x[5],
		V = x[6],
		Cy = x[7],
		N = x[8];
/* Конец блока решения системы уравнений */
//--------------------------------------------------------------
/* Начало блока интегрирования уравнений движения */
	Speed_gX = Speed * cos(Pitch);
	Speed_gY = Speed * sin(Pitch);
	Speed_gZ = 0;
	SquareMiddle = pow(V, 2.0 / 3.0);
	Length = pow(V, 1.0 / 3.0);

	Mass = (G_suh + G_He + G_fuel + m_pg * g) / g;
	
	ro_rg = cos(Yaw / 2) * cos(Pitch / 2) * cos(Roll / 2) - sin(Yaw / 2) * sin(Pitch / 2) * sin(Roll / 2);
	lamda_rg = sin(Yaw / 2) * sin(Pitch / 2) * cos(Roll / 2) + cos(Yaw / 2) * cos(Pitch / 2) * sin(Roll / 2);
	mu_rg = sin(Yaw / 2) * cos(Pitch / 2) * cos(Roll / 2) + cos(Yaw / 2) * sin(Pitch / 2) * sin(Roll / 2);
	nu_rg = cos(Yaw / 2) * sin(Pitch / 2) * cos(Roll / 2) - sin(Yaw / 2) * cos(Pitch / 2) * sin(Roll / 2);
	
	cout << ro_rg << " - ro_rg" << endl
		<< lamda_rg << " - lyambda_rg" << endl
		<< mu_rg << " - mu " << endl
		<< nu_rg << " - nu" << endl;

	ofstream file;
	file.open("file_result.txt");
	file << setw(20) << "Time" << " "
		<< setw(20) << "Mass" << " "
		<< setw(20) << "Point_gX" << " "
		<< setw(20) << "Point_gY" << " "
		<< setw(20) << "Point_gZ" << " "
		<< setw(20) << "Speed_gX" << " "
		<< setw(20) << "Speed_gY" << " "
		<< setw(20) << "Speed_gZ" << " " 
		<< setw(20) << "Pitch" << " "
		<< setw(20) << "Yaw" << " "
		<< setw(20) << "Roll" << " "
		<< setw(20) << "TETA" << " "
		<< setw(20) << "Alfa" << " "
		<< setw(20) << "Betta" << " "
		<< setw(20) << "Omegax" << " "
		<< setw(20) << "Omegay" << " "
		<< setw(20) << "Omegaz" << " "
		<< endl;
	part = 1;
	Time = 0;
	Step = 0.1;
	int j = 0;
	do 
	{
		if (abs(Time - j) < 0.01)
		{
			file << setw(20) << Time << " "
				<< setw(20) << Mass << " "
				<< setw(20) << Point_gX << " "
				<< setw(20) << Point_gY << " "
				<< setw(20) << Point_gZ << " "
				<< setw(20) << Speed_gX << " "
				<< setw(20) << Speed_gY << " "
				<< setw(20) << Speed_gZ << " "
				<< setw(20) << Pitch << " "
				<< setw(20) << Yaw << " "
				<< setw(20) << Roll << " "
				<< setw(20) << TETA << " "
				<< setw(20) << Alfa << " "
				<< setw(20) << Betta << " "
				<< setw(20) << omegaX << " "
				<< setw(20) << omegaY << " "
				<< setw(20) << omegaZ << " "
				<< endl;
			j++;
		}
		Rks4(14, Time);
		Time += Step;

	} while (Time < 1000);
	part = 0;
	Step = 0.1;
	do
	{
		if (abs(Time - j) < 0.01)
		{
			file << setw(20) << Time << " "
				<< setw(20) << Mass << " "
				<< setw(20) << Point_gX << " "
				<< setw(20) << Point_gY << " "
				<< setw(20) << Point_gZ << " "
				<< setw(20) << Speed_gX << " "
				<< setw(20) << Speed_gY << " "
				<< setw(20) << Speed_gZ << " "
				<< setw(20) << Pitch << " "
				<< setw(20) << Yaw << " "
				<< setw(20) << Roll << " "
				<< setw(20) << TETA << " "
				<< setw(20) << Alfa << " "
				<< setw(20) << Betta << " "
				<< setw(20) << omegaX << " "
				<< setw(20) << omegaY << " "
				<< setw(20) << omegaZ << " "
				<< endl;
			j++;
		}
		Rks4(14, Time);
		Time += Step;

	} while (Time < 10000);


	
	

delete[] el;
delete[] Left;



	file.close();
	

}

