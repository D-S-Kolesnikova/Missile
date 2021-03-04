// Missile_project_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// Dasha.cpp: определяет точку входа для консольного приложения.
//

#define _USE_MATH_DEFINES

#define	pr1Speed_X Left[0]
#define	pr1Speed_Y Left[1]
#define	pr1Speed_Z Left[2]
#define	pr1Point_X Left[3]
#define	pr1Point_Y Left[4]
#define	pr1Point_Z Left[5]
#define	pr1omegaX   Left[6]
#define	pr1omegaY   Left[7]
#define	pr1omegaZ   Left[8]
#define	pr1ro_rg    Left[9]
#define	pr1lamda_rg Left[10]
#define	pr1mu_rg    Left[11]
#define	pr1nu_rg    Left[12]
#define pr1B		Left[13]
#define pr1L		Left[14]
#define pr1H		Left[15]

#define	Speed_X el[0]
#define	Speed_Y el[1]
#define	Speed_Z el[2]
#define	Point_X el[3]
#define	Point_Y el[4]
#define	Point_Z el[5]
#define	omegaX   el[6]
#define	omegaY   el[7]
#define	omegaZ   el[8]
#define	ro_rg    el[9]
#define	lamda_rg el[10]
#define	mu_rg    el[11]
#define	nu_rg    el[12]
#define	B_    el[13]
#define	L_    el[14]
#define	H_    el[15]

#define TIME Time[0]

//#define ALFA(Vector) -atan2(Vector.Y, Vector.X) 
#define ALFA angle[0]
#define BETTA angle[1]
#define PITCH angle[2]
#define YAW angle[3]
#define ROLL angle[4]
#define TETA angle[5]
#define PSI angle[6]

double* el;
double* Left;
double* Time;
double* angle;

#include <conio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "macros.h"
#include "stdafx.h"
#include "dig.h"
#include "math_lib.h"
#include "Matrix+Vector.h"
#include "Atmosfera.h"
#include "Aerodynamics.h"
#include "Missile.h"
#include "structs.h"
#include "Initialisation_.h"
#include "Print.h"

using namespace std;



void SetAngle()
{
	if (fabs(2 * ro_rg * nu_rg + 2 * lamda_rg * mu_rg) <= 1.0)
		PITCH = asin(2 * ro_rg * nu_rg + 2 * lamda_rg * mu_rg);
	else
	{
		PITCH = asin((2 * ro_rg * nu_rg + 2 * lamda_rg * mu_rg) / fabs(2 * ro_rg * nu_rg + 2 * lamda_rg * mu_rg));
	}

	ROLL = atan2(2 * ro_rg * lamda_rg - 2 * nu_rg * mu_rg, ro_rg * ro_rg + mu_rg * mu_rg - nu_rg * nu_rg - lamda_rg * lamda_rg);
	YAW = atan2(2 * ro_rg * mu_rg - 2 * lamda_rg * nu_rg, ro_rg * ro_rg + lamda_rg * lamda_rg - mu_rg * mu_rg - nu_rg * nu_rg);

	Vector v;
	SET_VECTOR3(v, Speed_X, Speed_Y, Speed_Z);
	double velocity = v.GetLength();
	Matrix A_;
	A_.GetMatrix(ro_rg, lamda_rg, mu_rg, nu_rg);
	v = A_ * v;
	ALFA = -atan2(v.Y, v.X);
	if (fabs(v.Z / velocity) <= 1.0)
		BETTA = asin(v.Z / velocity);
	else
	{
		BETTA = asin((v.Z / velocity) / fabs(v.Z / velocity));
	}

	if (fabs(Speed_Y / velocity) <= 1.0)
		TETA = asin(Speed_Y / velocity);
	else
	{
		TETA = asin((Speed_Y / velocity) / fabs(Speed_Y / velocity));
	}
	PSI = atan2(-Speed_Z, Speed_X);

};




void RightPart(Object Rocket, Initial_Conditions InData, Earth_Struct Earth)
{
	
	FILE* fError;
	Missile Missile1;
	ParamAtmosferaStruct Atm1;
	SetAngle();
	Atm1 = ParamAtmos(Point_Y);
	Vector v(Speed_X, Speed_Y, Speed_Z);
	double velocity = v.GetLength();
	double q = Atm1.Ro * velocity * velocity / 2;
	double mahh = velocity / Atm1.a;
	
	Matrix A_;
	A_.GetMatrix(ro_rg, lamda_rg, mu_rg, nu_rg);
	Matrix A = A_;
	A = ~A; 
	
	double AlfaSpace = sqrt(pow((ALFA), 2) + pow((BETTA), 2));
	Aerodynamic objAerodynamic(mahh, AlfaSpace, ALFA, BETTA);
	Vector F, M, Gg;

	double SquareMiddle = M_PI * pow(Rocket.Dmiddle, 2) / 4;

	F.X = -Missile1.Cx(mahh, AlfaSpace) * q * SquareMiddle;
	F.Y = (Missile1.Cy(mahh, ALFA) * ALFA) * q * SquareMiddle;
	F.Z = (Missile1.Cz(mahh, BETTA) * BETTA) * q * SquareMiddle;

	M.X = (Missile1.Mx_wx(mahh, AlfaSpace) * omegaX * Rocket.LenghtLA/ velocity) * q * SquareMiddle * Rocket.LenghtLA;
	M.Y = (Missile1.My_wy(mahh, BETTA) * omegaY * Rocket.LenghtLA / velocity + Missile1.My_betta(mahh, BETTA) * BETTA) * q * SquareMiddle * Rocket.LenghtLA;
	M.Z = (Missile1.Mz_wz(mahh, ALFA) * omegaZ * Rocket.LenghtLA / velocity + Missile1.Mz_alfa(mahh, ALFA) * ALFA) * q * SquareMiddle * Rocket.LenghtLA;

	EarthModelIni(Earth.A_gd, Earth.E_gd, Earth.G_eq, Earth.Rate);
	Earth_Struct OutStruct = GetEarthParameters(H_, B_);
	Matrix MatrixGS_ST;
	set_GSSK_STSK_matrix33(MatrixGS_ST, InData.B42, InData.L42, PSI, B_, L_);
	//Vector3 Vector1 = mult_matrix33_vector3(MatrixGS_ST, OutStruct.GField);
	Vector Vector1 = MatrixGS_ST * OutStruct.GField;
	Gg = Vector1 * Rocket.mass;

	//Внесены изменения 

	Vector fg;
	//fg = A * F;
	fg = A * (F);
	fg = fg + Gg; //Записали все в нормальной земной
	//  TempPower=fg;

	Vector GS_V = MatrixGS_ST * v;

	pr1Speed_X = fg.X / Rocket.mass;
	pr1Speed_Y = fg.Y / Rocket.mass;
	pr1Speed_Z = fg.Z / Rocket.mass;
	pr1Point_X = Speed_X;
	pr1Point_Y = Speed_Y;
	pr1Point_Z = Speed_Z;
	pr1omegaX = M.X / Rocket.Ix0 - (Rocket.Iz0 - Rocket.Iy0) * omegaY * omegaZ / Rocket.Ix0;
	pr1omegaY = M.Y / Rocket.Iy0 - (Rocket.Ix0 - Rocket.Iz0) * omegaX * omegaZ / Rocket.Iy0;
	pr1omegaZ = M.Z / Rocket.Iz0 - (Rocket.Iy0 - Rocket.Ix0) * omegaX * omegaY / Rocket.Iz0;
	pr1ro_rg = -(omegaX * lamda_rg + omegaY * mu_rg + omegaZ * nu_rg) / 2;
	pr1lamda_rg = (omegaX * ro_rg - omegaY * nu_rg + omegaZ * mu_rg) / 2;
	pr1mu_rg = (omegaX * nu_rg + omegaY * ro_rg - omegaZ * lamda_rg) / 2;
	pr1nu_rg = (-omegaX * mu_rg + omegaY * lamda_rg + omegaZ * ro_rg) / 2;
	pr1B = GS_V.X / (InData.H42 + OutStruct.R_phi);
	pr1L = GS_V.Z / (cos(InData.B42) * (InData.H42 + OutStruct.R_lambda));
	pr1H = GS_V.Y;

	};

void Rks4(int Size, Object Rocket,Initial_Conditions InData, ModelParams_Struct Strct, Earth_Struct Earth)
{

	double a[5];
	double* elTemp;
	double* elh;

	elTemp = new double[Size];
	elh = new double[Size];
	double Step = Strct.Step;
	a[0] = a[1] = a[4] = Step / 2; a[2] = a[3] = Step;

	for (int i = 0; i < Size; i++)
		elTemp[i] = elh[i] = el[i];
	for (int i = 0; i < 4; i++)
	{
		RightPart(Rocket, InData, Earth);

		for (int i = 0; i < Size; i++)
			el[i] = elTemp[i];

		for (int j = 0; j < Size; j++)
		{
			elh[j] += a[i + 1] * Left[j] / 3.0;
			el[j] += a[i] * Left[j];
		}

	}

	if (elh[4] < 0)
	{
		Step = Step * (0 - elTemp[4]) / (elh[4] - elTemp[4]);
		a[0] = a[1] = a[4] = Step / 2; a[2] = a[3] = Step;
		for (int i = 0; i < Size; i++)
		{
			el[i] = elTemp[i];
		}
		for (int i = 0; i < Size; i++) //Цикл для элементов массива
			elh[i] = el[i];
		for (int i = 0; i < 4; i++)
		{
			RightPart(Rocket, InData, Earth);

			for (int j = 0; j < Size; j++)
			{
				elh[j] += a[i + 1] * Left[j] / 3;
				el[j] += a[i] * Left[j];
			}
		}
		for (int i = 0; i < Size; i++)
		{
			el[i] = elh[i];
		};

		double rg = sqrt(el[9] * el[9] + el[10] * el[10] + el[11] * el[11] + el[12] * el[12]);
		double rgh = sqrt(elh[9] * elh[9] + elh[10] * elh[10] + elh[11] * elh[11] + elh[12] * elh[12]);

		for (int j = 0; j < Size; j++)
		{
			el[j] = el[j] / rg;
			elh[j] = elh[j] / rgh;
		}

	}
	else {
		for (int i = 0; i < Size; i++)
		{
			el[i] = elh[i];
		};
		double rg = sqrt(el[9] * el[9] + el[10] * el[10] + el[11] * el[11] + el[12] * el[12]);
		double rgh = sqrt(elh[9] * elh[9] + elh[10] * elh[10] + elh[11] * elh[11] + elh[12] * elh[12]);

		for (int j = 0; j < Size; j++)
		{
			el[j] = el[j] / rg;
			elh[j] = elh[j] / rgh;
		}
	}

	delete[] elTemp;
	delete[] elh;

};


void main(void)
{
	using namespace std;

	el = new double[16];
	Left = new double[16];
	angle = new double[7];
	Time = new double[1];

	Initial_Conditions InData = {};
	Object Rocket = {};
	ModelParams_Struct ParamStr = {};
	Earth_Struct Earth = {};
	
	Initialisation(InData, Rocket, ParamStr, Earth);

	TIME = 0;

	Speed_X = InData.V0 * cos(InData.Pitch0);
	Speed_Y = InData.V0 * sin(InData.Pitch0);
	Speed_Z = 0;

	Point_X = InData.X0;
	Point_Y = InData.Y0;
	Point_Z = InData.Z0;

	omegaX = InData.Wx0;
	omegaY = InData.Wy0;
	omegaZ = InData.Wz0;

	ro_rg = cos(InData.Yaw0 / 2) * cos(InData.Pitch0 / 2) * cos(InData.Roll0 / 2) - sin(InData.Yaw0 / 2) * sin(InData.Pitch0 / 2) * sin(InData.Roll0 / 2);
	lamda_rg = sin(InData.Yaw0 / 2) * sin(InData.Pitch0 / 2) * cos(InData.Roll0 / 2) + cos(InData.Yaw0 / 2) * cos(InData.Pitch0 / 2) * sin(InData.Roll0 / 2);
	mu_rg = sin(InData.Yaw0 / 2) * cos(InData.Pitch0 / 2) * cos(InData.Roll0 / 2) + cos(InData.Yaw0 / 2) * sin(InData.Pitch0 / 2) * sin(InData.Roll0 / 2);
	nu_rg = cos(InData.Yaw0 / 2) * sin(InData.Pitch0 / 2) * cos(InData.Roll0 / 2) - sin(InData.Yaw0 / 2) * cos(InData.Pitch0 / 2) * sin(InData.Roll0 / 2);

	B_ = InData.B42;
	L_ = InData.L42;
	H_ = InData.H42;
	
	SetAngle();

	double j = 0, Step, PrintStep;
	Step = ParamStr.Step;
	PrintStep = ParamStr.PrintStep;

	
		for (TIME; abs(Point_Y) > 1E-8; TIME += Step)
		{
			
			if (abs(TIME - j) < 1E-8)
			{
					j += PrintStep;
					Print_();
			};
		
			Rks4(17, Rocket,InData, ParamStr, Earth);

			if (abs(Point_Y) < 1E-5)
			{
				Print_();		
			}

		}
			
}



// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
