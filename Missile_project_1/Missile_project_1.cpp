// Missile_project_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// Dasha.cpp: определяет точку входа для консольного приложения.
//

#define _USE_MATH_DEFINES

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
#include "Integration.h"
#include "Print.h"

using namespace std;
 
void Calculation_hit(Initial_Conditions InData, Object Rocket, ModelParams_Struct ParamStr, Earth_Struct Earth, AeroDH& ADH_KIT)
{

	Gss_Kr_to_SK42(InData.X_1, InData.Y_1, B42, L42);

	B_ = B42;
	L_ = L42;
	H_ = InData.H_1;

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


	SetAngle();
	PSI0 = PSI;

	DeltaPitch = 0;
	DeltaRoll = 0;
	DeltaYaw = 0;

	double j = 0, Step, PrintStep;
	Step = ParamStr.Step;
	PrintStep = ParamStr.PrintStep;

	for (TIME; Point_Y > 1E-5; TIME += Step)
	{

		if (abs(TIME - j) < 1E-8)
		{
			j += PrintStep;
			//cout << " ШАГ ЗАПИСАН В ФАЙЛ - " << j << endl;
			Print_();
		};
		Rks4(17, Rocket, InData, ParamStr, Earth, ADH_KIT);
	}
		Print_();
		MISS = sqrt(_SQR(Point_X - target_X) + _SQR(Point_Y - target_Y) + _SQR(Point_Z - target_Z));		
}


void main(void)
{
	setlocale(LC_ALL, "Russian");

	el = new double[16];
	Left = new double[16];
	angle = new double[10];
	Eiler = new double[6];
	Time = new double[1];
	target = new double[10];
	G = new double[3];
	delta = new double[4];

	Initial_Conditions InData = {};
	Object Rocket = {};
	ModelParams_Struct ParamStr = {};
	Earth_Struct Earth = {};
	DESIGNER ADH;
	AeroDH ADH_KIT;
	ADH_IN_PUT InPutADH;
	ParamAtmosferaStruct Atm1;
	
	Initialisation(InData, Rocket, ParamStr, Earth);
	Gss_Kr_to_SK42(InData.X_1, InData.Y_1, B42, L42);

	B_ = B42;
	L_ = L42;
	H_ = InData.H_1;

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


	SetAngle();
	PSI0 = PSI;

	DeltaPitch = 0;
	DeltaRoll = 0;
	DeltaYaw = 0;

	Vector V_ST(Speed_X, Speed_Y, Speed_Z);
	Atm1 = ParamAtmos(Point_Y);
	double velocity = V_ST.GetLength();
	double q = Atm1.Ro * velocity * velocity / 2;
	double mahh = velocity / Atm1.a;

	InPutADH.M__ = mahh;
	InPutADH.Alfa__ = ALFA;
	InPutADH.Betta__ = BETTA;
	InPutADH.H__ = Point_Y;

	InPutADH.Delta1__ = -delta_1;
	InPutADH.Delta2__ = -delta_2;
	InPutADH.Delta3__ = delta_3;
	InPutADH.Delta4__ = delta_4;

	InPutADH.n__ = 0;
	InPutADH.x_ct__ = Rocket.x_cm;
	InPutADH.V__ = velocity;

	InPutADH.S_har__ = (M_PI * pow(Rocket.Dmiddle,2)) / 4.0;
	InPutADH.l1_ar__ = Rocket.l1_ar;
	InPutADH.l2_ar__ = Rocket.l2_ar;
	InPutADH.l_har__ = Rocket.l_har;
	InPutADH.L_har__ = Rocket.L_har;

	InPutADH.wx__ = omegaX;
	InPutADH.wy__ = omegaY;
	InPutADH.wz__ = omegaZ;
	InPutADH.q__ = q;
	InPutADH.Px__ = 0;

	InPutADH.B_T = 1;

	
	const char* filename = "AHK.txt";
	ADH_KIT.SSet(filename);
	Calculation_hit(InData, Rocket, ParamStr, Earth, ADH_KIT);

	target_X = InData.X_target;
	target_Y = InData.Y_target;
	target_Z = InData.Z_target;
		
	


	/*
	Calculation_hit(InData, Rocket, ParamStr, Earth);
	double Step = 1;
	do {
		do {
			target_X += Step;
			Calculation_hit(InData, Rocket, ParamStr, Earth);
		} while (MISS < 5);
		target_X -= Step;
		Calculation_hit(InData, Rocket, ParamStr, Earth);
		Step = Step / 10.0;
	} while (abs(MISS - 5) > 1E-5);
	              
	double temp1 = target_X;
	double temp2 = target_Y;
	double temp3 = target_Z;
	double temp_miss;
	Calculation_hit(InData, Rocket, ParamStr, Earth);
	int n = 0;
	do {
		double fi = 0, r = 10, step_fi = 10;
		do {
			temp_miss = MISS;
			do {
				fi += step_fi * TO_RAD;
				target_X = r * cos(fi) + temp1;
				target_Z = r * sin(fi) + temp3;
				Calculation_hit(InData, Rocket, ParamStr, Earth);

			} while (MISS > 5);
			target_X -= r * cos(fi);
			target_Z -= r * sin(fi);
			fi -= step_fi * TO_RAD;
			step_fi = step_fi / 10.0;
			

		} while (abs(MISS - temp_miss) > 10E-5);

		target_X += r * cos(fi);
		target_Z += r * sin(fi);
		fi += step_fi * 10 * TO_RAD;
		n++;
		if (n > 60) {
			break;
		}
		 temp1 = target_X;
		 temp3 = target_Z;
		 Print_zone();
		 cout << "Точка под номером - " << n <<endl;
	} while (n < 100);*/

	cout << "Программа выполнена" << endl;

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
