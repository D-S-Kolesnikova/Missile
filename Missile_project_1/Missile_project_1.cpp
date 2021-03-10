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


void main(void)
{
	el = new double[16];
	Left = new double[16];
	angle = new double[9];
	Eiler = new double[6];
	Time = new double[1];


	Initial_Conditions InData = {};
	Object Rocket = {};
	ModelParams_Struct ParamStr = {};
	Earth_Struct Earth = {};
	
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

	DeltaPitch = 0;
	DeltaRoll = 0;
	DeltaYaw = 0;

	double j = 0, Step, PrintStep;
	Step = ParamStr.Step;
	PrintStep = ParamStr.PrintStep;

		for (TIME; Point_Y > 0; TIME += Step)
		{
			if (abs(TIME - j) < 1E-8)
			{
					j += PrintStep;
					Print_();
			};
			Rks4(17, Rocket,InData, ParamStr, Earth);
			if (Point_Y < 0)
			{
				break;
			}
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
