#pragma once
#include "structs.h"
//#include "EarthModel(ПЗ90).h"
#include "EarthModel.h"
#include "stdafx.h"
#include "dig.h"
#define _USE_MATH_DEFINES
#include "stdafx.h"
#include "dig.h"
#include <conio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>



double ToDegree()
{
	return 180 / M_PI;
}


void Initialisation(Initial_Conditions& Rocket1, Object& Rocket, ModelParams_Struct& ParamStr, Earth_Struct& Earth)
{
	FILE* fError;
	fError = fopen("fError.err", "w");
	FILE* FDATE;
	if ((FDATE = fopen("dasha.txt", "r")) == NULL) return;

	char cc[256];
	int num = sizeof(cc) / sizeof(char);

	while (fgets(cc, num, FDATE) != NULL)
	{
		/*<Начальное положение и начальная скорость ЛА>*/
		if (strstr(cc, "X0") != NULL)Rocket1.X0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Y0") != NULL)Rocket1.Y0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Z0") != NULL)Rocket1.Z0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "V0") != NULL)Rocket1.V0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "X_1") != NULL)Rocket1.X_1 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Y_1") != NULL)Rocket1.Y_1 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "H_1") != NULL)Rocket1.H_1 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "X_target") != NULL)Rocket1.X_target = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Y_target") != NULL)Rocket1.Y_target = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Z_target") != NULL)Rocket1.Z_target = StrFileToDouble(cc, fError);//el


		/*<Начальные угловое положение и угловые скорости относительно связанных осей координат>*/
		if (strstr(cc, "Alfa0") != NULL) Rocket1.Alfa0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Betta0") != NULL) Rocket1.Betta0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Pitch0") != NULL) Rocket1.Pitch0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Yaw0") != NULL) Rocket1.Yaw0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Roll0") != NULL) Rocket1.Roll0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Wx0") != NULL) Rocket1.Wx0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Wy0") != NULL) Rocket1.Wy0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Wz0") != NULL) Rocket1.Wz0 = StrFileToDouble(cc, fError);//el


		/*Массово-инерционные и геометрические характеристики*/
		if (strstr(cc, "m0") != NULL) Rocket.mass = StrFileToDouble(cc, fError);
		if (strstr(cc, "Ix0") != NULL) Rocket.Ix0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "Iy0") != NULL) Rocket.Iy0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "Iz0") != NULL) Rocket.Iz0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "LengthLA") != NULL) Rocket.LenghtLA = StrFileToDouble(cc, fError);
		if (strstr(cc, "Dm") != NULL) Rocket.Dmiddle = StrFileToDouble(cc, fError);
		if (strstr(cc, "l1_ar") != NULL) Rocket.l1_ar = StrFileToDouble(cc, fError);
		if (strstr(cc, "l2_ar") != NULL) Rocket.l2_ar = StrFileToDouble(cc, fError);
		if (strstr(cc, "l_har") != NULL) Rocket.l_har = StrFileToDouble(cc, fError);
		if (strstr(cc, "L_har") != NULL) Rocket.L_har = StrFileToDouble(cc, fError);
		if (strstr(cc, "x_cm") != NULL) Rocket.x_cm = StrFileToDouble(cc, fError);

		/*Параметры моделирования*/
		if (strstr(cc, "Step") != NULL) ParamStr.Step = StrFileToDouble(cc, fError);
		if (strstr(cc, "PrintData") != NULL) ParamStr.PrintStep = StrFileToDouble(cc, fError);

		//Исходные данные модели Земли("Эллипсоид Красовского)
		if (strstr(cc, "A_GD_1") != NULL) Earth.A_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "E_GD_1") != NULL)  Earth.E_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "G_EQUATOR_1") != NULL)  Earth.G_eq = StrFileToDouble(cc, fError);
		if (strstr(cc, "EARTH_RATE_1") != NULL) Earth.Rate = StrFileToDouble(cc, fError);

		/*Исходные данные модели Земли("ПЗ90)*/
		/*if (strstr(cc, "A_GD_2") != NULL) Earth.A_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "B_GD_2") != NULL) Earth.B_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "E_GD_2") != NULL)  Earth.E_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "G_EQUATOR_2") != NULL)  Earth.G_eq = StrFileToDouble(cc, fError);
		if (strstr(cc, "EARTH_RATE_2") != NULL) Earth.Rate = StrFileToDouble(cc, fError);
		if (strstr(cc, "fM") != NULL) Earth.fM = StrFileToDouble(cc, fError);
		if (strstr(cc, "R_sr") != NULL) Earth.R_sr = StrFileToDouble(cc, fError);
		if (strstr(cc, "alfa") != NULL) Earth.alpha = StrFileToDouble(cc, fError);
		if (strstr(cc, "gamma_a") != NULL) Earth.gamma_a = StrFileToDouble(cc, fError);
		if (strstr(cc, "gamma_b") != NULL) Earth.gamma_b = StrFileToDouble(cc, fError);*/

	}
}