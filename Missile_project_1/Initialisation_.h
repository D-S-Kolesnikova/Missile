#pragma once
#include "structs.h"
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
		/*<��������� ��������� � ��������� �������� ��>*/
		if (strstr(cc, "X0") != NULL)Rocket1.X0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Y0") != NULL)Rocket1.Y0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Z0") != NULL)Rocket1.Z0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "V0") != NULL)Rocket1.V0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "B42") != NULL)Rocket1.B42 = StrFileToDouble(cc, fError) / ToDegree();//el
		if (strstr(cc, "H42") != NULL)Rocket1.H42 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "L42") != NULL)Rocket1.L42 = StrFileToDouble(cc, fError) / ToDegree();//el

		/*<��������� ������� ��������� � ������� �������� ������������ ��������� ���� ���������>*/
		if (strstr(cc, "Alfa0") != NULL) Rocket1.Alfa0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Betta0") != NULL) Rocket1.Betta0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Pitch0") != NULL) Rocket1.Pitch0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Yaw0") != NULL) Rocket1.Yaw0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Roll0") != NULL) Rocket1.Roll0 = StrFileToDouble(cc, fError) / ToDegree();
		if (strstr(cc, "Wx0") != NULL) Rocket1.Wx0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Wy0") != NULL) Rocket1.Wy0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "Wz0") != NULL) Rocket1.Wz0 = StrFileToDouble(cc, fError);//el
		if (strstr(cc, "A0") != NULL) Rocket1.A0 = StrFileToDouble(cc, fError);

		/*�������-����������� � �������������� ��������������*/
		if (strstr(cc, "m0") != NULL) Rocket.mass = StrFileToDouble(cc, fError);
		if (strstr(cc, "Ix0") != NULL) Rocket.Ix0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "Iy0") != NULL) Rocket.Iy0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "Iz0") != NULL) Rocket.Iz0 = StrFileToDouble(cc, fError);
		if (strstr(cc, "LengthLA") != NULL) Rocket.LenghtLA = StrFileToDouble(cc, fError);
		if (strstr(cc, "Dm") != NULL) Rocket.Dmiddle = StrFileToDouble(cc, fError);

		/*��������� �������������*/
		if (strstr(cc, "Step") != NULL) ParamStr.Step = StrFileToDouble(cc, fError);
		if (strstr(cc, "PrintData") != NULL) ParamStr.PrintStep = StrFileToDouble(cc, fError);

		/*�������� ������ ������ �����*/
		if (strstr(cc, "A_GD") != NULL) Earth.A_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "E_GD") != NULL)  Earth.E_gd = StrFileToDouble(cc, fError);
		if (strstr(cc, "G_EQUATOR") != NULL)  Earth.G_eq = StrFileToDouble(cc, fError);
		if (strstr(cc, "EARTH_RATE") != NULL) Earth.Rate = StrFileToDouble(cc, fError);
	}
}