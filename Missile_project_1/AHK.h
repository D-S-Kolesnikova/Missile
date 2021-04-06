#pragma once
//Полное определение аэродинамических характеристик
#ifndef _ADH_h_
#define _ADH_h_
//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include "Project1.h"
#include "math_lib.h"



//---------------------------------------------------------------------------
const double ToDegreee = 180 / 3.141592654;
//----------------------------------------------------------------------------------
typedef	struct _Aerodynamics_Struct
{
	Vector3	R;			//Равнодействующая аэродинамических сил
	Vector3 M;			//Равнодействующая аэродинамических моментов

	unsigned char N;	//Число компонентов типа double для подсчёта размера структуры
						//Vector3 состоит из 3-х double
						//Matrix33 состоит из 9-ти double
} Aerodynamics_Struct;
//структура входных данных
typedef struct
{
	//параметры для расчета комплекта АДХ в конкретной точке траектории
	double M__;           //мах
	double H__;           //высота полета
	double Alfa__;        //угол атаки
	double n__;           //нерасчетность
	double x_ct__;        //центр тяжести
	double Betta__;       //угол скольжения
	double Delta1__;      //углы поворота каждого из рулей по отдельности
	double Delta2__;      //в системе по тангажу
	double Delta3__;
	double Delta4__;
	double l2_ar__;       //плечо управляющей силы в канале крена
	double l_har__;       //характерная длина в канале крена
	double S_har__;       //характерная площадь
	double L_har__;       //характерная длина в канатах тангажа и рысканья
	double l1_ar__;       //плечо управляющей силы в каналах тангажа и рысканья
	double V__;           //модуль скорости
	double wx__;          //угловые скорости в связанной СК
	double wy__;
	double wz__;
	double q__;           //скоростной напор
	double Px__;			 //тяга двигателя
	char GSS__;           //1 ГСС работает, 0 ГСС выключена
	char B_T;			 //1 боевая ракета, 0 телеметрическая ракета
} ADH_IN_PUT;


typedef struct
{
	double Fxk_;
	double Fyk_;
	double Fzk_;

	double Mxk_;
	double Myk_;
	double Mzk_;

	double Mx_wx_;
	double My_wy_;
	double Mz_wz_;

	double Fx_con_;
	double Fy_con_;
	double Fz_con_;

	double Mx_con_;
	double My_con_;
	double Mz_con_;
	double AlfaSpace;
	double Fi_p;
} DESIGNER;

typedef struct TAeroDH
{
	void SSet(const char* FileAddress);
	double alur(double, double, double, double, double);
	struct KOEFFICIENT Cxp_NOMINAL;
	struct KOEFFICIENT CxpI_NOMINAL;

	struct KOEFFICIENT Delta_Cx_tr_NOMINAL;
	struct KOEFFICIENT Cy_kr_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al00_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al05_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al10_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al15_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al20_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al25_NOMINAL;
	struct KOEFFICIENT Cy_ar_Al30_NOMINAL;
	struct KOEFFICIENT f_ar_NOMINAL;
	struct KOEFFICIENT Cx_don_NOMINAL;
	struct KOEFFICIENT Cxp_H_NOMINAL;
	//	struct KOEFFICIENT Cy_Alfa_al00_NOMINAL;
	struct KOEFFICIENT Mz_wz_Al00_NOMINAL;
	struct KOEFFICIENT Mz_wz_NOMINAL;
	struct KOEFFICIENT Mx_wx_NOMINAL;
	struct KOEFFICIENT Cy_M05_NOMINAL;
	struct KOEFFICIENT Cyp_NOMINAL;
	struct KOEFFICIENT CypI_NOMINAL;
	struct KOEFFICIENT Cyp_H_NOMINAL;
	struct KOEFFICIENT Cd_M05_NOMINAL;
	struct KOEFFICIENT Cdp_NOMINAL;
	struct KOEFFICIENT CdpI_NOMINAL;
	struct KOEFFICIENT Cdp_H_NOMINAL;
	struct KOEFFICIENT Mx_ko_NOMINAL;
	struct KOEFFICIENT dCx_T_NOMINAL;
	struct KOEFFICIENT dCd_T_NOMINAL;

	double Cyp_(double, double, double, double);
	double Cdp_(double, double, double, double, double, char);
	double Cx_(double, double, double, double, double, char);
	double Cy_ar_(double, double, double);
	double Km_(double, double, double);

	DESIGNER KitADH(ADH_IN_PUT);
} AeroDH;

//--------------------------------------------------------------------------------
#endif
