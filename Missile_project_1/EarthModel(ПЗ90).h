#pragma once
#include "Matrix+Vector.h"


typedef struct {
	double A_gd;
	double G_eq;
	double Rate;
	double E_gd;
	double R_phi;
	double R_lambda;
	double R_gd;
	double phi_gr;
	Vector AcclCp_GS;
	Vector GField;
}  Earth_Struct;

void Gss_Kr_to_SK42(double x, double y, double& B42_, double& L42_);
void			EarthModelIni(double a_gd, double e_gd, double g_eq, double e_rate, double R_sr, double fM_, double alpha_, double gamma_a_, double gamma_b_);
Earth_Struct	GetEarthParameters(double H, double B, double L);  //Функция расчёта параметров модели Земли