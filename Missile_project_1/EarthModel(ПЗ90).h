/*#pragma once
#include "Matrix+Vector.h"
#include "EarthModel.h"

typedef struct {
	double A_gd;
	double B_gd;
	double G_eq;
	double Rate;
	double E_gd;
	double R_sr;
	double fM;
	double alpha;
	double gamma_a;
	double gamma_b;
	double R_phi;
	double R_lambda;
	double R_gd;
	double phi_gr;
	Vector AcclCp_GS;
	Vector GField;
}  Earth_Struct;

void Gss_Kr_to_SK42(double x, double y, double& B42_, double& L42_);
void matr_sph_grin(const double grin_coord[], double sph_coord[]);
void matr_grin_sph(const double sph_coord[], double grin_coord[]);
void matr_grin_geod(const double geod_coord[], double grin_coord[]);
void SF_GR(Vector GR, Vector& SF);
void GR_GDZ(Vector GDZ, Vector& GR);
void GR_SF(Vector SF, Vector& GR);
double delt_g_anom(double B, double L, double H);
double g_norm(double B, double H);
void	EarthModelIni(double a_gd, double b_gd, double e_gd, double g_eq, double e_rate, double R_sr, double fM_, double alpha_, double gamma_a_, double gamma_b_);
Earth_Struct GetEarthParameters(double H, double B, double L);  //Функция расчёта параметров модели Земли
*/
