#pragma once
#define _USE_MATH_DEFINES

#include <math.h>


class Missile {
	double ToDegree()
	{
		return 180 / M_PI;
	}
public:
	double SpeedSound(double height);
	double DensityAir(double height);
	double Cx(double M, double alfa_);
	double Cy(double M, double alfa_);
	double Cy_delta(double M, double alfa_);
	double Cz(double M, double alfa_);
	double Cz_delta(double M, double alfa_);
	double Mx_wx(double M, double alfa_);
	double Mz_wz(double M, double alfa_);
	double Mz_alfa(double M, double alfa_);
	double Mz_delta(double M, double alfa_);
	double My_wy(double M, double alfa_);
	double My_betta(double M, double alfa_);
	double My_delta(double M, double betta_);
	double Mxvr(double M, double alfa_);
};

