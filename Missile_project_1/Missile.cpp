#include <math.h>

#include "Missile.h"


double Missile::SpeedSound(double height)
{
	if (height <= 0.0)
		return 340.0;

	else if (height > 0.0 && height <= 11000.0)
		return sqrt(299.79e-9 * height * height - 2.6135 * height + 115.80e+3);

	else if (height > 11000.0 && height <= 25000.0)
		return 295.07;

	else if (height > 25000.0 && height <= 46000.0)
		return sqrt(-179.80e-9 * height * height + 1.1098 * height + 59.433e+3);

	else
		return 331.82;
};

double Missile::DensityAir(double height)
{
	if (height <= 0.0)
		return 1.225;

	else if (height > 0.0 && height <= 13000.0)
		return 3.1088e-9 * height * height - 114.04e-6 * height + 1.223;

	else if (height > 13000.0 && height <= 26000.0)
		return log(-298.95e-12 * height * height + (12.780e+3) / height + 28.677e-6 * height);

	else if (height > 26000.0 && height <= 39000.0)
		return exp(-7.9310 * log(height) - 94590 / height + 80.891);

	else if (height > 39000.0 && height <= 49000.0)
		return exp(1.4632e-9 * height * height + 264.69e-3 * log(height) - 266.74e-6 * height);

	else
		return 0.0;
};

double Missile::Cx(double M, double alfa_)
{
	return 1.0 / (73.211 / exp(M) - 47.483 / M + 16.878);
	//return 0.3;
};

double Missile::Cy(double M, double alfa_)
{
	double Ds = 1.86 * (11.554 / exp(M) - 2.5191e-3 * M * M - 5.024 / M + 52.836e-3 * M + 4.112);
	if (Ds >= 0)
		return sqrt(Ds);
	else
		return 1.86 * 1.039;
};

double Missile::Cy_delta(double M, double alfa_)
{
	double p1 = 0.0, p2 = 0.0, p3 = 0.0, p4 = 0.0;
	alfa_ = fabs(alfa_ * ToDegree());
	p1 = 1 / (243.84e-3 / exp(alfa_) + 74.309e-3);
	p2 = log(1.9773 * alfa_ * alfa_ - 25.587 * alfa_ + 83.354);
	p3 = 18.985 * alfa_ * alfa_ - 375.76 * alfa_ + 1471;
	p4 = -51.164e-3 * alfa_ * alfa_ + 805.52e-3 * alfa_ + 1.8929;
	return (-p1 * 1e-6 * M * M + p2 * 1e-12 * exp(M) - p3 * 1e-6 * M - p4 * 1e-3) * 2;
};

double Missile::Cz(double M, double alfa_)
{
	return -Cy(M, alfa_);
};


double Missile::Cz_delta(double M, double alfa_)
{
	return -Cy_delta(M, alfa_);
};

double Missile::Mx_wx(double M, double alfa_)
{
	return -0.005 * 0.6786;
};



double Missile::Mz_wz(double M, double alfa_)
{

	return 1.89 * (146.79e-6 * M * M - 158.98e-3 / M - 7.6396e-3 * M - 68.195e-3);
};

double Missile::Mz_alfa(double M, double alfa_)
{
	return (-766.79e-3 / exp(M) + 438.74e-3 / M + 5.8822e-3 * M - 158.34e-3);
};

double Missile::Mz_delta(double M, double alfa_)
{
	double k1 = 0.0, k2 = 0.0;
	alfa_ = fabs(alfa_ * ToDegree());
	k1 = exp(-19.488e-3 * alfa_ * alfa_ - 378.62e-3 * alfa_ + 6.7518);
	k2 = exp(-21.234e-3 * alfa_ * alfa_ - 635.84e-6 * exp(alfa_) - 98.296e-3 * alfa_ + 2.5938);
	return 1.89 * sqrt(k1 * 1e-9 * M * M + k2 * 1e-6);
};

double Missile::My_wy(double M, double alfa_)
{
	return Mz_wz(M, alfa_);
};

double Missile::My_betta(double M, double alfa_)
{
	return Missile::Mz_alfa(M, alfa_);
};

double Missile::My_delta(double M, double betta_)
{
	return Mz_delta(M, betta_);
};


double Missile::Mxvr(double M, double alfa_)
{

	return 0.0;

};