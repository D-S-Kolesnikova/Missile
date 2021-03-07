
#include "macros.h"
#include "math_lib.h"
#include "Matrix+Vector.h"
#include "Atmosfera.h"
#include "Missile.h"
#include "structs.h"
#include "Initialisation_.h"

double* el;
double* Left;
double* Time;
double* angle;

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
	Vector F, M, Gg;

	double SquareMiddle = M_PI * pow(Rocket.Dmiddle, 2) / 4;

	F.X = -Missile1.Cx(mahh, AlfaSpace) * q * SquareMiddle;
	F.Y = (Missile1.Cy(mahh, ALFA) * ALFA) * q * SquareMiddle;
	F.Z = (Missile1.Cz(mahh, BETTA) * BETTA) * q * SquareMiddle;

	M.X = (Missile1.Mx_wx(mahh, AlfaSpace) * omegaX * Rocket.LenghtLA / velocity) * q * SquareMiddle * Rocket.LenghtLA;
	M.Y = (Missile1.My_wy(mahh, BETTA) * omegaY * Rocket.LenghtLA / velocity + Missile1.My_betta(mahh, BETTA) * BETTA) * q * SquareMiddle * Rocket.LenghtLA;
	M.Z = (Missile1.Mz_wz(mahh, ALFA) * omegaZ * Rocket.LenghtLA / velocity + Missile1.Mz_alfa(mahh, ALFA) * ALFA) * q * SquareMiddle * Rocket.LenghtLA;

	EarthModelIni(Earth.A_gd, Earth.E_gd, Earth.G_eq, Earth.Rate);
	Earth_Struct OutStruct = GetEarthParameters(H_, B_);
	/*Matrix33 MatrixGS_ST;
		set_GSSK_STSK_matrix33(MatrixGS_ST, InData.B42, InData.L42, PSI, B_, L_);
		Vector3 Vector1 = mult_matrix33_vector3(MatrixGS_ST, OutStruct.GField);
	 - Закомментировала 07.03 перешла на определения через Vec_Matr из класса*/

	Matrix GS_ST;
	GS_ST.GSSK_STSK(B42, L42, PSI, B_, L_);
	Vector G_GS;
	G_GS = OutStruct.GField;

	Vector V_ST = GS_ST * G_GS;
	Gg = Rocket.mass * V_ST;

	Vector fg;
	fg = A * F;
	fg = fg + Gg;

	/*	Vector3 Speed = { Speed_X, Speed_Y, Speed_Z };
		Vector3 GS_V = mult_matrix33_vector3(MatrixGS_ST, Speed);
	- Закомментировала 07.03 перешла на определения через Vec_Matr из класса*/

	Vector GS_V = GS_ST * v;

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
	pr1B = GS_V.X / (H_ + OutStruct.R_phi);
	pr1L = GS_V.Z / (cos(B_) * (H_ + OutStruct.R_lambda));
	pr1H = GS_V.Y;

};

void Rks4(int Size, Object Rocket, Initial_Conditions InData, ModelParams_Struct Strct, Earth_Struct Earth)
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
