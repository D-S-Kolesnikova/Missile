
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
double* Eiler;

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

/*------------------Расчет матриц направляющих косинусов-----------------------*/
// Матрица перехода от географической СК к стартовой СК
	Matrix GS_ST;
	GS_ST.GSSK_STSK(B42, L42, PSI, B_, L_);
// Матрица перехода от стартовой СК к географической СК
	Matrix ST_GS = ~GS_ST;
//Матрица перехода от связанной СК к норм земной
	Matrix A_;
	A_.GetMatrix(ro_rg, lamda_rg, mu_rg, nu_rg);
//Матрица перехода от норм земной к стартовой СК
	Matrix A = A_;
	A = ~A;
//Матрица перехода от географического трехгранника к свзанному
	Matrix GS_SV;
	GS_SV.GSSK_SVSK(YAW, PITCH, ROLL);
/*------------------Расчет линейных и угловых скоростей-----------------------*/
//Линейная скорость в стартовой СК
	Vector V_ST(Speed_X, Speed_Y, Speed_Z);
//Линейная скорость в географической СК
	Vector V_GS = GS_ST * V_ST;
//Угловая скорость вращения СК42 в проекциях на оси географической СК
	Vector RateEarth_GS(Earth.Rate * cos(B_), Earth.Rate * sin(B_), 0);
//Угловая скорость вращения географической СК относительно осей СК42 в проекциях на оси стартовой СК
	Vector RateEarth_ST = GS_ST * RateEarth_GS;
//Угловая скорость вращения географической СК относительно осей СК42 в проекциях на оси связанной СК
	Vector RateEarth_SV = GS_SV * RateEarth_GS;
//Угловая скорость вращения географической СК относительно осей СК42 в проекциях на оси географической СК
	Vector RateGS_GS(V_GS.Z / (H_ + Earth.R_lambda), V_GS.Z * tan(B_) / (H_ + Earth.R_lambda), -V_GS.X / (H_ + Earth.R_phi));
//Угловая скорость вращения осей связанной СК относительно осей географической СК в проекциях на оси связанной СК	
	Vector RateGS_SV = GS_SV * RateGS_GS;
//Угловая скорость вращения Земли в проекциях на оси связанной СК
	Vector RateSV_SV(omegaX, omegaY, omegaZ);
//Абсолютная угловая скорость ЛА в проекциях на связанной СК
	Vector RateAbs_SV = RateSV_SV + RateEarth_SV + RateGS_SV;
//Функция расчета углов
	SetAngle();

	Atm1 = ParamAtmos(Point_Y);
	double velocity = V_ST.GetLength();
	double q = Atm1.Ro * velocity * velocity / 2;
	double mahh = velocity / Atm1.a;
	double AlfaSpace = sqrt(pow((ALFA), 2) + pow((BETTA), 2));
	double SquareMiddle = M_PI * pow(Rocket.Dmiddle, 2) / 4;

	double K_pitch1 = 2;
	double K_yaw1 = 5;
	double K_roll1 = 3;
	double K_roll2 = 1;

	pr1Pitch = omegaY * sin(ROLL) + omegaZ * cos(ROLL);
	pr1Yaw = (omegaY * cos(ROLL) - omegaZ * sin(ROLL)) / cos(PITCH);
	pr1Roll = omegaX - tan(PITCH) * (omegaY * cos(ROLL) - omegaZ * sin(ROLL));

	DeltaPitch = -K_pitch1 * pr1Pitch;
	DeltaYaw = -K_yaw1 * pr1Yaw;
	DeltaRoll = -K_roll1 * ROLL - K_roll2 * pr1Roll;


/*-------Расчет внешних сил и моментов, действующих на ЛА (кроме силы тяжести), в проекциях на СвСК----*/
	Vector F, M, Gg;
	
	F.X = -Missile1.Cx(mahh, AlfaSpace) * q * SquareMiddle;
	F.Y = (Missile1.Cy(mahh, ALFA) * ALFA + Missile1.Cy_delta(mahh, ALFA) * DeltaPitch) * q * SquareMiddle;
	F.Z = (Missile1.Cz(mahh, BETTA) * BETTA + Missile1.Cz_delta(mahh,ALFA) * DeltaYaw) * q * SquareMiddle;

	M.X = (Missile1.Mx_wx(mahh, AlfaSpace) * omegaX * Rocket.LenghtLA / velocity) * q * SquareMiddle * Rocket.LenghtLA;
	M.Y = (Missile1.My_wy(mahh, BETTA) * omegaY * Rocket.LenghtLA / velocity + Missile1.My_betta(mahh, BETTA) * BETTA + Missile1.My_delta(mahh, BETTA) * DeltaYaw) * q * SquareMiddle * Rocket.LenghtLA;
	M.Z = (Missile1.Mz_wz(mahh, ALFA) * omegaZ * Rocket.LenghtLA / velocity + Missile1.Mz_alfa(mahh, ALFA) * ALFA + Missile1.Mz_delta(mahh, ALFA) * DeltaPitch) * q * SquareMiddle * Rocket.LenghtLA;
	
//Инициализация числовых костанст для гравитационной модели
	EarthModelIni(Earth.A_gd, Earth.E_gd, Earth.G_eq, Earth.Rate);
//Расчет параметров гравитационной модели
	Earth_Struct OutStruct = GetEarthParameters(H_, B_);
//Проекция ускорения силы тяжести на оси географической СК
	Vector G_GS;
	G_GS = OutStruct.GField;
//Проекция ускорения силы тяжести на оси стартовой СК
	Vector G_ST = GS_ST * G_GS;
//Проекция силы тяжести на оси стартовой СК
	Gg = Rocket.mass * G_ST;
//Кажущееся ускорение в проекциях на оси связанной СК	
	Vector fg;
	fg = A * F;
	fg = fg + Gg;
//Ускорение Кориолиса в проекциях на оси стартовой СК
	Vector Accl_Kor = 2.0 * (RateEarth_ST % V_ST);
//Абсолютное ускорение в проекциях на оси стартовой СК
	Vector F_ST = fg * pow(Rocket.mass, -1.) - Accl_Kor;
//Тензор инерции на оси связанной СК
	Matrix J_SV;
	J_SV.TestGetMatrix(Rocket.Ix0, 0, 0, 0, Rocket.Iy0, 0, 0, 0, Rocket.Iz0);
//Производная тензора инерции на оси связанной СК
	Matrix pr1_J_SV;
	pr1_J_SV.InitNull();
//Обратный тензор инерции на оси связанной СК
	Matrix J_SV_inverse = J_SV;
	J_SV_inverse.Inverse_Matrix();
//Динамическое уравнение Эйлера
	Vector pr1_RateAbs_SV;
	pr1_RateAbs_SV = J_SV_inverse * (M - pr1_J_SV * RateAbs_SV - RateAbs_SV % (J_SV * RateAbs_SV));

	
	pr1Speed_X = F_ST.X;
	pr1Speed_Y = F_ST.Y;
	pr1Speed_Z = F_ST.Z;
	pr1Point_X = Speed_X;
	pr1Point_Y = Speed_Y;
	pr1Point_Z = Speed_Z;
	pr1omegaX = pr1_RateAbs_SV.X;
	pr1omegaY = pr1_RateAbs_SV.Y;
	pr1omegaZ = pr1_RateAbs_SV.Z;
	//pr1omegaX = M.X / Rocket.Ix0 - (Rocket.Iz0 - Rocket.Iy0) * omegaY * omegaZ / Rocket.Ix0;
	//pr1omegaY = M.Y / Rocket.Iy0 - (Rocket.Ix0 - Rocket.Iz0) * omegaX * omegaZ / Rocket.Iy0;
	//pr1omegaZ = M.Z / Rocket.Iz0 - (Rocket.Iy0 - Rocket.Ix0) * omegaX * omegaY / Rocket.Iz0;
	pr1ro_rg = -(omegaX * lamda_rg + omegaY * mu_rg + omegaZ * nu_rg) / 2;
	pr1lamda_rg = (omegaX * ro_rg - omegaY * nu_rg + omegaZ * mu_rg) / 2;
	pr1mu_rg = (omegaX * nu_rg + omegaY * ro_rg - omegaZ * lamda_rg) / 2;
	pr1nu_rg = (-omegaX * mu_rg + omegaY * lamda_rg + omegaZ * ro_rg) / 2;
	pr1B = V_GS.X / (H_ + OutStruct.R_phi);
	pr1L = V_GS.Z / (cos(B_) * (H_ + OutStruct.R_lambda));
	pr1H = V_GS.Y;
	
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
