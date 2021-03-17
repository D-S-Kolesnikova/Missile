
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

//Инициализация числовых костанст для гравитационной модели
	EarthModelIni(Earth.A_gd, Earth.E_gd, Earth.G_eq, Earth.Rate);
//Расчет параметров гравитационной модели
	Earth_Struct OutStruct = GetEarthParameters(H_, B_);
/*------------------Расчет матриц направляющих косинусов-----------------------*/
// Матрица перехода от географической СК к стартовой СК
	Matrix GS_ST;
	GS_ST.GSSK_STSK(B42, L42, PSI0, B_, L_);
// Матрица перехода от стартовой СК к географической СК
	//Matrix ST_GS = ~GS_ST;
//Матрица перехода от стартовой СК к связанной СК
	Matrix ST_SV;
	ST_SV.GetMatrix(ro_rg, lamda_rg, mu_rg, nu_rg);
//Матрица перехода от связанной СК к стартовой СК
	Matrix SV_ST = ST_SV;
	SV_ST = ~SV_ST;
//Матрица перехода от географического трехгранника к свзанному
	Matrix GS_SV;
	GS_SV.GSSK_SVSK(YAW, PITCH, ROLL);
//Матрица перехода от связанного трехгранника к географическому
	Matrix SV_GS = GS_SV;
	~SV_GS;
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
	Vector RateGS_GS(V_GS.Z / (H_ + OutStruct.R_lambda), V_GS.Z * tan(B_) / (H_ + OutStruct.R_lambda), -V_GS.X / (H_ + OutStruct.R_phi));
//Угловая скорость вращения осей связанной СК относительно осей географической СК в проекциях на оси связанной СК	
	Vector RateGS_SV = GS_SV * RateGS_GS;
//Угловая скорость вращения ЛА в проекциях на оси связанной СК
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
//Динамические уравнения Эйлера
	pr1Pitch = omegaY * sin(ROLL) + omegaZ * cos(ROLL);
	pr1Yaw = (omegaY * cos(ROLL) - omegaZ * sin(ROLL)) / cos(PITCH);
	pr1Roll = omegaX - tan(PITCH) * (omegaY * cos(ROLL) - omegaZ * sin(ROLL));
/*-------Расчет частных производных внешних сил и моментов, действующих на ЛА (кроме силы тяжести), в проекциях на СвСК----*/
	double Y_alfa = Missile1.Cy(mahh, ALFA)  * q * SquareMiddle;
	double Y_delta =Missile1.Cy_delta(mahh, ALFA)* q* SquareMiddle;
	double Z_betta = Missile1.Cz(mahh, BETTA)  * q * SquareMiddle;
	double Z_delta = Missile1.Cz_delta(mahh, BETTA) * q * SquareMiddle;
	double Mx_omgX = (Missile1.Mx_wx(mahh, AlfaSpace) * Rocket.LenghtLA / velocity) * q * SquareMiddle * Rocket.LenghtLA;
	double Mx_delta = Missile1.Mx_delta(q, SquareMiddle, Rocket.LenghtLA) * q * SquareMiddle * Rocket.LenghtLA;
	double My_omgY = (Missile1.My_wy(mahh, BETTA) * Rocket.LenghtLA / velocity)* q * SquareMiddle * Rocket.LenghtLA;
	double My_betta = Missile1.My_betta(mahh, BETTA) * q * SquareMiddle * Rocket.LenghtLA;
	double My_delta = Missile1.My_delta(mahh, BETTA) * q * SquareMiddle * Rocket.LenghtLA;
	double Mz_omgZ = Missile1.Mz_wz(mahh, ALFA) * Rocket.LenghtLA / velocity * q * SquareMiddle * Rocket.LenghtLA;
	double Mz_alfa = Missile1.Mz_alfa(mahh, ALFA) * q * SquareMiddle * Rocket.LenghtLA;
	double Mz_delta = Missile1.Mz_delta(mahh, ALFA) * q * SquareMiddle * Rocket.LenghtLA;
/*-------Расчет внешних сил и моментов, действующих на ЛА (кроме силы тяжести), в проекциях на СвСК----*/
	Vector F, M, Gg;
	
	F.X = -Missile1.Cx(mahh, AlfaSpace) * q * SquareMiddle;
	F.Y = Y_alfa * ALFA + Y_delta * DeltaPitch;
	F.Z = Z_betta * BETTA + Z_delta * DeltaYaw;

	M.X = Mx_omgX * omegaX + Mx_delta * DeltaRoll;
	M.Y = My_omgY * omegaY+ My_betta * BETTA + My_delta * DeltaYaw;
	M.Z = Mz_omgZ * omegaZ + Mz_alfa * ALFA + Mz_delta * DeltaPitch;
	
//Расчет динамических коэффициентов	
	double a11 = - Mz_omgZ / Rocket.Iz0;
	double a12 = - Mz_alfa / Rocket.Iz0;
	double a13 = - Mz_delta/ Rocket.Iz0;
	double a42 = Y_alfa/ (Rocket.mass * velocity);
	double a43 = Y_delta/ (Rocket.mass * velocity);
	double b11 = - My_omgY / Rocket.Iy0;
	double b12 = - My_betta/ Rocket.Iy0;
	double b13 = - My_delta/ Rocket.Iy0;
	double b42 = - Z_betta/ (Rocket.mass * velocity);
	double b43 = - Z_delta/ (Rocket.mass * velocity);//??
	double c11 = - Mx_omgX /  Rocket.Ix0;
	double c13 = - Mx_delta/ Rocket.Ix0;
//Расчет параметров СС
	double K_p = (a12 * a43 - a13 * a42) / (a12 + a11 * a42);
	double T1_p = -a13 / (a13 * a42 - a12 * a43);
	double T_p = 1 / sqrt(a12 + a11 +a42);
	double E_p = (a11 + a42) / (2 * sqrt(a12 + a11 * a42));
	double K_y = (b12 * b43 - b13 * b42) / (b12 + b11 * b42);
	double T1_y = -b13 / (b13 * b42 - b12 * b43);
	double T_y = 1 / sqrt(b12 + b11 + b42);
	double E_y = (b11 + b42) / (2 * sqrt(b12 + b11 * b42));
	double K_r = -c13 / c11;
	double T_r = 1 / c11;
	double Kcc_p = 0.95;
	double Kcc_y = 0.95;
	double Tcc_r = 0.01;
	double Ecc_p = 0.35;
	double Ecc_y = 0.35;
	double Ecc_r = 0.35;
	double K2_p = -2 * T_p * (E_p * T1_p - Ecc_p * Ecc_p * T_p - sqrt(pow(Ecc_p, 4) * _SQR(T_p) - 2 * E_p * _SQR(Ecc_p) * T1_p * T_p + _SQR(T1_p * Ecc_p)))/(K_p * T1_p * T1_p);
	double K2_y = -2 * T_y *(E_p * T1_y - Ecc_y * Ecc_y * T_y - sqrt(pow(Ecc_y, 4) * _SQR(T_y) - 2 * E_y * _SQR(Ecc_y) * T1_y * T_y + _SQR(T1_y * Ecc_y))) / (K_y * T1_y * T1_y);
	double K1_r = (2 * Ecc_r * T_r - Tcc_r) / (K_r * Tcc_r);
	double K2_r = T_r /( K_r * Tcc_r * Tcc_r);
//Расчет углов наклона рулей 
	DeltaPitch = -K2_p * pr1Pitch;
	DeltaYaw = -K2_y * pr1Yaw;
	DeltaRoll = -K2_r * ROLL - K1_r * pr1Roll;
//Проекция ускорения силы тяжести на оси географической СК
	Vector G_GS;
	G_GS = OutStruct.GField;
//Проекция ускорения силы тяжести на оси стартовой СК
	Vector G_ST = GS_ST * G_GS;
//Проекция силы тяжести на оси стартовой СК
	Gg = Rocket.mass * G_ST;
//Кажущееся ускорение в проекциях на оси связанной СК	
	Vector fg;
	fg = SV_ST * F;
	fg = fg + Gg;
//Ускорение Кориолиса в проекциях на оси стартовой СК
	Vector Accl_Kor = 2.0 * (RateEarth_ST % V_ST);
//Абсолютное ускорение в проекциях на оси стартовой СК
	Vector F_ST = fg * pow(Rocket.mass, -1.) - Accl_Kor;
	//Vector F_ST = fg * pow(Rocket.mass, -1.);
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

	pr1_RateAbs_SV = J_SV_inverse * (M - (pr1_J_SV * RateAbs_SV + RateAbs_SV % (J_SV * RateAbs_SV)));
	
//Уравнение Пуассона
	Matrix RateSV_I;
	RateSV_I.TestGetMatrix(0.0, RateAbs_SV.Z, RateAbs_SV.Y,
		RateAbs_SV.Z, 0.0, -RateAbs_SV.X,
		-RateAbs_SV.Y, RateAbs_SV.X, 0.0);
	Matrix RateGS_I;
	RateGS_I.TestGetMatrix(0.0, -RateEarth_GS.Z - RateGS_GS.Z, RateEarth_GS.Y + RateGS_GS.Y,
		RateEarth_GS.Z + RateGS_GS.Z, 0.0, -RateEarth_GS.X - RateGS_GS.X,
		-RateEarth_GS.Y - RateGS_GS.Y, RateEarth_GS.X + RateGS_GS.X, 0.0);
	Matrix Res1 = SV_GS * RateGS_I;
	Matrix Res2 = RateGS_I * SV_GS;
	Matrix pr1Matrix_SV_GS = Res1 - Res2;

//Уравнение Пуассона для СтСК
	Matrix RateST_I;
	RateST_I.TestGetMatrix(0.0, -RateEarth_ST.Z, RateEarth_ST.Y,
		RateEarth_ST.Z, 0.0, -RateEarth_ST.X,
		-RateEarth_ST.Y, RateEarth_ST.X, 0.0);
	 Res1 = SV_ST * RateGS_I;
	 Res2 = RateST_I * SV_ST;
	Matrix pr1Matrix_SV_ST = Res1 - Res2;

	pr1Speed_X = F_ST.X;
	pr1Speed_Y = F_ST.Y;
	pr1Speed_Z = F_ST.Z;
	pr1Point_X = Speed_X;
	pr1Point_Y = Speed_Y;
	pr1Point_Z = Speed_Z;
	pr1omegaX = pr1_RateAbs_SV.X;
	pr1omegaY = pr1_RateAbs_SV.Y;
	pr1omegaZ = pr1_RateAbs_SV.Z;
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

	if (Point_Y < 0)
		
	  do{
			Step = Step * 0.5;
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

	  } while (Point_Y < 0);
		  

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



