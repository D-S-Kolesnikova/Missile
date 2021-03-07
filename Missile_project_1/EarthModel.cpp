#define _USE_MATH_DEFINES

#include <math.h>

#include "EarthModel.h"
#include "macros.h"

//Параметры модели Земли, необходимые для расчёта (числовые значения в скобках указаны для эллипсоида Красовского)
static	double	A_GD, B_GD;				//Большая (6378245.0 м) и малая (6356863.0 м) полуоси земного эллипсоида
static	double	E_GD, E2_GD;			//Первый (0.08181337) и второй эксцентриситеты земного эллипсоида
static	double	EARTH_RATE;				//Угловая скорость вращения Земли (7.292116E-5 рад/сек)
static	double	G_EQUATOR;				//Величина ускорения свободного падения на экваторе (9.78049)
static	double	Q_GR;					//Отношение центробежной силы, возникающей вследствие вращения Земли, к силе тяжести на экваторе
static	double	K1, K2, K3, K4, K5, K6;	//Вспомогательные коэффициенты



void EarthModelIni(double a_gd, double e_gd, double g_eq, double e_rate)
{
	A_GD = a_gd;
	E_GD = e_gd;
	EARTH_RATE = e_rate;
	G_EQUATOR = g_eq;
	B_GD = A_GD * sqrt(1 - _SQR(E_GD));
	E2_GD = sqrt(_SQR(A_GD) - _SQR(B_GD)) / B_GD;
	Q_GR = _SQR(EARTH_RATE) * A_GD / G_EQUATOR;

	K1 = G_EQUATOR * 0.5 * (Q_GR - _SQR(E_GD)) * (1.0 + _SQR(E_GD) * (7.0 * _SQR(E_GD) - 30.0 * Q_GR) / (14.0 * (Q_GR - _SQR(E_GD))));
	K2 = K1 * _SQR(E_GD) * (30.0 * Q_GR - 21.0 * _SQR(E_GD)) / (14.0 * (Q_GR - _SQR(E_GD)));
	K3 = K1 * _SQR(E_GD) * (7.0 * _SQR(E_GD) - 10.0 * Q_GR) / (2.0 * (Q_GR - _SQR(E_GD)));
	K4 = -G_EQUATOR * (1.0 - 0.5 * _SQR(E_GD) * (1.0 + 0.25 * _SQR(E_GD)) + 1.5 * Q_GR * (1.0 - 5.0 / 14.0 * _SQR(E_GD)));
	K5 = -G_EQUATOR * (_SQR(E_GD) * (1.0 - 0.5 * _SQR(E_GD)) - Q_GR * (1.0 - 15.0 / 7.0 * _SQR(E_GD))) / 2.0;
	K6 = G_EQUATOR * (_SQR(E_GD) - Q_GR - _SQR(E_GD) * (0.5 * _SQR(E_GD) - 15.0 / 7.0 * Q_GR)) * 1.5;

}
Earth_Struct	GetEarthParameters(double H, double B)
{
	//Функция расчёта выходных параметров
	//H - высота над поверхностью референц-эллипсоида Красовского
	//B - геодезическая широта
	Earth_Struct OutStruct;
	double  W_gd, V_gd;				//Основные сфероидические функции
	double  R_phi, R_lambda;		//Главные радиусы кривизны поверхности эллипсоида вращения
	double  R_gd0;					//Вспомогательный радиус-вектор на поверхности Земли в ГцССК
	double  phi_gr0, phi_gr;		//Геоцентрическая широта
	Vector F1_GCS;					//Проекции вектора напряженности поля тяготения на ГцССК
	Vector F1_GS;					//Проекции вектора напряженности поля тяготения на ГССК
	Vector AcclCp_GS;				//Проекции вектора центростремительного ускорения на ГССК
	Vector F2_GS;					//Проекции вектора напряженности поля силы тяжести на ГССК

	double R_gd;

	W_gd = sqrt(1 - _SQR(E_GD * sin(B)));
	V_gd = sqrt(1 - _SQR(E2_GD * cos(B)));
	R_phi = A_GD * (1 - _SQR(E_GD)) / (_SQR(W_gd) * W_gd);
	R_lambda = A_GD / W_gd;
	R_gd0 = sqrt(_SQR(R_lambda * cos(B)) + _SQR(R_lambda * (1 - _SQR(E_GD)) * sin(B)));
	phi_gr0 = atan(_SQR(B_GD) / _SQR(A_GD) * tan(B));
	R_gd = sqrt(_SQR(R_gd0 * cos(phi_gr0) + H * cos(B)) + _SQR(R_gd0 * sin(phi_gr0) + H * sin(B)));
	phi_gr = asin((R_gd0 * sin(phi_gr0) + H * sin(B)) / R_gd);

	F1_GCS.X = _SQR(A_GD / R_gd) * _SQR(A_GD / R_gd) * sin(2 * phi_gr) * (K1 + _SQR(A_GD / R_gd) * (K2 + K3 * _SQR(sin(phi_gr))));
	F1_GCS.Y = _SQR(A_GD / R_gd) * (K4 + _SQR(A_GD / R_gd) * (K5 + K6 * _SQR(sin(phi_gr))));
	F1_GCS.Z = 0.;

	F1_GS.X = F1_GCS.X * cos(B - phi_gr) - F1_GCS.Y * sin(B - phi_gr);
	F1_GS.Y = F1_GCS.X * sin(B - phi_gr) + F1_GCS.Y * cos(B - phi_gr);
	F1_GS.Z = 0.;

	AcclCp_GS.X = _SQR(EARTH_RATE) * R_gd * (sin(B) * cos(B) * cos(B - phi_gr) + _SQR(sin(B)) * sin(B - phi_gr));
	AcclCp_GS.Y = _SQR(EARTH_RATE) * R_gd * (-_SQR(cos(B)) * cos(B - phi_gr) - sin(B) * cos(B) * sin(B - phi_gr));
	AcclCp_GS.Z = 0.;

	F2_GS = F1_GS - AcclCp_GS;

	OutStruct.A_gd = A_GD;
	OutStruct.G_eq = G_EQUATOR;
	OutStruct.Rate = EARTH_RATE;
	OutStruct.E_gd = E_GD;
	OutStruct.R_phi = R_phi;
	OutStruct.R_lambda = R_lambda;
	OutStruct.R_gd = R_gd;
	OutStruct.phi_gr = phi_gr;
	OutStruct.AcclCp_GS = F1_GS;	 //Удалить, когда закончу с проверкой модели Земли
	OutStruct.GField = F2_GS;

	return OutStruct;
}

void Gss_Kr_to_SK42(double x, double y, double& B42_, double& L42_)
{
	double beta = x / 6367558.4968;
	double B0 = beta + sin(2 * beta) * (0.00252588685 - 0.00001491860 * _SQR(sin(beta)) + 0.00000011904 * _SQR(_SQR(sin(beta))));
	int n = y * 1E-6;
	double z0 = (y - (10.0 * n + 5.0) * 1E5) / (6378245 * cos(B0));
	double k = pow(sin(B0),6);
	double b1 = (0.01672 - 0.0063 * _SQR(sin(B0)) + 0.01188 * _SQR(_SQR(sin(B0))) - 0.00328 * pow(sin(B0), 6));
	double b2 = 0.042858 - 0.025318 * _SQR(sin(B0)) + 0.014346 * _SQR(_SQR(sin(B0))) - 0.001264 * pow(sin(B0), 6);
	double b3 = 0.10500614 - 0.04559916 * _SQR(sin(B0)) + 0.00228901 * _SQR(_SQR(sin(B0))) - 0.00002987 * pow(sin(B0), 6);
	double b4 = 0.251684631 - 0.003369263 * _SQR(sin(B0)) + 0.000011276 * _SQR(_SQR(sin(B0)));
	double deltaB = -_SQR(z0) * sin(2 * B0) * (b4 - _SQR(z0) * (b3 - _SQR(z0) * (b2 - _SQR(z0)* b1)));
	double l1 = 0.0038 + 0.0524 * _SQR(sin(B0)) + 0.0482 * _SQR(_SQR(sin(B0))) + 0.0032 * pow(sin(B0), 6);
	double l2 = 0.01225 + 0.09477 * _SQR(sin(B0)) + 0.03282 * _SQR(_SQR(sin(B0))) - 0.00034 * pow(sin(B0), 6);
	double l3 = 0.0420025 + 0.1487407 * _SQR(sin(B0)) + 0.005942 * _SQR(_SQR(sin(B0))) - 0.000015 * pow(sin(B0), 6);
	double l4 = 0.16778975 + 0.16273586 * _SQR(sin(B0)) - 0.0005249 * _SQR(_SQR(sin(B0))) - 0.00000846 * pow(sin(B0), 6);
	double l5 = 1 - 0.0033467108 * _SQR(sin(B0)) - 0.0000056002 * _SQR(_SQR(sin(B0))) - 0.0000000187 * pow(sin(B0), 6);
	double l = z0 * (l5 - _SQR(sin(z0)) * (l4 - _SQR(sin(z0)) * (l3 - _SQR(sin(z0)) * (l2 - _SQR(sin(z0)) * l1))));
	B42_ = B0 + deltaB;
	L42_ = 6 * (n - 0.5) / 57.29577951 + l;
}

/*Earth_Struct	GetEarthParameters(double H, double B) {
	//Функция расчёта выходных параметров
	//H - высота над поверхностью референц-эллипсоида Красовского
	//B - геодезическая широта
	Earth_Struct OutStruct;
	double  e_gd, e2_gd;  			//Первый и второй эксцентриситеты земного эллипсоида
	double  W_gd, V_gd;				//Основные сфероидические функции
	double  R_phi, R_lambda;		//Главные радиусы кривизны поверхности эллипсоида вращения
	double  R_gd0;								//Вспомогательный радиус-вектор на поверхности Земли в ГцССК
	double  phi_gr0, phi_gr;					//Геоцентрическая широта
	double  q_gr;					//Отношение центробежной силы, возникающей вследствие вращения Земли, к силе тяжести на экваторе
	Vector3 F1_GCS;					//Проекции вектора напряженности поля тяготения на ГцССК
	Vector3 F1_GS;					//Проекции вектора напряженности поля тяготения на ГССК
	Vector3 AcclCp_GS;				//Проекции вектора центростремительного ускорения на ГССК
	Vector3 F2_GS;					//Проекции вектора напряженности поля силы тяжести на ГССК
	double K1, K2, K3, K4, K5, K6;	//Вспомогательные коэффициенты
	double R_gd;					//Радиус-вектор положения центра масс ЛА в ГцССК
	double u_gd;
	double l_gr;
	double l1_gr;
	double Dmu_gr;
	double a1, a2, a3;
	double P0_gr;

	double R_gd, e_gd, e2_gd, u_gd, l_gr, l1_gr, Dmu_gr, P0_gr, C_gr, Nu_gr, q_gr;
	double a1, a2, a3;//Радиус-вектор положения центра масс ЛА в ГцССК
	e_gd = sqrt(_SQR(A_GD) - _SQR(B_GD)) / A_GD; //Переместить в раздел инициализации
	e2_gd = sqrt(_SQR(A_GD) - _SQR(B_GD)) / B_GD; //Переместить в раздел инициализации
	W_gd = sqrt(1 - _SQR(e_gd * sin(B)));
	V_gd = sqrt(1 + _SQR(e_gd * sin(B)));
	R_phi = A_GD * (1 - _SQR(e_gd)) / (_SQR(W_gd) * W_gd);
	R_lambda = A_GD / W_gd;
	//R_gd0 = sqrt(_SQR(R_lambda * cos(B)) + _SQR(R_lambda * (1 - _SQR(e_gd)) * sin(B)));
	//phi_gr0 = atan(_SQR(B_GD) / _SQR(A_GD) * tan(B));
	R_gd = sqrt(_SQR(R_lambda + H) * _SQR(cos(B)) + _SQR(R_lambda * (1 - _SQR(e_gd)) + H) * _SQR(sin(B)));
	u_gd = atan(B_GD / A_GD * tan(B));
	phi_gr = atan(B_GD / A_GD * tan(u_gd));
	l_gr = sqrt(_SQR(E_GD) / (1 - _SQR(E_GD)));
	a1 = 2 * (_SQR(A_GD) - _SQR(B_GD));
	a2 = _SQR(R_gd) - (_SQR(A_GD) - _SQR(B_GD));
	a3 = sqrt(pow(R_gd, 4) + _SQR(_SQR(A_GD) - _SQR(B_GD)) - 2 * R_gd * R_gd * (_SQR(A_GD) - _SQR(B_GD)) * cos(2 * phi_gr));
	l1_gr = sqrt((2 * (_SQR(A_GD) - _SQR(B_GD))) / (_SQR(R_gd) - (_SQR(A_GD) - _SQR(B_GD)) + sqrt(pow(R_gd, 4) + _SQR(_SQR(A_GD) - _SQR(B_GD)) - 2 * R_gd * R_gd * (_SQR(A_GD) - _SQR(B_GD)) * cos(2 * phi_gr))));
	Dmu_gr = _SQR(EARTH_RATE) * pow(l_gr, 3) / (2 * M_PI * ((3 + l_gr * l_gr) * atan(l_gr) - 3 * l_gr));
	P0_gr = 2 * M_PI * Dmu_gr * (1 + l_gr * l_gr) / pow(l_gr, 3) * (atan(l_gr) - l_gr / (1 + l_gr * l_gr));
	C_gr = (G_EQUATOR + _SQR(EARTH_RATE) * A_GD - P0_gr * A_GD) / (2 * M_PI * Dmu_gr * A_GD);
	Nu_gr = (_SQR(A_GD) - _SQR(B_GD)) / _SQR(l1_gr) - _SQR(B_GD);
	q_gr = _SQR(EARTH_RATE) * A_GD / G_EQUATOR;

	F_GCS.X = 2 * M_PI * Dmu_gr * sin(phi_gr) * cos(phi_gr) * (_SQR(A_GD) * B_GD * R_gd * pow((_SQR(A_GD) - _SQR(B_GD)), -1.5) *
		(3* atan(l1_gr) - l1_gr/( 1 + _SQR(l1_gr)) - 2 * l1_gr ) - C_gr * pow(A_GD,3) * B_GD * (A_GD/R_gd) * _SQR(e_gd) *
		sqrt(_SQR(B_GD) + Nu_gr) / (_SQR(_SQR(B_GD) + Nu_gr) * _SQR(cos(phi_gr)) + _SQR(_SQR(A_GD) + Nu_gr) * _SQR(sin(phi_gr))));
	F_GCS.Y = -2 * M_PI * Dmu_gr * (_SQR(A_GD) * B_GD * R_gd * pow((_SQR(A_GD) - _SQR(B_GD)), -1.5) * (2 * (l1_gr - atan(l1_gr)) + _SQR(cos(phi_gr)) *
		(3 * atan(l1_gr) - l1_gr / (1 + _SQR(l1_gr)) - 2 * l1_gr)) + C_gr * A_GD * B_GD * (A_GD / R_gd)) * sqrt(_SQR(B_GD) + Nu_gr) *
		((_SQR(B_GD) + Nu_gr) * _SQR(cos(phi_gr)) + (_SQR(A_GD) + Nu_gr) * _SQR(sin(phi_gr)))/ (_SQR(_SQR(B_GD) + Nu_gr) *_SQR(cos(phi_gr)) +
		_SQR(_SQR(A_GD) + Nu_gr) * _SQR(sin(phi_gr)));
	F_GCS.Z = 0;

	F1_GCS.X = _SQR(A_GD / R_gd) * _SQR(A_GD / R_gd) * sin(2 * phi_gr) * (K1 + _SQR(A_GD / R_gd) * (K2 + K3 * _SQR(sin(phi_gr))));
	F1_GCS.Y = _SQR(A_GD / R_gd) * (K4 + _SQR(A_GD / R_gd) * (K5 + K6 * _SQR(sin(phi_gr))));
	F1_GCS.Z = 0.;

	//F1_GCS.X = _SQR(A_GD / R_gd) * _SQR(A_GD / R_gd) * sin(2 * phi_gr) * (K1 + _SQR(A_GD / R_gd) * (K2 + K3 * _SQR(sin(phi_gr))));
	//F1_GCS.Y = _SQR(A_GD / R_gd) * (K4 + _SQR(A_GD / R_gd) * (K5 + K6 * _SQR(sin(phi_gr))));
	//F1_GCS.Z = 0.;

	//F1_GS.X = F1_GCS.X * cos(B - phi_gr) - F1_GCS.Y * sin(B - phi_gr);
	//F1_GS.Y = F1_GCS.X * sin(B - phi_gr) + F1_GCS.Y * cos(B - phi_gr);
	//F1_GS.Z = 0.;

	//AcclCp_GS.X = _SQR(EARTH_RATE) * R_gd * (sin(B) * cos(B) * cos(B - phi_gr) + _SQR(sin(B)) * sin(B - phi_gr));
	//AcclCp_GS.Y = _SQR(EARTH_RATE) * R_gd * (-_SQR(cos(B)) * cos(B - phi_gr) - sin(B) * cos(B) * sin(B - phi_gr));
	//AcclCp_GS.Z = 0.;

	//F2_GS = sub_vector3(F1_GS, AcclCp_GS);

	OutStruct.A_gd = A_GD;
	OutStruct.Rate = EARTH_RATE;
	OutStruct.E_gd = e_gd;
	OutStruct.R_phi = R_phi;
	OutStruct.R_lambda = R_lambda;
	OutStruct.R_gd = R_gd;
	OutStruct.phi_gr = phi_gr;
	//OutStruct.AcclCp_GS = F1_GS;	
	//OutStruct.GField = F2_GS;

	return OutStruct;
}

//double Fx = G_EQUATOR * 0.5 * (Q_GR - _SQR(E_GD)) * _SQR(A_GD / R_gd) * _SQR(A_GD / R_gd) * sin(2 * phi_gr)*
//		(1.0 + _SQR(E_GD) * (7.0 * _SQR(E_GD) - 30.0 * Q_GR) / (14.0 * (Q_GR - _SQR(E_GD))))+
//		(1.0 + ((30.0 * Q_GR - 21.0 * _SQR(E_GD)) / (14.0 * (Q_GR - _SQR(E_GD))) + _SQR(sin(phi_gr))*(7.0 * _SQR(E_GD) - 10.0 * Q_GR) / (2.0 * (Q_GR - _SQR(E_GD)))) * _SQR(E_GD) * _SQR(A_GD/R_gd));*/