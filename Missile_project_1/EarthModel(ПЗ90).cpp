/*#include <math.h>
#include "EarthModel(ПЗ90).h"
#include "macros.h"

static	double	A_GD, B_GD;					//Большая полуось земного эллипсоида
static	double	E_GD,E2_GD;				//Первый  и второй эксцентриситеты земного эллипсоида
static	double	EARTH_RATE = 7.292115E-5;		//Угловая скорость вращения Земли 
static	double	G_EQUATOR = 9.7803284;			//Величина ускорения свободного падения на экваторе 

//Параметры модели ОЗЭ
static double R;                             //средний радиус Земли
static double fM;                            //гравитационная постоянная 
static double alpha;						 //сжатие ОЗЭ				
static double gamma_a;                     //ускорение свободного падения на экваторе
static double gamma_b;                     //ускорение свободного падения на полюсе

void EarthModelIni(double a_gd, double b_gd, double e_gd, double g_eq, double e_rate, double R_sr, double fM_, double alpha_, double gamma_a_, double gamma_b_)
{
	A_GD = a_gd;
	B_GD = b_gd;
	E_GD = e_gd;
	EARTH_RATE = e_rate;
	G_EQUATOR = g_eq;
	R = R_sr;
	fM = fM_;
	alpha = alpha_;
	gamma_a = gamma_a_;
	gamma_b = gamma_b_;
}

/*void matr_grin_sph(Vector SF, Vector GR) {

	GR.X = SF.X * cos(SF.Y) * sin(SF.Z);
	GR.Y = SF.X * sin(SF.Y) * sin(SF.Z);
	GR.Z = SF.X * cos(SF.Z);

}

void matr_grin_geod(Vector GDZ, Vector GR) {

	//массив входной: B, L, H

	double N = A_GD / sqrt(1 - E_GD * E_GD * pow(sin(GDZ.X), 2));

	GR.X = (N + GDZ.Z) * cos(GDZ.X) * cos(GDZ.Y);
	GR.Y = (N + GDZ.Z) * cos(GDZ.X) * sin(GDZ.Y);
	GR.Z = ((1 - E_GD * E_GD) * N + GDZ.Z) * sin(GDZ.X);

}

void matr_grin_geod(const double geod_coord[], double grin_coord[]) {

	//массив входной: B, L, H

	double N = A_GD / sqrt(1 - _SQR(E_GD) * pow(sin(geod_coord[0]), 2));

	grin_coord[0] = (N + geod_coord[2]) * cos(geod_coord[0]) * cos(geod_coord[1]);
	grin_coord[1] = (N + geod_coord[2]) * cos(geod_coord[0]) * sin(geod_coord[1]);
	grin_coord[2] = ((1 - _SQR(E_GD)) * N + geod_coord[2]) * sin(geod_coord[0]);

}


void matr_grin_sph(const double sph_coord[], double grin_coord[]) {

	//ro fi lambda
	grin_coord[0] = sph_coord[0] * cos(sph_coord[1]) * sin(sph_coord[2]);
	grin_coord[1] = sph_coord[0] * sin(sph_coord[1]) * sin(sph_coord[2]);
	grin_coord[2] = sph_coord[0] * cos(sph_coord[2]);

}


void matr_sph_grin(const double grin_coord[], double sph_coord[]) {

	sph_coord[0] = sqrt(pow(grin_coord[0], 2) + pow(grin_coord[1], 2) + pow(grin_coord[2], 2));
	sph_coord[1] = atan(grin_coord[1] / grin_coord[0]);
	sph_coord[2] = atan((sqrt(pow(grin_coord[0], 2) + pow(grin_coord[1], 2))) / (grin_coord[2]));
}

/*void matr_sph_grin(Vector GR, Vector SF) {

	SF.X = sqrt(pow(GR.X, 2) + pow(GR.Y, 2) + pow(GR.Z, 2));
	SF.Y = atan(GR.Y / GR.X);
	SF.Z = atan((sqrt(pow(GR.X, 2) + pow(GR.Y, 2))) / (GR.Z));

}


double g_norm(double B, double H) {

	double beta = (gamma_b - gamma_a) / (gamma_a);
	double beta1 = 0.25 * alpha * beta + 0.125 * pow(alpha, 2);
	double gamma_B_H0 = gamma_a * (1 + beta * pow(sin(B), 2) - beta1 * pow(sin(2 * B), 2));

	double m = pow(EARTH_RATE, 2) * A_GD / gamma_a;
	double f1 = (-2 * (1 + alpha + m)) / (A_GD);
	double f2 = 4 * alpha / A_GD;
	double f3 = 3 / pow(A_GD, 2);

	double delt_gamma_a = -0.87 * exp(-0.116 * pow(H, 1.047));

	double gamma_B_H = gamma_B_H0 * (1 + f1 * H + f2 * H * pow(sin(B), 2) + f3 * pow(H, 2)) + delt_gamma_a; //размерность в мГал

	double g = gamma_B_H / 100000;//перевод размерности в м/c^2
	return g;

}
double delt_g_anom(double B, double L, double H) {

	double matr_dot_mass[60][4] =
	{
		//eps * 10^10, X [км], Y[км], Z[км]
		{    -1917861.343    ,    -1597.53455    ,    3389.08854    ,    -1206.07844    }    ,
		{    -7649.811    ,    5243.88105    ,    2173.09105    ,    31.67769    }    ,
		{    -23204717.367    ,    -694.74764    ,    -2543.95209    ,    3010.82934    }    ,
		{    -8601525.203    ,    597.26083    ,    4124.29275    ,    3032.72188    }    ,
		{    -39262968.108    ,    1472.01984    ,    2496.1355    ,    2305.68074    }    ,
		{    -37656928.613    ,    -1395.72710    ,    522.57926    ,    1219.8023    }    ,
		{    55417.233    ,    -4729.43278    ,    241.0622    ,    -1386.75430    }    ,
		{    -1122388.574    ,    1613.74542    ,    -1389.46971    ,    -4742.37198    }    ,
		{    11553321.883    ,    -737.33302    ,    4223.69658    ,    -89.84823    }    ,
		{    -3212166.768    ,    -3498.93844    ,    2123.53095    ,    -2780.28445    }    ,
		{    -127570.128    ,    2443.27977    ,    -2075.43783    ,    -3115.11405    }    ,
		{    537718.408    ,    1399.38281    ,    -2714.24422    ,    -1999.26660    }    ,
		{    159431.526    ,    -1807.95598    ,    4038.16007    ,    -1984.43463    }    ,
		{    -9621707.493    ,    500.99426    ,    2417.51767    ,    2081.23547    }    ,
		{    14954047.019    ,    564.56372    ,    4169.6212    ,    3107.28151    }    ,
		{    -1454978.047    ,    -3606.22095    ,    2202.63232    ,    2555.27655    }    ,
		{    12173.463    ,    3494.94521    ,    -2225.41218    ,    3504.56634    }    ,
		{    6235318.028    ,    2819.10268    ,    1215.19399    ,    3031.84264    }    ,
		{    -4022993.036    ,    830.78873    ,    1187.36734    ,    2665.51499    }    ,
		{    -13233537.348    ,    -2407.50792    ,    -987.75326    ,    1168.97751    }    ,
		{    -837955.460    ,    2745.40696    ,    3312.38455    ,    561.61201    }    ,
		{    66944.487    ,    474.28039    ,    1205.35801    ,    -4546.03756    }    ,
		{    9533469.996    ,    -1675.28817    ,    -797.14875    ,    -3048.67319    }    ,
		{    13126.101    ,    1951.39223    ,    -5027.13641    ,    -1609.06073    }    ,
		{    133490487.677    ,    -1131.17750    ,    -625.93077    ,    614.81653    }    ,
		{    11202.018    ,    4552.11518    ,    1485.4265    ,    -1997.43144    }    ,
		{    -8306318.912    ,    32.85983    ,    1869.79667    ,    -2877.64543    }    ,
		{    1563860.673    ,    -3588.59466    ,    2192.00472    ,    2532.32304    }    ,
		{    -4110249.075    ,    -1688.88208    ,    -2509.59158    ,    997.3978    }    ,
		{    -708432.828    ,    142.71007    ,    -2365.17339    ,    3559.90481    }    ,
		{    115758.347    ,    -2583.62219    ,    3876.78111    ,    -1530.49081    }    ,
		{    47117.768    ,    -4410.11412    ,    -1491.17613    ,    1951.16572    }    ,
		{    3060709.563    ,    -3512.99151    ,    2136.54197    ,    -2789.00902    }    ,
		{    125694.63    ,    2887.19498    ,    -1520.35002    ,    -2929.00588    }    ,
		{    -19909835.138    ,    -1221.99104    ,    -1430.21334    ,    2886.63494    }    ,
		{    -25600026.158    ,    -795.89525    ,    4094.06397    ,    -148.43103    }    ,
		{    8998492.144    ,    -859.93100    ,    -2675.58249    ,    3049.08749    }    ,
		{    6460340.213    ,    134.47321    ,    -3777.76920    ,    -737.29613    }    ,
		{    22635646.931    ,    2030.1321    ,    2786.99103    ,    1094.12609    }    ,
		{    14963823.854    ,    -894.96828    ,    3946.71389    ,    -241.43049    }    ,
		{    104200607.824    ,    1424.84214    ,    2420.76271    ,    1954.22532    }    ,
		{    -29025.051    ,    3640.95418    ,    2890.38428    ,    -1299.38754    }    ,
		{    -25923827.511    ,    -1666.65137    ,    -1236.38163    ,    -177.16683    }    ,
		{    4710.121    ,    -2926.62377    ,    1305.39588    ,    4683.40996    }    ,
		{    -72296468.140    ,    1704.3058    ,    2570.01705    ,    1437.9415    }    ,
		{    11488188.566    ,    -1329.71266    ,    -1480.67585    ,    3122.43463    }    ,
		{    30871.376    ,    2993.13075    ,    -3221.87772    ,    2310.62018    }    ,
		{    -6821341.595    ,    539.11732    ,    4204.18554    ,    3168.98432    }    ,
		{    4053215.384    ,    520.10547    ,    -3574.66707    ,    -271.75382    }    ,
		{    2364380.549    ,    -1696.20051    ,    1611.26783    ,    2166.50206    }    ,
		{    -5107799.845    ,    -2806.42940    ,    90.032    ,    -437.20618    }    ,
		{    1152914.935    ,    1613.2997    ,    -1398.25374    ,    -4728.63714    }    ,
		{    14079228.613    ,    1331.61419    ,    -1428.22398    ,    234.32978    }    ,
		{    24148.39    ,    2240.6393    ,    2260.34015    ,    -3797.81205    }    ,
		{    14981240.134    ,    -523.57471    ,    -2437.64320    ,    3040.87357    }    ,
		{    -7746811.230    ,    2751.33446    ,    1231.49383    ,    2975.90618    }    ,
		{    -52402113.809    ,    490.2787    ,    -1222.02040    ,    299.60165    }    ,
		{    -9702053.744    ,    214.47365    ,    -3717.91494    ,    -621.72796    }    ,
		{    6998066.913    ,    43.26692    ,    1953.11104    ,    -2943.24995    }    ,
		{    -11022424.436    ,    -1649.80829    ,    -755.81754    ,    -2992.38699    }

	};

	double matr_dot_mas_sph[60][4] = {};
	for (int i = 0; i <= 59; i++) {
		matr_dot_mas_sph[i][0] = matr_dot_mass[i][0] / pow(10, 10);

		double dot_mas_grin_coord_temp[3] = {};
		double dot_mas_sph_coord_temp[3] = {};

		dot_mas_grin_coord_temp[0] = matr_dot_mass[i][1] * 1000;
		dot_mas_grin_coord_temp[1] = matr_dot_mass[i][2] * 1000;
		dot_mas_grin_coord_temp[2] = matr_dot_mass[i][3] * 1000;

		matr_sph_grin(dot_mas_grin_coord_temp, dot_mas_sph_coord_temp);

		matr_dot_mas_sph[i][1] = dot_mas_sph_coord_temp[0];
		matr_dot_mas_sph[i][2] = dot_mas_sph_coord_temp[1];
		matr_dot_mas_sph[i][3] = dot_mas_sph_coord_temp[2];

	}

	double geod_coord_la[3] = {};
	geod_coord_la[0] = B;
	geod_coord_la[1] = L;
	geod_coord_la[2] = H;

	//координаты ЛА в геоцентрических прямоугольных координатах
	double grin_coord_la[3] = {};
	matr_grin_geod(geod_coord_la, grin_coord_la);

	//координаты ЛА в геоцентрических сферических координатах
	double sph_coord_la[3] = {};
	matr_sph_grin(grin_coord_la, sph_coord_la);

	double delt_g = 0;

	for (int i = 0; i <= 59; i++) {

		double ro = sph_coord_la[0];
		double fi = sph_coord_la[1];
		double lambda = sph_coord_la[2];

		double eps_i = matr_dot_mas_sph[i][0];
		double ro_i = matr_dot_mas_sph[i][1];
		double fi_i = matr_dot_mas_sph[i][2];
		double lambda_i = matr_dot_mas_sph[i][3];

		double cos_psi = sin(fi_i) * sin(fi) + cos(fi_i) * cos(fi) * cos(lambda_i - lambda);
		double r_i = sqrt(pow(ro, 2) + pow(ro_i, 2) - 2 * ro * ro_i * cos_psi);

		delt_g = delt_g + eps_i * (pow(R, 2) - pow(ro_i, 2) - 3 * pow(r_i, 2)) / (2 * R * pow(r_i, 3));

	}

	delt_g = delt_g * fM;

	return delt_g;
}


/*double delt_g_anom(double B, double L, double H) {

	double matr_dot_mass[60][4] =
	{
		//eps * 10^10, X [км], Y[км], Z[км]
		{    -1917861.343    ,    -1597.53455    ,    3389.08854    ,    -1206.07844    }    ,
		{    -7649.811    ,    5243.88105    ,    2173.09105    ,    31.67769    }    ,
		{    -23204717.367    ,    -694.74764    ,    -2543.95209    ,    3010.82934    }    ,
		{    -8601525.203    ,    597.26083    ,    4124.29275    ,    3032.72188    }    ,
		{    -39262968.108    ,    1472.01984    ,    2496.1355    ,    2305.68074    }    ,
		{    -37656928.613    ,    -1395.72710    ,    522.57926    ,    1219.8023    }    ,
		{    55417.233    ,    -4729.43278    ,    241.0622    ,    -1386.75430    }    ,
		{    -1122388.574    ,    1613.74542    ,    -1389.46971    ,    -4742.37198    }    ,
		{    11553321.883    ,    -737.33302    ,    4223.69658    ,    -89.84823    }    ,
		{    -3212166.768    ,    -3498.93844    ,    2123.53095    ,    -2780.28445    }    ,
		{    -127570.128    ,    2443.27977    ,    -2075.43783    ,    -3115.11405    }    ,
		{    537718.408    ,    1399.38281    ,    -2714.24422    ,    -1999.26660    }    ,
		{    159431.526    ,    -1807.95598    ,    4038.16007    ,    -1984.43463    }    ,
		{    -9621707.493    ,    500.99426    ,    2417.51767    ,    2081.23547    }    ,
		{    14954047.019    ,    564.56372    ,    4169.6212    ,    3107.28151    }    ,
		{    -1454978.047    ,    -3606.22095    ,    2202.63232    ,    2555.27655    }    ,
		{    12173.463    ,    3494.94521    ,    -2225.41218    ,    3504.56634    }    ,
		{    6235318.028    ,    2819.10268    ,    1215.19399    ,    3031.84264    }    ,
		{    -4022993.036    ,    830.78873    ,    1187.36734    ,    2665.51499    }    ,
		{    -13233537.348    ,    -2407.50792    ,    -987.75326    ,    1168.97751    }    ,
		{    -837955.460    ,    2745.40696    ,    3312.38455    ,    561.61201    }    ,
		{    66944.487    ,    474.28039    ,    1205.35801    ,    -4546.03756    }    ,
		{    9533469.996    ,    -1675.28817    ,    -797.14875    ,    -3048.67319    }    ,
		{    13126.101    ,    1951.39223    ,    -5027.13641    ,    -1609.06073    }    ,
		{    133490487.677    ,    -1131.17750    ,    -625.93077    ,    614.81653    }    ,
		{    11202.018    ,    4552.11518    ,    1485.4265    ,    -1997.43144    }    ,
		{    -8306318.912    ,    32.85983    ,    1869.79667    ,    -2877.64543    }    ,
		{    1563860.673    ,    -3588.59466    ,    2192.00472    ,    2532.32304    }    ,
		{    -4110249.075    ,    -1688.88208    ,    -2509.59158    ,    997.3978    }    ,
		{    -708432.828    ,    142.71007    ,    -2365.17339    ,    3559.90481    }    ,
		{    115758.347    ,    -2583.62219    ,    3876.78111    ,    -1530.49081    }    ,
		{    47117.768    ,    -4410.11412    ,    -1491.17613    ,    1951.16572    }    ,
		{    3060709.563    ,    -3512.99151    ,    2136.54197    ,    -2789.00902    }    ,
		{    125694.63    ,    2887.19498    ,    -1520.35002    ,    -2929.00588    }    ,
		{    -19909835.138    ,    -1221.99104    ,    -1430.21334    ,    2886.63494    }    ,
		{    -25600026.158    ,    -795.89525    ,    4094.06397    ,    -148.43103    }    ,
		{    8998492.144    ,    -859.93100    ,    -2675.58249    ,    3049.08749    }    ,
		{    6460340.213    ,    134.47321    ,    -3777.76920    ,    -737.29613    }    ,
		{    22635646.931    ,    2030.1321    ,    2786.99103    ,    1094.12609    }    ,
		{    14963823.854    ,    -894.96828    ,    3946.71389    ,    -241.43049    }    ,
		{    104200607.824    ,    1424.84214    ,    2420.76271    ,    1954.22532    }    ,
		{    -29025.051    ,    3640.95418    ,    2890.38428    ,    -1299.38754    }    ,
		{    -25923827.511    ,    -1666.65137    ,    -1236.38163    ,    -177.16683    }    ,
		{    4710.121    ,    -2926.62377    ,    1305.39588    ,    4683.40996    }    ,
		{    -72296468.140    ,    1704.3058    ,    2570.01705    ,    1437.9415    }    ,
		{    11488188.566    ,    -1329.71266    ,    -1480.67585    ,    3122.43463    }    ,
		{    30871.376    ,    2993.13075    ,    -3221.87772    ,    2310.62018    }    ,
		{    -6821341.595    ,    539.11732    ,    4204.18554    ,    3168.98432    }    ,
		{    4053215.384    ,    520.10547    ,    -3574.66707    ,    -271.75382    }    ,
		{    2364380.549    ,    -1696.20051    ,    1611.26783    ,    2166.50206    }    ,
		{    -5107799.845    ,    -2806.42940    ,    90.032    ,    -437.20618    }    ,
		{    1152914.935    ,    1613.2997    ,    -1398.25374    ,    -4728.63714    }    ,
		{    14079228.613    ,    1331.61419    ,    -1428.22398    ,    234.32978    }    ,
		{    24148.39    ,    2240.6393    ,    2260.34015    ,    -3797.81205    }    ,
		{    14981240.134    ,    -523.57471    ,    -2437.64320    ,    3040.87357    }    ,
		{    -7746811.230    ,    2751.33446    ,    1231.49383    ,    2975.90618    }    ,
		{    -52402113.809    ,    490.2787    ,    -1222.02040    ,    299.60165    }    ,
		{    -9702053.744    ,    214.47365    ,    -3717.91494    ,    -621.72796    }    ,
		{    6998066.913    ,    43.26692    ,    1953.11104    ,    -2943.24995    }    ,
		{    -11022424.436    ,    -1649.80829    ,    -755.81754    ,    -2992.38699    }

	};

	double matr_dot_mas_sph[60][4] = {};
	for (int i = 0; i <= 59; i++) {
		matr_dot_mas_sph[i][0] = matr_dot_mass[i][0] / pow(10, 10);

		Vector tempGR, tempSF;

		tempGR.X = matr_dot_mass[i][1] * 1000;
		tempGR.Y = matr_dot_mass[i][2] * 1000;
		tempGR.Z = matr_dot_mass[i][3] * 1000;

		matr_sph_grin(tempGR, tempSF);

		matr_dot_mas_sph[i][1] = tempSF.X;
		matr_dot_mas_sph[i][2] = tempSF.Y;
		matr_dot_mas_sph[i][3] = tempSF.Z;
	}

	double geod_coord_la[3] = {};
	Vector GDZ(B, L, H), GR, G;
	geod_coord_la[0] = B;
	geod_coord_la[1] = L;
	geod_coord_la[2] = H;

	//координаты ЛА в геоцентрических прямоугольных координатах
	matr_grin_geod(GDZ, GR);

	//координаты ЛА в геоцентрических сферических координатах
	Vector SF;
	matr_sph_grin(GR, SF);


	double delt_g = 0;

	for (int i = 0; i <= 59; i++) {
		double ro = SF.X;
		double fi = SF.Y;
		double lambda = SF.Z;

		double eps_i = matr_dot_mas_sph[i][0];
		double ro_i = matr_dot_mas_sph[i][1];
		double fi_i = matr_dot_mas_sph[i][2];
		double lambda_i = matr_dot_mas_sph[i][3];

		double cos_psi = sin(fi_i) * sin(fi) + cos(fi_i) * cos(fi) * cos(lambda_i - lambda);
		double r_i = sqrt(pow(ro, 2) + pow(ro_i, 2) - 2 * ro * ro_i * cos_psi);

		delt_g = delt_g + eps_i * (pow(R, 2) - pow(ro_i, 2) - 3 * pow(r_i, 2)) / (2 * R * pow(r_i, 3));
	}
	
	delt_g = delt_g * fM;
	return delt_g;
}

Earth_Struct	GetEarthParameters(double H, double B, double L)
{
	//Функция расчёта выходных параметров

	double R_gd;
	Earth_Struct OutStruct;
	double  W_gd, V_gd;				//Основные сфероидические функции
	double  R_phi, R_lambda;		//Главные радиусы кривизны поверхности эллипсоида вращения
	double  R_gd0;					//Вспомогательный радиус-вектор на поверхности Земли в ГцССК
	double  phi_gr0, phi_gr;		//Геоцентрическая широта
				
	W_gd = sqrt(1 - _SQR(E_GD * sin(B)));
	V_gd = sqrt(1 - _SQR(E2_GD * cos(B)));
	R_phi = A_GD * (1 - _SQR(E_GD)) / (_SQR(W_gd) * W_gd);
	R_lambda = A_GD / W_gd;
	phi_gr0 = atan(_SQR(B_GD) / _SQR(A_GD) * tan(B));
	R_gd0 = sqrt(_SQR(R_lambda * cos(B)) + _SQR(R_lambda * (1 - _SQR(E_GD)) * sin(B)));
	R_gd = sqrt(_SQR(R_gd0 * cos(phi_gr0) + H * cos(B)) + _SQR(R_gd0 * sin(phi_gr0) + H * sin(B)));
	phi_gr = asin((R_gd0 * sin(phi_gr0) + H * sin(B)) / R_gd);
	
	Vector F_G;//Проекции вектора напряженности поля тяготения на норм земную СК

	double g0 = g_norm(B, H) + delt_g_anom(B, L, H);
	double fi_g = atan((1 - _SQR(E_GD)) * tan(B));
	double R1 = sqrt(_SQR(A_GD) * (1 - _SQR(E_GD)) / (_SQR(sin(fi_g)) + (1 - _SQR(E_GD)) * _SQR(cos(fi_g))));
	double R2 = sqrt(_SQR(R1) + _SQR(H) + 2 * R1 * H * cos(B - fi_g));
	double D = (H * cos(B) + R1 * cos(fi_g) * (1 - _SQR(R1 / R2)));
	F_G.X = _SQR(EARTH_RATE) * D * sin(B);
	F_G.Y = - (g0 * _SQR(R1 / R2) - _SQR(EARTH_RATE) * D * cos(B));
	F_G.Z = 0;

	OutStruct.A_gd = A_GD;
	OutStruct.G_eq = gamma_a;
	OutStruct.Rate = EARTH_RATE;
	OutStruct.E_gd = E_GD;
	OutStruct.R_phi = R_phi;
	OutStruct.R_lambda = R_lambda;
	OutStruct.R_gd = R_gd;
	OutStruct.phi_gr = phi_gr;
	OutStruct.AcclCp_GS = 0;
	OutStruct.GField = F_G;

	return OutStruct;
}

void Gss_Kr_to_SK42(double x, double y, double& B42_, double& L42_)
{
	double beta = x / 6367558.4968;
	double B0 = beta + sin(2 * beta) * (0.00252588685 - 0.00001491860 * _SQR(sin(beta)) + 0.00000011904 * _SQR(_SQR(sin(beta))));
	int n = y * 1E-6;
	double z0 = (y - (10.0 * n + 5.0) * 1E5) / (6378245 * cos(B0));
	double k = pow(sin(B0), 6);
	double b1 = (0.01672 - 0.0063 * _SQR(sin(B0)) + 0.01188 * _SQR(_SQR(sin(B0))) - 0.00328 * pow(sin(B0), 6));
	double b2 = 0.042858 - 0.025318 * _SQR(sin(B0)) + 0.014346 * _SQR(_SQR(sin(B0))) - 0.001264 * pow(sin(B0), 6);
	double b3 = 0.10500614 - 0.04559916 * _SQR(sin(B0)) + 0.00228901 * _SQR(_SQR(sin(B0))) - 0.00002987 * pow(sin(B0), 6);
	double b4 = 0.251684631 - 0.003369263 * _SQR(sin(B0)) + 0.000011276 * _SQR(_SQR(sin(B0)));
	double deltaB = -_SQR(z0) * sin(2 * B0) * (b4 - _SQR(z0) * (b3 - _SQR(z0) * (b2 - _SQR(z0) * b1)));
	double l1 = 0.0038 + 0.0524 * _SQR(sin(B0)) + 0.0482 * _SQR(_SQR(sin(B0))) + 0.0032 * pow(sin(B0), 6);
	double l2 = 0.01225 + 0.09477 * _SQR(sin(B0)) + 0.03282 * _SQR(_SQR(sin(B0))) - 0.00034 * pow(sin(B0), 6);
	double l3 = 0.0420025 + 0.1487407 * _SQR(sin(B0)) + 0.005942 * _SQR(_SQR(sin(B0))) - 0.000015 * pow(sin(B0), 6);
	double l4 = 0.16778975 + 0.16273586 * _SQR(sin(B0)) - 0.0005249 * _SQR(_SQR(sin(B0))) - 0.00000846 * pow(sin(B0), 6);
	double l5 = 1 - 0.0033467108 * _SQR(sin(B0)) - 0.0000056002 * _SQR(_SQR(sin(B0))) - 0.0000000187 * pow(sin(B0), 6);
	double l = z0 * (l5 - _SQR(sin(z0)) * (l4 - _SQR(sin(z0)) * (l3 - _SQR(sin(z0)) * (l2 - _SQR(sin(z0)) * l1))));
	B42_ = B0 + deltaB;
	L42_ = 6 * (n - 0.5) / 57.29577951 + l;
}*/