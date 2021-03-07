#pragma once
typedef struct Initial_Conditions
{
	double X0;
	double Y0;
	double Z0;
	double X_1;
	double Y_1;
	double H_1;
	double V0;
	double Alfa0;
	double Betta0;
	double Pitch0;
	double Yaw0;
	double Roll0;
	double Wx0;
	double Wy0;
	double Wz0;
}Initial_Conditions;

typedef struct Object
{
	double mass;
	double Ix0;
	double Iy0;
	double Iz0;
	double Dmiddle;
	double LenghtLA;
	
} Object;



typedef struct _ModelParams_Struct
{
	double Step;    		//Шаг моделирования полёта
	double PrintStep;		//Шаг печати в файл
	
} ModelParams_Struct;

typedef	struct _Matrix_Struct
{
	Matrix GS_SV;		//Матрица направляющих косинусов между ГССК и СвСК
	Matrix GS_ST;		//Матрица направляющих косинусов между ГССК и СтСК
	Matrix ST_GS;		//Матрица направляющих косинусов между СтСК и ГССК
	Matrix ST_SV;		//Матрица направляющих косинусов между СтСК и СвСК

	Matrix RateSV_I;	//Косocимметрическая матрица угловых скоростей СвСК
	Matrix RateGS_I;	//Косocимметрическая матрица угловых скоростей ГССК
	Matrix RateST_I;  //Кососимметрическая матрица угловых скоростей СтСК

	unsigned char N;	//Число компонентов типа double для подсчёта размера структуры
						//Vector3 состоит из 3-х double
						//Matrix33 состоит из 9-ти double
} Matrix_Struct;