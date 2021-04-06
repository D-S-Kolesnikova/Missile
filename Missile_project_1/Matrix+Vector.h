//----------------------------Array.h (begin)--------------------------
//Определения операций между матрицами, операций между 
//векторами и операций между матрицами и векторами
#include <stdio.h>

#ifndef TYPES_H
#define TYPES_H
struct Vector

{
public:

	double X;
	double Y;
	double Z;


	Vector(double vx = 0.0, double vy = 0.0, double vz = 0.0) // Это конструктор класса
	{
		X = vx; Y = vy; Z = vz;
	};

	void TestGetVector(double a1, double a2, double a3);
	double GetLength();
	double GetPitch();
	double GetYaw();
	void ShowV(Vector pos);
	void fShowV(Vector pos, FILE* fMatrix);
	friend Vector  operator + (Vector pos1, Vector pos2);
	friend Vector  operator - (Vector pos1, Vector pos2);
	friend Vector  operator * (double number, Vector pos);
	friend Vector  operator * (Vector pos, double number);
	friend Vector  operator % (Vector a, Vector b);

};
class Matrix
{
protected:
		double m[3][3];

public:
	double DET_Matrix();
	void Inverse_Matrix();
	void GetMatrix(double r, double l, double mm, double n);
	void GetMatrix(double nu, double mu);
	void GetMatrix(double tetta, double psi, double gamma);
	void TestGetMatrix(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9);
	void GSSK_STSK(double B0, double L0, double Psi0, double B, double L);
	void GSSK_SVSK(double psi, double tetta, double gamma);
	void Matrix_B(double fi_, double hi_, double r_);
	void InitNull();
	friend Matrix& operator ~ (Matrix& pos);
	friend Matrix  operator + (Matrix pos1, Matrix pos2);
	friend Vector  operator * (Vector pos1, Vector pos2);
	friend Vector  operator / (Vector pos1, Vector pos2);
	friend Matrix  operator - (Matrix pos1, Matrix pos2);
	friend Matrix  operator * (Matrix pos1, Matrix pos2);
	friend Vector  operator * (Vector v, Matrix pos);
	friend Vector  operator * (Matrix pos, Vector v);
	friend Matrix  operator * (double number, Matrix pos);
	friend Matrix  operator * (Matrix pos, double number);
	friend Matrix  operator / (Matrix pos, double number);
	friend Vector  operator / (Vector pos, double number);
	friend Vector  operator % (Vector a, Vector b);
	void ShowM(Matrix pos);
	void fShowM(Matrix pos, FILE* fMatrix);
};
#endif

/*void set_GSSK_STSK_matrix33(Matrix Matrix, double B0, double L0, double PSI0, double B, double L) {

	/* Матрицы перехода 
	//Матрица перехода от ГССК к СтСК
	Matrix.m[0][0] = cos(B0) * cos(PSI0) * cos(B) + sin(B) * (cos(PSI0) * cos(L0 - L) * sin(B0) + sin(PSI0) * sin(L0 - L));
	Matrix.m[0][1] = cos(PSI0) * (-cos(B) * cos(L0 - L) * sin(B0) + cos(B0) * sin(B)) - cos(B) * sin(PSI0) * sin(L0 - L);
	Matrix.m[0][2] = cos(L0 - L) * sin(PSI0) - cos(PSI0) * sin(B0) * sin(L0 - L);
	Matrix.m[1][0] = cos(B) * sin(B0) - cos(B0) * cos(L0 - L) * sin(B);
	Matrix.m[1][1] = cos(B0) * cos(B) * cos(L0 - L) + sin(B0) * sin(B);
	Matrix.m[1][2] = cos(B0) * sin(L0 - L);
	Matrix.m[2][0] = -sin(PSI0) * (cos(B0) * cos(B) + cos(L0 - L) * sin(B0) * sin(B)) + cos(PSI0) * sin(B) * sin(L0 - L);
	Matrix.m[2][1] = -cos(B0) * sin(PSI0) * sin(B) + cos(B) * (cos(L0 - L) * sin(B0) * sin(PSI0) - cos(PSI0) * sin(L0 - L));
	Matrix.m[2][2] = cos(PSI0) * cos(L0 - L) + sin(B0) * sin(PSI0) * sin(L0 - L);
	return;}*/