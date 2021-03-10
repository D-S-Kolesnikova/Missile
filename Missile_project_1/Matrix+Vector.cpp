#include <math.h>
#include <iostream>
#include "Matrix+Vector.h"



void Vector::TestGetVector(double a1, double a2, double a3)
{
	X = a1; Y = a2; Z = a3;
};

double Vector::GetLength()
{
	return sqrt(X * X + Y * Y + Z * Z);
};

double Vector::GetPitch()
{
	double length = GetLength();
	if (length == 0) return 0;
	return asin(Y / length);
};

double Vector::GetYaw()
{
	return -atan2(Z, X);
};

void Vector::ShowV(Vector pos)
{
	printf("\n%f\n%f\n%f\n", pos.X, pos.Y, pos.Z);
};

void Vector::fShowV(Vector pos, FILE* fMatrix)
{
	fprintf(fMatrix, "%f\n%f\n%f\n", pos.X, pos.Y, pos.Z);
};


Vector  operator + (Vector pos1, Vector pos2)
{
	Vector pos_new;
	pos_new.X = pos1.X + pos2.X;
	pos_new.Y = pos1.Y + pos2.Y;
	pos_new.Z = pos1.Z + pos2.Z;
	return pos_new;
};

Vector  operator - (Vector pos1, Vector pos2)
{
	Vector pos_new;
	pos_new.X = pos1.X - pos2.X;
	pos_new.Y = pos1.Y - pos2.Y;
	pos_new.Z = pos1.Z - pos2.Z;
	return pos_new;
};

Vector  operator * (double number, Vector pos)
{
	Vector pos_new;
	pos_new.X = pos.X * number;
	pos_new.Y = pos.Y * number;
	pos_new.Z = pos.Z * number;
	return pos_new;
};

Vector  operator * (Vector pos, double number)
{
	Vector pos_new;
	pos_new.X = pos.X * number;
	pos_new.Y = pos.Y * number;
	pos_new.Z = pos.Z * number;
	return pos_new;
};


double Matrix:: DET_Matrix() {

	/*3.7. Вычисление определителя матрицы*/
	return   m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
		- m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
		+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

}

void Matrix::Inverse_Matrix() {

	/*3.8. Обращение матрицы, после этой операции матрица Мatrix будет обращенной
	Функция возвращает 0 при успешном обращении и 1 при вырожденной матрице*/

	double det; int i,j;

	det = DET_Matrix();

	// проверка на вырожденность матрицы
	if (det == 0.) {
		std::cout << "Матрица вырожденная!";
	}

	Matrix k;
	k.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
	k.m[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / det;
	k.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;

	k.m[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / det;
	k.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
	k.m[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / det;

	k.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
	k.m[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / det;
	k.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			m[i][j] = k.m[i][j];
		}
	}
}

void Matrix::GetMatrix(double tetta, double psi, double gamma)
{
	m[0][0] = cos(tetta) * cos(psi);
	m[0][1] = -sin(tetta) * cos(psi) * cos(gamma) + sin(psi) * sin(gamma);
	m[0][2] = sin(tetta) * cos(psi) * sin(gamma) + sin(psi) * cos(gamma);

	m[1][0] = sin(tetta);
	m[1][1] = cos(tetta) * cos(gamma);
	m[1][2] = -cos(tetta) * sin(gamma);

	m[2][0] = -cos(tetta) * sin(psi);
	m[2][1] = sin(tetta) * sin(psi) * cos(gamma) + cos(psi) * sin(gamma);
	m[2][2] = -sin(tetta) * sin(psi) * sin(gamma) + cos(psi) * cos(gamma);
};

void Matrix::GetMatrix(double r, double l, double mm, double n)
{
	m[0][0] = r * r + l * l - mm * mm - n * n;
	m[0][1] = 2 * (r * n + l * mm);
	m[0][2] = 2 * (-r * mm + l * n);

	m[1][0] = 2 * (-r * n + l * mm);
	m[1][1] = r * r + mm * mm - n * n - l * l;
	m[1][2] = 2 * (r * l + n * mm);

	m[2][0] = 2 * (r * mm + n * l);
	m[2][1] = 2 * (-r * l + n * mm);
	m[2][2] = r * r + n * n - l * l - mm * mm;

};


void Matrix::GetMatrix(double nu, double mu)
{
	m[0][0] = cos(nu) * cos(mu);
	m[0][1] = cos(nu) * sin(mu);
	m[0][2] = -sin(nu);

	m[1][0] = -sin(mu);
	m[1][1] = cos(mu);
	m[1][2] = 0.0;

	m[2][0] = sin(nu) * cos(mu);
	m[2][1] = sin(nu) * sin(mu);
	m[2][2] = cos(nu);
};

void Matrix::TestGetMatrix(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9)
{
	m[0][0] = a1;
	m[0][1] = a2;
	m[0][2] = a3;

	m[1][0] = a4;
	m[1][1] = a5;
	m[1][2] = a6;

	m[2][0] = a7;
	m[2][1] = a8;
	m[2][2] = a9;
};

void Matrix::GSSK_STSK(double B0, double L0, double PSI0, double B, double L)
{
	m[0][0] = cos(B0) * cos(PSI0) * cos(B) + sin(B) * (cos(PSI0) * cos(L0 - L) * sin(B0) + sin(PSI0) * sin(L0 - L));
	m[0][1] = cos(PSI0) * (-cos(B) * cos(L0 - L) * sin(B0) + cos(B0) * sin(B)) - cos(B) * sin(PSI0) * sin(L0 - L);
	m[0][2] = cos(L0 - L) * sin(PSI0) - cos(PSI0) * sin(B0) * sin(L0 - L);
	m[1][0] = cos(B) * sin(B0) - cos(B0) * cos(L0 - L) * sin(B);
	m[1][1] = cos(B0) * cos(B) * cos(L0 - L) + sin(B0) * sin(B);
	m[1][2] = cos(B0) * sin(L0 - L);
	m[2][0] = -sin(PSI0) * (cos(B0) * cos(B) + cos(L0 - L) * sin(B0) * sin(B)) + cos(PSI0) * sin(B) * sin(L0 - L);
	m[2][1] = -cos(B0) * sin(PSI0) * sin(B) + cos(B) * (cos(L0 - L) * sin(B0) * sin(PSI0) - cos(PSI0) * sin(L0 - L));
	m[2][2] = cos(PSI0) * cos(L0 - L) + sin(B0) * sin(PSI0) * sin(L0 - L);
}

void Matrix::GSSK_SVSK(double psi, double tetta, double gamma)
{
	m[0][0] = cos(tetta) * cos(psi) + sin(gamma) * sin(tetta)* sin(psi);
	m[0][1] = cos(gamma) * sin(tetta);
	m[0][2] = sin(tetta) * cos(psi) * sin(gamma) - sin(psi) * cos(gamma);

	m[1][0] = -cos(psi) * sin(tetta) + sin(gamma) * cos(tetta) * sin(psi);
	m[1][1] = cos(tetta) * cos(gamma);
	m[1][2] = cos(psi) * cos(tetta) * sin(gamma) + sin(tetta) * sin(psi);

	m[2][0] = cos(gamma) * sin(psi);
	m[2][1] = -sin(gamma);
	m[2][2] = cos(psi) * cos(gamma);
}

void Matrix::InitNull()
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) m[i][j] = 0.0;
};





Matrix& operator ~ (Matrix& pos)
{
	double memory;
	for (int i = 0; i < 2; i++)
		for (int j = i + 1; j < 3; j++)
		{
			memory = pos.m[i][j];
			pos.m[i][j] = pos.m[j][i];
			pos.m[j][i] = memory;
		}
	return pos;
};

 Matrix  operator + (Matrix pos1, Matrix pos2)
{
	Matrix pos_new;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos_new.m[i][j] = pos1.m[i][j] + pos2.m[i][j];
	return pos_new;
};




Vector  operator * (Vector pos1, Vector pos2)
{
	Vector pos_new;
	pos_new.X = pos1.X * pos2.X;
	pos_new.Y = pos1.Y * pos2.Y;
	pos_new.Z = pos1.Z * pos2.Z;
	return pos_new;
};

Vector  operator / (Vector pos1, Vector pos2)
{
	Vector pos_new;
	pos_new.X = pos1.X / pos2.X;
	pos_new.Y = pos1.Y / pos2.Y;
	pos_new.Z = pos1.Z / pos2.Z;
	return pos_new;
};

Matrix  operator - (Matrix pos1, Matrix pos2)
{
	Matrix pos_new;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos_new.m[i][j] = pos1.m[i][j] - pos2.m[i][j];
	return pos_new;
};

Matrix  operator * (Matrix pos1, Matrix pos2)
{
	Matrix pos_new;
	pos_new.InitNull();
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				pos_new.m[i][j] += pos1.m[i][k] * pos2.m[k][j];
	return pos_new;
};

 Vector  operator * (Vector v, Matrix pos)
{
	Vector v_new;
	v_new.X = v.X * pos.m[0][0] + v.Y * pos.m[1][0] + v.Z * pos.m[2][0];
	v_new.Y = v.X * pos.m[0][1] + v.Y * pos.m[1][1] + v.Z * pos.m[2][1];
	v_new.Z = v.X * pos.m[0][2] + v.Y * pos.m[1][2] + v.Z * pos.m[2][2];
	return v_new;
};

Vector  operator * (Matrix pos, Vector v)
{
	Vector v_new;
	v_new.X = v.X * pos.m[0][0] + v.Y * pos.m[0][1] + v.Z * pos.m[0][2];
	v_new.Y = v.X * pos.m[1][0] + v.Y * pos.m[1][1] + v.Z * pos.m[1][2];
	v_new.Z = v.X * pos.m[2][0] + v.Y * pos.m[2][1] + v.Z * pos.m[2][2];
	return v_new;
};

Matrix  operator * (double number, Matrix pos)
{
	Matrix pos_new;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) pos_new.m[i][j] = pos.m[i][j] * number;
	return pos_new;
};


Matrix  operator * (Matrix pos, double number)
{
	Matrix pos_new;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) pos_new.m[i][j] = pos.m[i][j] * number;
	return pos_new;
};


Matrix  operator / (Matrix pos, double number)
{
	Matrix pos_new;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) pos_new.m[i][j] = pos.m[i][j] / number;
	return pos_new;
};
 Vector  operator / (Vector pos, double number)
{
	Vector pos_new;
	pos_new.X = pos.X / number;
	pos_new.Y = pos.Y / number;
	pos_new.Z = pos.Z / number;
	return pos_new;
};


 //Векторное умножение (обозначается так [a,b] или так a x b)
 Vector  operator % (Vector a, Vector b)
{
	Vector pos_new;
	pos_new.X = a.Y * b.Z - a.Z * b.Y;
	pos_new.Y = a.Z * b.X - a.X * b.Z;
	pos_new.Z = a.X * b.Y - a.Y * b.X;
	return pos_new;
};

void Matrix::ShowM(Matrix pos)
{
	for (int i = 0; i < 3; i++)
	{
		if (i != 0) printf("\n");
		for (int j = 0; j < 3; j++)
			printf("%f\t", pos.m[i][j]);
	};
	printf("\n");
};

void Matrix::fShowM(Matrix pos, FILE* fMatrix)
{
	for (int i = 0; i < 3; i++)
	{
		if (i != 0) fprintf(fMatrix, "\n");
		for (int j = 0; j < 3; j++)
			fprintf(fMatrix, "%f\t", pos.m[i][j]);
	};
	fprintf(fMatrix, "\n");
};



