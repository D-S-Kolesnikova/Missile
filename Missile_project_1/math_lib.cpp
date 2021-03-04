/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ���������� �������� ����� ���������, �������� ����� ��������� � �������� ����� ��������� �          //
// ���������                                                                                           //
//                                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <conio.h>

#include "math_lib.h"
#include "macros.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ����������                                                                                          //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------
Vector3 set_vector3(double X_value, double Y_value, double Z_value) {
	
	/*1.1. ������������ �������� ��������� ������� - ������������� �������*/
	
	Vector3 Vector;

	Vector.X = X_value;
	Vector.Y = Y_value;
	Vector.Z = Z_value;
	
	return Vector;
	
}

//-------------------------------------------------------------------------------------------------------
void set_element_vector3(Vector3 *Vector, double Value, unsigned short int Index) {

	/* 1.2. ���������� �������� Value �������� Index =(_X,_Y,_Z) == (0,1,2) ������� Vector */
	
	((double *)Vector)[Index] = Value;
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
double get_length_vector3(Vector3 Vector) {

	/*1.3.��������� ����� �������*/
	
	return sqrt(_SQR(Vector.X) + _SQR(Vector.Y) + _SQR(Vector.Z));
	
}

//-------------------------------------------------------------------------------------------------------
double get_lambda_vector3(Vector3 Vector) {

	/*1.4. ����������� ����������� �����*/
	// ������ �������
	
	return atan2(Vector.Y,sqrt(_SQR(Vector.X) + _SQR(Vector.Z)));
	
}

//-------------------------------------------------------------------------------------------------------
double get_fi_vector3(Vector3 Vector) {
	
	/*1.4. ����������� ����������� �����*/
	// ������ ������
	
	return -atan2(Vector.Z, Vector.X);
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 set_zero_vector3(void) {

	/*1.5. ������������� �������� ������*/

	Vector3 vector;

	vector.X = 0.;
	vector.Y = 0.;
	vector.Z = 0.;

	return vector;

}

//-------------------------------------------------------------------------------------------------------
Vector3 set_unitary_vector3(void) {

	/*1.6. ������������� ���������� ������*/

	Vector3 vector;

	vector.X = 1.;
	vector.Y = 1.;
	vector.Z = 1.;

	return vector;

}

//-------------------------------------------------------------------------------------------------------
double get_element_vector3(Vector3 Vector, unsigned short int Index) {

	/*1.7. ��������� �������� �������� ������� Vector � �������� Index = (_X,_Y,_Z) == (0,1,2) */

	return ((double *)&Vector)[Index];

}

//-------------------------------------------------------------------------------------------------------
void set_matrix33(Matrix33 Matrix, double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33) {

	/*2.1. ������������ �������� ���� ��������� ������� - ������������� �������
	   Matrix - ������� ���� Matrix33
	   a_ij - (i,j)-�� ������� ������� Matrix*/

	Matrix[0][0] = a11;
	Matrix[0][1] = a12;
	Matrix[0][2] = a13;
	Matrix[1][0] = a21;
	Matrix[1][1] = a22;
	Matrix[1][2] = a23;
	Matrix[2][0] = a31;
	Matrix[2][1] = a32;
	Matrix[2][2] = a33;
	   
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void set_element_matrix33(Matrix33 Matrix, unsigned short int i, unsigned short int j, double value) {

	/*2.2. ������������ �������� value (i,j)-��� �������� ������� Matrix
	   Matrix - ������� ���� Matrix33*/

	Matrix[i][j] = value;
	   
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void set_svsk_measure_matrix33(Matrix33 Matrix, double nu, double mu) {
	
	/*2.3. ������� ��������*/
	//������� �������� �� ��������� � ������������� �.�.
	
	Matrix[0][0] =   cos(nu)*cos(mu);
	Matrix[0][1] =   cos(nu)*sin(mu);
	Matrix[0][2] = - sin(nu);
	Matrix[1][0] = - sin(mu);
	Matrix[1][1] =   cos(mu);
	Matrix[1][2] =   0.;
	Matrix[2][0] =   sin(nu)*cos(mu);
	Matrix[2][1] =   sin(nu)*sin(mu);
	Matrix[2][2] =   cos(nu);

	return;
	
}

//-------------------------------------------------------------------------------------------------------
void set_nearth_stsk_matrix33(Matrix33 Matrix_GS_ST_l, double B42_l, double L42_l, double B42_ST_l, double L42_ST_l, double Ag_l) {

	/*2.3. ������� ��������*/
	//������� �������� �� ���������� ������ � ��������� �.�.

	Matrix33 Matrix_GS_SK42u, Matrix_SK42u_ST;
	set_matrix33(Matrix_GS_SK42u,
								-sin(B42_l)*cos(L42_l), cos(B42_l)*cos(L42_l), -sin(L42_l),
								-sin(B42_l)*sin(L42_l), cos(B42_l)*sin(L42_l),  cos(L42_l),
											cos(B42_l),            sin(B42_l),        0.0);
	set_matrix33(Matrix_SK42u_ST,
								-sin(B42_ST_l)*cos(L42_ST_l)*cos(Ag_l)-sin(L42_ST_l)*sin(Ag_l), -sin(B42_ST_l)*sin(L42_ST_l)*cos(Ag_l)+cos(L42_ST_l)*sin(Ag_l),   cos(B42_ST_l)*cos(Ag_l),
																   cos(B42_ST_l)*cos(L42_ST_l),                                    cos(B42_ST_l)*sin(L42_ST_l),             sin(B42_ST_l),
								-sin(L42_ST_l)*cos(Ag_l)+sin(B42_ST_l)*cos(L42_ST_l)*sin(Ag_l),  cos(L42_ST_l)*cos(Ag_l)+sin(B42_ST_l)*sin(L42_ST_l)*sin(Ag_l), -cos(B42_ST_l)*sin(Ag_l));

	mult_matrix33(Matrix_GS_ST_l,Matrix_SK42u_ST,Matrix_GS_SK42u);
	return;

}

//-------------------------------------------------------------------------------------------------------
void set_svsk_nearth_matrix33(Matrix33 Matrix, double theta, double psi, double gamma) {

	/*2.3. ������� ��������*/
	//������� �������� �� ��������� � ���������� ������ �.�.(������������������ ���������: ���� (�� ������� �������) -> ���� -> ������)

	Matrix[0][0] =   cos(theta)*cos(psi) - sin(theta)*sin(gamma)*sin(psi);
	Matrix[0][1] = - sin(theta)*cos(psi) - cos(theta)*sin(gamma)*sin(psi);  
	Matrix[0][2] = - cos(gamma)*sin(psi);
	Matrix[1][0] =   sin(theta)*cos(gamma);
	Matrix[1][1] =   cos(theta)*cos(gamma);
	Matrix[1][2] = - sin(gamma);  
	Matrix[2][0] =   sin(theta)*sin(gamma)*cos(psi) + cos(theta)*sin(psi);
	Matrix[2][1] =   cos(theta)*sin(gamma)*cos(psi) - sin(theta)*sin(psi);
	Matrix[2][2] =	 cos(gamma)*cos(psi);

	return;

}

//-------------------------------------------------------------------------------------------------------
void set_sksk_svsk_matrix33(Matrix33 Matrix, double alpha, double beta) {

	/*2.3 ������� ��������*/
	//������� �������� �� ���������� � ��������� �.�.

	Matrix[0][0] =  cos(alpha)*cos(beta);
	Matrix[0][1] =  sin(alpha);
	Matrix[0][2] = -cos(alpha)*sin(beta);
	Matrix[1][0] = -sin(alpha)*cos(beta);
	Matrix[1][1] =  cos(alpha);
	Matrix[1][2] =  sin(alpha)*sin(beta);
	Matrix[2][0] =  sin(beta);
	Matrix[2][1] =  0;
	Matrix[2][2] =  cos(beta);

	return;

}

//-------------------------------------------------------------------------------------------------------
void set_nearth_svsk_rg_matrix33(Matrix33 Matrix, double r, double l, double mm, double n) {

	/*2.3. ������� ��������*/
	//������� �������� �� ���������� ������ � ��������� �.�.(� ���������� �������-����������)

	Matrix[0][0] = (r)*(r) + (l)*(l) - (mm)*(mm) - (n)*(n);
	Matrix[0][1] = 2*((r)*(n) + (l)*(mm));
	Matrix[0][2] = 2*(- (r)*(mm) + (l)*(n));
	Matrix[1][0] = 2*(- (r)*(n)+(l)*(mm));
	Matrix[1][1] = (r)*(r) + (mm)*(mm) - (n)*(n) - (l)*(l);
	Matrix[1][2] = 2*((r)*(l) + (n)*(mm));
	Matrix[2][0] = 2*((r)*(mm) + (n)*(l));
	Matrix[2][1] = 2*(- (r)*(l) + (n)*(mm));
	Matrix[2][2] = (r)*(r) + (n)*(n) - (l)*(l) - (mm)*(mm);
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void set_zero_matrix33(Matrix33 Matrix) {

	/*2.4. ������������� ������� �������*/

	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			Matrix[i][j] = 0.;
		}
	}

	return;

}

//-------------------------------------------------------------------------------------------------------
void set_unitary_matrix33(Matrix33 Matrix) {

	/*2.5. ������������� ��������� �������*/

	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			Matrix[i][j] = 0.;
		}
		Matrix[i][i] = 1.;
	}

	return;

}

//-------------------------------------------------------------------------------------------------------
Vector3 add_vector3(Vector3 vector_1, Vector3 vector_2) {

	/*3.1. c������� ��������
	res_vector = vector_1 + vector_2 */
	
	Vector3 vector_res;

	vector_res.X = vector_1.X + vector_2.X;
	vector_res.Y = vector_1.Y + vector_2.Y;
	vector_res.Z = vector_1.Z + vector_2.Z;
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 sub_vector3(Vector3 vector_1, Vector3 vector_2) {

	/*3.2. ��������� ��������
	res_vector = vector_1 - vector_2 */
	
	Vector3 vector_res;

	vector_res.X = vector_1.X - vector_2.X;
	vector_res.Y = vector_1.Y - vector_2.Y;
	vector_res.Z = vector_1.Z - vector_2.Z;
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 mult_vector3(Vector3 vector_1, Vector3 vector_2) {

	/*3.3. ������������ ��������� ��������
	res_vector = vector_1 * vector_2 */
	
	Vector3 vector_res;

	vector_res.X = vector_1.X * vector_2.X;
	vector_res.Y = vector_1.Y * vector_2.Y;
	vector_res.Z = vector_1.Z * vector_2.Z;
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 div_vector3(Vector3 vector_1, Vector3 vector_2) {

	/*3.4. ������������ ������� ��������
	res_vector = vector_1 / vector_2 */
	
	Vector3 vector_res;

	vector_res.X = vector_1.X / vector_2.X;
	vector_res.Y = vector_1.Y / vector_2.Y;
	vector_res.Z = vector_1.Z / vector_2.Z;
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
void equate_matrix33(Matrix33 matrix_dest, Matrix33 matrix_source) {
	
	/*3.5. ������������� ���� ������
	Matrix_dest = Matrix_source */

	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			matrix_dest[i][j] = matrix_source[i][j];
		}
	}
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void transpose_matrix33(Matrix33 Matrix) {

	/*3.6. ���������������� �������, ����� ���� �������� ������� �atrix ����� ���������������*/
	
	unsigned short int i,j;
	double memory;

	for(i = 0; i < 2; i++) {
		for(j = i + 1; j < 3; j++) {
			memory = Matrix[i][j];
			Matrix[i][j] = Matrix[j][i];
			Matrix[j][i] = memory;
		}
	}

	return;
}

//-------------------------------------------------------------------------------------------------------
double det_matrix33(Matrix33 Matrix) {

	/*3.7. ���������� ������������ �������*/
		
	return   Matrix[0][0]*(Matrix[1][1]*Matrix[2][2] - Matrix[1][2]*Matrix[2][1])
		   - Matrix[0][1]*(Matrix[1][0]*Matrix[2][2] - Matrix[1][2]*Matrix[2][0])
		   + Matrix[0][2]*(Matrix[1][0]*Matrix[2][1] - Matrix[1][1]*Matrix[2][0]);
	
}

//-------------------------------------------------------------------------------------------------------
unsigned short int inverse_matrix33(Matrix33 Matrix) {

	/*3.8. ��������� �������, ����� ���� �������� ������� �atrix ����� ����������
	������� ���������� 0 ��� �������� ��������� � 1 ��� ����������� �������*/

	unsigned short int i,j;
	double det;
	Matrix33 temp_matrix;

	det = det_matrix33(Matrix);

	// �������� �� ������������� �������
	if (det == 0.) {
		return 1;
	}

	temp_matrix[0][0] =   (Matrix[1][1]*Matrix[2][2] - Matrix[1][2]*Matrix[2][1])/det;
	temp_matrix[1][0] = - (Matrix[1][0]*Matrix[2][2] - Matrix[1][2]*Matrix[2][0])/det;
	temp_matrix[2][0] =   (Matrix[1][0]*Matrix[2][1] - Matrix[1][1]*Matrix[2][0])/det;

	temp_matrix[0][1] = - (Matrix[0][1]*Matrix[2][2] - Matrix[0][2]*Matrix[2][1])/det;
	temp_matrix[1][1] =   (Matrix[0][0]*Matrix[2][2] - Matrix[0][2]*Matrix[2][0])/det;
	temp_matrix[2][1] = - (Matrix[0][0]*Matrix[2][1] - Matrix[0][1]*Matrix[2][0])/det;

	temp_matrix[0][2] =   (Matrix[0][1]*Matrix[1][2] - Matrix[0][2]*Matrix[1][1])/det;
	temp_matrix[1][2] = - (Matrix[0][0]*Matrix[1][2] - Matrix[0][2]*Matrix[1][0])/det;
	temp_matrix[2][2] =   (Matrix[0][0]*Matrix[1][1] - Matrix[0][1]*Matrix[1][0])/det;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			Matrix[i][j] = temp_matrix[i][j];
		}
	}
	
	return 0;
	
}

//-------------------------------------------------------------------------------------------------------
void add_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2) {

	/*3.9. �������� ������ 
	res_matrix = matrix_1 + matrix_2 */

	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			res_matrix[i][j] = matrix_1[i][j] + matrix_2[i][j];
		}
	}
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void sub_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2) {

	/*3.10. ��������� ������ 
	res_matrix = matrix_1 - matrix_2 */

	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			res_matrix[i][j] = matrix_1[i][j] - matrix_2[i][j];
		}
	}
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void mult_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2) {

	/*3.11. ��������� ������ 
	res_matrix = matrix_1 * matrix_2 */

	unsigned short int i,j,k;

	res_matrix[0][0] = 0.;
	res_matrix[0][1] = 0.;
	res_matrix[0][2] = 0.;
	res_matrix[1][0] = 0.;
	res_matrix[1][1] = 0.;
	res_matrix[1][2] = 0.;
	res_matrix[2][0] = 0.;
	res_matrix[2][1] = 0.;
	res_matrix[2][2] = 0.;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			for(k = 0; k < 3; k++) {
				res_matrix[i][j] += matrix_1[i][k] * matrix_2[k][j];
			}
		}
	}
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 mult_vector3_matrix33(Vector3 vector, Matrix33 matrix) {

	/*3.12. ��������� �������-������ �� �������. ��������� - ������-������ 
	res_vector = vector*matrix */
	
	Vector3 vector_res;

	vector_res.X = vector.X*matrix[0][0] + vector.Y*matrix[1][0] + vector.Z*matrix[2][0];
	vector_res.Y = vector.X*matrix[0][1] + vector.Y*matrix[1][1] + vector.Z*matrix[2][1];
	vector_res.Z = vector.X*matrix[0][2] + vector.Y*matrix[1][2] + vector.Z*matrix[2][2];
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 mult_matrix33_vector3(Matrix33 matrix, Vector3 vector) {
	/*3.13. ��������� ������� �� ������-�������. ��������� - ������-�������
	res_vector = matrix*vector */
	
	Vector3 vector_res;

	vector_res.X = vector.X*matrix[0][0] + vector.Y*matrix[0][1] + vector.Z*matrix[0][2];
	vector_res.Y = vector.X*matrix[1][0] + vector.Y*matrix[1][1] + vector.Z*matrix[1][2];
	vector_res.Z = vector.X*matrix[2][0] + vector.Y*matrix[2][1] + vector.Z*matrix[2][2];
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
void mult_num_matrix33(Matrix33 matrix_res, double number, Matrix33 matrix) {
	
	/*3.14. ��������� ����� �� ������� 
	res_matrix = number*matrix */

	unsigned short int i,j;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			matrix_res[i][j] = matrix[i][j]*number;
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 mult_num_vector3(double number, Vector3 vector) {

	/*3.15. ��������� ����� �� ������ 
	res_vector = vector*number */
	
	Vector3 vector_res;

	vector_res.X = vector.X*number;
	vector_res.Y = vector.Y*number;
	vector_res.Z = vector.Z*number;
	
	return vector_res;
	
}

//-------------------------------------------------------------------------------------------------------
Vector3 multvect_vector3(Vector3 vector_1, Vector3 vector_2) {

	/*3.16. ��������� ��������� (������������ ��� [vector_1,vector_1] ��� ��� vector_1 x vector_1) */
	
	Vector3 vector_res;

	vector_res.X = vector_1.Y*vector_2.Z - vector_1.Z*vector_2.Y;
	vector_res.Y = vector_1.Z*vector_2.X - vector_1.X*vector_2.Z;
	vector_res.Z = vector_1.X*vector_2.Y - vector_1.Y*vector_2.X;
	
	return vector_res;
	
}


Vector3 add_num_vector3(double number, Vector3 vector) {

	/*3.17. �������������� �������� ������� � �����*/

	Vector3 vector_res;

	vector_res.X = vector.X + number;
	vector_res.Y = vector.Y + number;
	vector_res.Z = vector.Z + number;
	
	return vector_res;

}

//-------------------------------------------------------------------------------------------------------
void print_console_vector3(Vector3 Vector) {

	/*4.1. ������ ������� �����������*/
	// ������ �� ������� (�����)
	
	printf("\n%f\n%f\n%f\n", Vector.X, Vector.Y, Vector.Z);
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void print_file_vector3(Vector3 Vector, FILE *file_descriptor) {
	
	/*4.1. ������ ������� �����������*/
	// ������ � ���� � ������������ file_descriptor
	
	fprintf(file_descriptor,"%f\n%f\n%f\n", Vector.X, Vector.Y, Vector.Z);
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void print_console_matrix33(Matrix33 matrix) {
	
	/*4.2. ������ ������� ����������� � ���������������*/
	// ������ �� ������� (�����)
	
	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		if(i) {
			printf("\n");
		}
		for(j = 0; j < 3; j++) {
			printf("%f\t",matrix[i][j]);
		}
	}
	printf("\n");
	
	return;
	
}

//-------------------------------------------------------------------------------------------------------
void print_file_matrix33(Matrix33 matrix, FILE *file_descriptor) {
	
	/*4.2. ������ ������� ����������� � ���������������*/
	// ������ � ���� � ������������ file_descriptor
	
	unsigned short int i,j;

	for(i = 0; i < 3; i++) {
		if(i) {
			fprintf(file_descriptor,"\n");
		}
		for(j = 0; j < 3; j++) {
			fprintf(file_descriptor,"%f\t",matrix[i][j]);
		}
	}
	
	return;
	
}

double pow10(int number) {
	/*5.1. ���������� ����� 10 � ����� �������*/
	double Pow10 = 1.0;
	if (number > 0)
		for (int i = 0; i < number; i++)
			Pow10 *= 10.0;
	else
		for (int i = 0; i > number; i--)
			Pow10 /= 10.0;
	return Pow10;
}

void set_GSSK_STSK_matrix33(Matrix33 Matrix, double B0, double L0, double PSI0, double B, double L) {

	/* ������� �������� */
	//������� �������� �� ���� � ����

	Matrix[0][0] = cos(B0) * cos(PSI0) * cos(B)+ sin(B)*(cos(PSI0)*cos(L0-L)*sin(B0)+sin(PSI0)*sin(L0-L));
	Matrix[0][1] = cos(PSI0) * (-cos(B) * cos(L0 - L) * sin(B0) + cos(B0) * sin(B)) - cos(B) * sin(PSI0) * sin(L0 - L);
	Matrix[0][2] = cos(L0 - L) * sin(PSI0) - cos(PSI0) * sin(B0) * sin(L0 - L);
	Matrix[1][0] = cos(B) * sin(B0) - cos(B0) * cos(L0 - L)*sin(B);
	Matrix[1][1] = cos(B0) * cos(B) * cos(L0 - L) + sin(B0) * sin(B);
	Matrix[1][2] = cos(B0) * sin(L0 - L);
	Matrix[2][0] = -sin(PSI0) * (cos(B0) * cos(B) + cos(L0 - L) * sin(B0) * sin(B)) + cos(PSI0) * sin(B) * sin(L0 - L);
	Matrix[2][1] = -cos(B0) * sin(PSI0) * sin(B) + cos(B) * (cos(L0 - L) * sin(B0) * sin(PSI0) - cos(PSI0) * sin(L0 - L));
	Matrix[2][2] = cos(PSI0) * cos(L0 - L) + sin(B0) * sin(PSI0) * sin(L0 - L);

	return;

}