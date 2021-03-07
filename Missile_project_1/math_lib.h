/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ����������� �������� ����� ���������, �������� ����� ��������� � �������� ����� ��������� �         //
// ���������                                                                                           //
//                                                                                                     //
// ���������� ����������� � ������� ��������, ���� �������. ��� ������ ������������ �������            //
// ����������� ���������������� ������ _USE_MACROS, �������� �������� ���� ����. �������������         //
// ������������ ����� � ������� �������� �������� ����������                                           //
//                                                                                                     //
// �������� �� 02.05.2012                                                                              //
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _MATH_LIB_H_
#define _MATH_LIB_H_

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// �������� ��������� � �������                                                                        //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//���������� �� math.h, ������� ��������� ������������� ���. �������� (� ���������� ����� ��)
#define _USE_MATH_DEFINES	

#include <stdio.h>
#include <math.h>

#define PI M_PI
#define TO_DEG (180/PI)
#define TO_RAD (PI/180)

/* ������ ��� �������� ����������. ��������� �������, ������������-�� � ��������� ������� ��������,
���� ������� */
#define _USE_MACROS 0 // ��� ��������� 0 ������������ �������, ��� ������ ������������� �������� - �������

/* �������������� ������� */
#define _SQR(w) ((w)*(w)) // ������� �����

/* ������� �������� */

// ��������� ������ � ��������� - ������ � ��������� ������� �� �������������� ����������
#define _Index_X 0
#define _Index_Y 1
#define _Index_Z 2

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ����������� �����                                                                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*������ � 3 ������������*/
struct vector3 {
	double X;
	double Y;
	double Z;
};

// ��� ������� �� ���� ���������
typedef struct vector3 Vector3;
// ��� ������� 3�3
typedef double Matrix33[3][3];

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ���������� � ������� ���������� //                                                                  //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// !!! � ����������� �������� ������������ �������� ����������������� ���������� "�������", ��� ���������
// �������� � ������������� ������ �������� � ������, ���� �������� ����������

/*1. ������������� �������� � ��������� �� �������������*/
/*1.1. ������������ �������� ��������� ������� - ������������� �������*/
#define SET_VECTOR3(Vector,X_value,Y_value,Z_value) (Vector.X = X_value, \
													 Vector.Y = Y_value, \
													 Vector.Z = Z_value)
/* 1.2. ���������� �������� Value �������� Index =(_X,_Y,_Z) == (0,1,2) ������� Vector */
#define SET_ELEMENT_VECTOR3(Vector,Value,Index) ((double *)Vector)[Index] = Value
/*1.3.��������� ����� �������*/
#define GET_LENGTH_VECTOR3(Vector) sqrt(_SQR(Vector.X) + _SQR(Vector.Y) + _SQR(Vector.Z))
/*1.4. ����������� ����������� �����*/
#define	GET_LAMBDA_VECTOR3(Vector)   atan2(Vector.Y,sqrt(_SQR(Vector.X) + _SQR(Vector.Z))) // ������ �������
#define	GET_FI_VECTOR3(Vector)     - atan2(Vector.Z,Vector.X)                              // ������ ������
/*1.5. ������������� �������� �������*/
#define SET_ZERO_VECTOR3(vector) SET_VECTOR3(vector,0.,0.,0.)
/*1.6. ������������� ���������� �������*/
#define SET_UNITARY_VECTOR3(vector) SET_VECTOR3(vector,1.,1.,1.)
/*1.7. ��������� �������� �������� ������� Vector � �������� Index = (_X,_Y,_Z) == (0,1,2) */
#define GET_ELEMENT_VECTOR3(Vector, Index) ((double *)Vector)[Index]

/*2. ������������� ������*/
/*2.1. ������������ �������� ���� ��������� ������� - ������������� �������
	   Matrix - ������� ���� Matrix33
	   a_ij - (i,j)-�� ������� ������� Matrix*/
#define	SET_MATRIX33(Matrix,a11,a12,a13,a21,a22,a23,a31,a32,a33) (Matrix[0][0] = a11, \
																  Matrix[0][1] = a12, \
																  Matrix[0][2] = a13, \
																  Matrix[1][0] = a21, \
																  Matrix[1][1] = a22, \
																  Matrix[1][2] = a23, \
																  Matrix[2][0] = a31, \
																  Matrix[2][1] = a32, \
																  Matrix[2][2] = a33)
/*2.2. ������������ �������� value (i,j)-��� �������� ������� Matrix
	   Matrix - ������� ���� Matrix33*/
#define SET_ELEMENT_MATRIX33(Matrix,i,j,value) Matrix[i][j] = value
/*2.3. ������� ��������*/
//������� �������� �� ��������� � ������������� �.�.
#define	SET_SVSK_MEASURE_MATRIX33(Matrix,nu,mu)	(Matrix[0][0] =   cos(nu)*cos(mu), \
												 Matrix[0][1] =   cos(nu)*sin(mu), \
												 Matrix[0][2] = - sin(nu),         \
												 Matrix[1][0] = - sin(mu),         \
												 Matrix[1][1] =   cos(mu),         \
												 Matrix[1][2] =   0.0,             \
												 Matrix[2][0] =   sin(nu)*cos(mu), \
												 Matrix[2][1] =   sin(nu)*sin(mu), \
												 Matrix[2][2] =   cos(nu))
//������� �������� �� ���������� � ��������� �.�.
#define SET_SKSK_SVSK_MATRIX33(Matrix,alpha,beta) (Matrix[0][0] =  cos(alpha)*cos(beta), \
												   Matrix[0][1] =  sin(alpha) 		   , \
												   Matrix[0][2] = -cos(alpha)*sin(beta), \
												   Matrix[1][0] = -sin(alpha)*cos(beta), \
												   Matrix[1][1] =  cos(alpha)		   , \
												   Matrix[1][2] =  sin(alpha)*sin(beta), \
												   Matrix[2][0] =  sin(beta)		   , \
												   Matrix[2][1] =  0				   , \
												   Matrix[2][2] =  cos(beta))
//������� �������� �� ��������� � ���������� ������ �.�.(� ����� ������)
#define SET_SVSK_NEARTH_MATRIX33(Matrix,theta,psi,gamma) (Matrix[0][0] =   cos(theta)*cos(psi),                                  \
														  Matrix[0][1] = - sin(theta)*cos(psi)*cos(gamma) + sin(psi)*sin(gamma), \
														  Matrix[0][2] =   sin(theta)*cos(psi)*sin(gamma) + sin(psi)*cos(gamma), \
														  Matrix[1][0] =   sin(theta),                                           \
														  Matrix[1][1] =   cos(theta)*cos(gamma),                                \
														  Matrix[1][2] = - cos(theta)*sin(gamma),                                \
														  Matrix[2][0] = - cos(theta)*sin(psi),                                  \
														  Matrix[2][1] =   sin(theta)*sin(psi)*cos(gamma) + cos(psi)*sin(gamma), \
														  Matrix[2][2] = - sin(theta)*sin(psi)*sin(gamma) + cos(psi)*cos(gamma))
//������� �������� �� ���������� ������ � ��������� �.�.(� ���������� �������-����������)
#define	SET_NEARTH_SVSK_RG_MATRIX33(Matrix,r,l,mm,n) (Matrix[0][0] = (r)*(r) + (l)*(l) - (mm)*(mm) - (n)*(n),  \
													  Matrix[0][1] = 2*((r)*(n) + (l)*(mm)),                   \
													  Matrix[0][2] = 2*(- (r)*(mm) + (l)*(n)),                 \
													  Matrix[1][0] = 2*(- (r)*(n)+(l)*(mm)),                   \
													  Matrix[1][1] = (r)*(r) + (mm)*(mm) - (n)*(n) - (l)*(l),  \
													  Matrix[1][2] = 2*((r)*(l) + (n)*(mm)),                   \
													  Matrix[2][0] = 2*((r)*(mm) + (n)*(l)),                   \
													  Matrix[2][1] = 2*(- (r)*(l) + (n)*(mm)),                 \
													  Matrix[2][2] = (r)*(r) + (n)*(n) - (l)*(l) - (mm)*(mm))
/*2.4. ������������� ������� �������*/
#define SET_ZERO_MATRIX33(Matrix) do {                           \
									  unsigned short int i,j;    \
									  for(i = 0; i < 3; i++)     \
										  for(j = 0; j < 3; j++) \
											  Matrix[i][j] = 0.; \
								  } while(0) // �������� do-while ��������, ����� �������� ���������
											 // "����� � �������" � �������� ����������, ���������� "else".
											 // ����������� ���� ���
/*2.5. ������������� ��������� �������*/
#define SET_UNITARY_MATRIX33(Matrix) do {                           \
										 unsigned short int i,j;    \
										 for(i = 0; i < 3; i++) {   \
											 for(j = 0; j < 3; j++) \
												Matrix[i][j] = 0.;  \
											Matrix[i][i] = 1.;      \
										 }                          \
									 } while(0) // �������� do-while ��������, ����� �������� ���������
												// "����� � �������" � �������� ����������, ���������� "else".
												// ����������� ���� ���
/*3. �������� � ��������� � ���������*/
/*3.1. c������� ��������
res_vector = vector_1 + vector_2 */
#define	ADD_VECTOR3(res_vector,vector_1,vector_2) (res_vector.X = vector_1.X + vector_2.X, \
												   res_vector.Y = vector_1.Y + vector_2.Y, \
												   res_vector.Z = vector_1.Z + vector_2.Z)
/*3.2. ��������� ��������
res_vector = vector_1 - vector_2 */
#define	SUB_VECTOR3(res_vector,vector_1,vector_2) (res_vector.X = vector_1.X - vector_2.X, \
												   res_vector.Y = vector_1.Y - vector_2.Y, \
												   res_vector.Z = vector_1.Z - vector_2.Z)
/*3.3. ������������ ��������� ��������
res_vector = vector_1 * vector_2 */
#define	MULT_VECTOR3(res_vector,vector_1,vector_2) (res_vector.X = vector_1.X * vector_2.X, \
												    res_vector.Y = vector_1.Y * vector_2.Y, \
												    res_vector.Z = vector_1.Z * vector_2.Z)
/*3.4. ������������ ������� ��������
res_vector = vector_1 / vector_2 */
#define	DIV_VECTOR3(res_vector,vector_1,vector_2) (res_vector.X = vector_1.X / vector_2.X, \
												   res_vector.Y = vector_1.Y / vector_2.Y, \
												   res_vector.Z = vector_1.Z / vector_2.Z)
/*3.5. ������������� ���� ������
Matrix_dest = Matrix_source */
#define	EQUATE_MATRIX33(Matrix_dest,Matrix_source) do {                                                 \
													   unsigned short int i,j;                          \
													   for(i = 0; i < 3; i++)                           \
														   for(j = 0; j < 3; j++)                       \
															   Matrix_dest[i][j] = Matrix_source[i][j]; \
												   } while(0) // �������� do-while ��������, ����� �������� ���������
															  // "����� � �������" � �������� ����������, ���������� "else".
															  // ����������� ���� ���

/*3.6. ���������������� �������, ����� ���� �������� ������� �atrix ����� ���������������*/
#define	TRANSPOSE_MATRIX33(Matrix) do {                                     \
									   double memory;                       \
									   unsigned short int i,j;              \
									   for(i = 0; i < 2; i++)               \
										   for(j = i + 1; j < 3; j++) {     \
											   memory = Matrix[i][j];       \
											   Matrix[i][j] = Matrix[j][i]; \
											   Matrix[j][i] = memory;       \
										   }                                \
								   } while(0) // �������� do-while ��������, ����� �������� ���������
											  // "����� � �������" � �������� ����������, ���������� "else".
											  // ����������� ���� ���
/*3.7. ���������� ������������ �������*/
#define DET_MATRIX33(Matrix)   Matrix[0][0]*(Matrix[1][1]*Matrix[2][2] - Matrix[1][2]*Matrix[2][1]) \
							 - Matrix[0][1]*(Matrix[1][0]*Matrix[2][2] - Matrix[1][2]*Matrix[2][0]) \
							 + Matrix[0][2]*(Matrix[1][0]*Matrix[2][1] - Matrix[1][1]*Matrix[2][0])

/*3.8. ��������� �������, ����� ���� �������� ������� �atrix ����� ����������
	������� ���������� 0 ��� �������� ��������� � 1 ��� ����������� �������*/
#define INVERSE_MATRIX33(Matrix) do {                                                                                 \
								   double det;                                                                        \
								   Matrix33 temp_matrix;                                                              \
								   det = DET_MATRIX33(Matrix);                                                        \
								   if (det == 0.)  {                                                                  \
									   break;                                                                         \
								   }                                                                                  \
								   temp_matrix[0][0] =   (Matrix[1][1]*Matrix[2][2] - Matrix[1][2]*Matrix[2][1])/det; \
								   temp_matrix[1][0] = - (Matrix[1][0]*Matrix[2][2] - Matrix[1][2]*Matrix[2][0])/det; \
								   temp_matrix[2][0] =   (Matrix[1][0]*Matrix[2][1] - Matrix[1][1]*Matrix[2][0])/det; \
								   temp_matrix[0][1] = - (Matrix[0][1]*Matrix[2][2] - Matrix[0][2]*Matrix[2][1])/det; \
								   temp_matrix[1][1] =   (Matrix[0][0]*Matrix[2][2] - Matrix[0][2]*Matrix[2][0])/det; \
								   temp_matrix[2][1] = - (Matrix[0][0]*Matrix[2][1] - Matrix[0][1]*Matrix[2][0])/det; \
								   temp_matrix[0][2] =   (Matrix[0][1]*Matrix[1][2] - Matrix[0][2]*Matrix[1][1])/det; \
								   temp_matrix[1][2] = - (Matrix[0][0]*Matrix[1][2] - Matrix[0][2]*Matrix[1][0])/det; \
								   temp_matrix[2][2] =   (Matrix[0][0]*Matrix[1][1] - Matrix[0][1]*Matrix[1][0])/det; \
								   EQUATE_MATRIX33(Matrix,temp_matrix);                                               \
							   } while (0)
							   // �������� do-while ��������, ����� �������� ���������
							   // "����� � �������" � �������� ����������, ���������� "else".
							   // ����������� ���� ���
/*3.9. �������� ������ 
res_matrix = matrix_1 + matrix_2 */
#define	ADD_MATRIX33(res_matrix,matrix_1,matrix_2) do {                                                            \
													   unsigned short int i,j;                                     \
													   for(i = 0; i < 3; i++)                                      \
														   for(j = 0; j < 3; j++)                                  \
															   res_matrix[i][j] = matrix_1[i][j] + matrix_2[i][j]; \
												   } while(0) // �������� do-while ��������, ����� �������� ���������
															  // "����� � �������" � �������� ����������, ���������� "else".
															  // ����������� ���� ���
/*3.10. ��������� ������ 
res_matrix = matrix_1 - matrix_2 */
#define	SUB_MATRIX33(res_matrix,matrix_1,matrix_2) do {                                                            \
													   unsigned short int i,j;                                     \
													   for(i = 0; i < 3; i++)                                      \
														   for(j = 0; j < 3; j++)                                  \
															   res_matrix[i][j] = matrix_1[i][j] - matrix_2[i][j]; \
												   } while(0) // �������� do-while ��������, ����� �������� ���������
															  // "����� � �������" � �������� ����������, ���������� "else".
															  // ����������� ���� ���
/*3.11. ��������� ������ 
res_matrix = matrix_1 * matrix_2 */
#define	MULT_MATRIX33(res_matrix,matrix_1,matrix_2) do {                                                                  \
														unsigned short int i,j,k;                                         \
														SET_ZERO_MATRIX33(res_matrix);                                    \
														for(i = 0; i < 3; i++)                                            \
															for(j = 0; j < 3; j++)                                        \
																for(k = 0; k < 3; k++)                                    \
																	res_matrix[i][j] += matrix_1[i][k] * matrix_2[k][j];  \
													} while(0) // �������� do-while ��������, ����� �������� ���������
															   // "����� � �������" � �������� ����������, ���������� "else".
															   // ����������� ���� ���
/*3.12. ��������� �������-������ �� �������. ��������� - ������-������ 
res_vector = vector*matrix */
#define	MULT_VECTOR3_MATRIX33(res_vector,vector,matrix) (res_vector.X = vector.X*matrix[0][0] + vector.Y*matrix[1][0] + vector.Z*matrix[2][0], \
														 res_vector.Y = vector.X*matrix[0][1] + vector.Y*matrix[1][1] + vector.Z*matrix[2][1], \
														 res_vector.Z = vector.X*matrix[0][2] + vector.Y*matrix[1][2] + vector.Z*matrix[2][2])
/*3.13. ��������� ������� �� ������-�������. ��������� - ������-�������
res_vector = matrix*vector */
#define	MULT_MATRIX33_VECTOR3(res_vector,matrix,vector) (res_vector.X = vector.X*matrix[0][0] + vector.Y*matrix[0][1] + vector.Z*matrix[0][2], \
														 res_vector.Y = vector.X*matrix[1][0] + vector.Y*matrix[1][1] + vector.Z*matrix[1][2], \
														 res_vector.Z = vector.X*matrix[2][0] + vector.Y*matrix[2][1] + vector.Z*matrix[2][2])
/*3.14. ��������� ����� �� ������� 
res_matrix = number*matrix */
#define MULT_NUM_MATRIX33(res_matrix,number,matrix) do {                                                  \
														unsigned short int i,j;                           \
														for(i = 0; i < 3; i++)                            \
															for(j = 0; j < 3; j++)                        \
																res_matrix[i][j] = matrix[i][j]*(number); \
													} while(0) // �������� do-while ��������, ����� �������� ���������
															   // "����� � �������" � �������� ����������, ���������� "else".
															   // ����������� ���� ���
/*3.15. ��������� ����� �� ������ 
res_vector = vector*number */
#define MULT_NUM_VECTOR3(res_vector,number,vector) (res_vector.X = vector.X*(number), \
													res_vector.Y = vector.Y*(number), \
													res_vector.Z = vector.Z*(number))
/*3.16. ��������� ��������� (������������ ��� [vector_1,vector_1] ��� ��� vector_1 x vector_1) */
#define MULTVECT_VECTOR3(res_vector,vector_1,vector_2) (res_vector.X = vector_1.Y*vector_2.Z - vector_1.Z*vector_2.Y, \
														res_vector.Y = vector_1.Z*vector_2.X - vector_1.X*vector_2.Z, \
														res_vector.Z = vector_1.X*vector_2.Y - vector_1.Y*vector_2.X)
/*3.17. �������������� �������� ������� � �����*/
#define ADD_NUM_VECTOR3(res_vector, number, vector) (res_vector.X = vector.X + (number), \
													 res_vector.Y = vector.Y + (number), \
													 res_vector.Z = vector.Z + (number))
/*3.18. ������������� ��������*/
#define EQUATE_VECTOR3(dest_vector, source_vector) dest_vector = source_vector

/*������� �����-������ �������� � ������*/
/*4.1. ������ ������� �����������*/
#define PRINT_CONSOLE_VECTOR3(Vector) printf("\n%f\n%f\n%f\n", Vector.X, Vector.Y, Vector.Z)                              // ������ �� �������
#define PRINT_FILE_VECTOR3(Vector, file_descriptor) fprintf(file_descriptor, "%f\n%f\n%f\n", Vector.X, Vector.Y, Vector.Z) // ������ � ���� � ������������
																													       // file_descriptor
/*4.2. ������ ������� ����������� � ���������������*/
// ������ �� �������
#define PRINT_CONSOLE_MATRIX33(matrix) do {                                     \
										   unsigned short int i,j;              \
										   for(i = 0; i < 3; i++){              \
											   if(i) printf("\n");              \
											   for(j = 0; j < 3; j++)           \
												   printf("%f\t",matrix[i][j]); \
										   }                                    \
									   } while(0) // �������� do-while ��������, ����� �������� ���������
											      // "����� � �������" � �������� ����������, ���������� "else".
											      // ����������� ���� ���
// ������ � ���� � ������������ file_descriptor
#define PRINT_FILE_MATRIX33(matrix,file_descriptor) do {                                                      \
														unsigned short int i,l;                               \
														for(i = 0; i < 3; i++) {                              \
															if(i) fprintf(file_descriptor,"\n");              \
															for(j = 0; j < 3; j++)                            \
																fprintf(file_descriptor,"%f\t",matrix[i][j]); \
														}                                                     \
													} while(0) // �������� do-while ��������, ����� �������� ���������
															   // "����� � �������" � �������� ����������, ���������� "else".
															   // ����������� ���� ���

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ���������� � ������� �������                                                                        //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*1. ������������� �������� � ��������� �� �������������*/
/*1.1. ������������ �������� ��������� ������� - ������������� �������*/
Vector3 set_vector3(double X_value, double Y_value, double Z_value);
/* 1.2. ���������� �������� Value �������� Index = (_X,_Y,_Z) == (0,1,2) ������� Vector */
void set_element_vector3(Vector3 *Vector, double Value, unsigned short int Index);
/*1.3.��������� ����� �������*/
double get_length_vector3(Vector3 Vector);
/*1.4. ����������� ����������� �����*/
double get_lambda_vector3(Vector3 Vector); // ������ �������
double get_fi_vector3(Vector3 Vector);     // ������ ������
/*1.5. ������������� �������� �������*/
Vector3 set_zero_vector3(void);
/*1.6. ������������� ���������� �������*/
Vector3 set_unitary_vector3(void);
/*1.7. ��������� �������� �������� ������� Vector � �������� Index = (_X,_Y,_Z) == (0,1,2) */
double get_element_vector3(Vector3 Vector, unsigned short int Index);

/*2. ������������� ������*/
/*2.1. ������������ �������� ���� ��������� ������� - ������������� �������
	   Matrix - ������� ���� Matrix33
	   a_ij - (i,j)-�� ������� ������� Matrix*/
void set_matrix33(Matrix33 Matrix, double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);
/*2.2. ������������ �������� value (i,j)-��� �������� ������� Matrix
	   Matrix - ������� ���� Matrix33*/
void set_element_matrix33(Matrix33 Matrix, unsigned short int i, unsigned short int j, double value);
/*2.3. ������� ��������*/
//������� �������� �� ��������� � ������������� �.�.
void set_svsk_measure_matrix33(Matrix33 Matrix, double nu, double mu);
//������� �������� �� ���������� ������ � ��������� �.�.
void set_nearth_stsk_matrix33(Matrix33 Matrix_GS_ST_l, double B42_l, double L42_l, double B42_ST_l, double L42_ST_l, double Ag_l);
//������� �������� �� ��������� � ���������� ������ �.�.(� ����� ������)
void set_svsk_nearth_matrix33(Matrix33 Matrix, double theta, double psi, double gamma);
//������� �������� �� ���������� ������ � ��������� �.�.(� ���������� �������-����������)
void set_nearth_svsk_rg_matrix33(Matrix33 Matrix, double r, double l, double mm, double n);
//������� �������� �� ���� � ����

/*2.4. ������������� ������� �������*/
void set_zero_matrix33(Matrix33 Matrix);
/*2.5. ������������� ��������� �������*/
void set_unitary_matrix33(Matrix33 Matrix);

/*3. �������� � ��������� � ���������*/
/*3.1. c������� ��������
res_vector = vector_1 + vector_2 */
Vector3 add_vector3(Vector3 vector_1, Vector3 vector_2);
/*3.2. ��������� ��������
res_vector = vector_1 - vector_2 */
Vector3 sub_vector3(Vector3 vector_1, Vector3 vector_2);
/*3.3. ������������ ��������� ��������
res_vector = vector_1 * vector_2 */
Vector3 mult_vector3(Vector3 vector_1, Vector3 vector_2);
/*3.4. ������������ ������� ��������
res_vector = vector_1 / vector_2 */
Vector3 div_vector3(Vector3 vector_1, Vector3 vector_2);
/*3.5. ������������� ���� ������
Matrix_dest = Matrix_source */
void equate_matrix33(Matrix33 matrix_dest, Matrix33 matrix_source);
/*3.6. ���������������� �������, ����� ���� �������� ������� �atrix ����� ���������������*/
void transpose_matrix33(Matrix33 Matrix);
/*3.7. ���������� ������������ �������*/
double det_matrix33(Matrix33 Matrix);
/*3.8. ��������� �������, ����� ���� �������� ������� �atrix ����� ����������
������� ���������� 0 ��� �������� ��������� � 1 ��� ����������� �������*/
unsigned short int inverse_matrix33(Matrix33 Matrix);
/*3.9. �������� ������ 
res_matrix = matrix_1 + matrix_2 */
void add_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2);
/*3.10. ��������� ������ 
res_matrix = matrix_1 - matrix_2 */
void sub_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2);
/*3.11. ��������� ������ 
res_matrix = matrix_1 * matrix_2 */
void mult_matrix33(Matrix33 res_matrix, Matrix33 matrix_1, Matrix33 matrix_2);
/*3.12. ��������� �������-������ �� �������. ��������� - ������-������ 
res_vector = vector*matrix */
Vector3 mult_vector3_matrix33(Vector3 vector, Matrix33 matrix);
/*3.13. ��������� ������� �� ������-�������. ��������� - ������-�������
res_vector = matrix*vector */
Vector3 mult_matrix33_vector3(Matrix33 matrix, Vector3 vector);
/*3.14. ��������� ����� �� ������� 
res_matrix = number*matrix */
void mult_num_matrix33(Matrix33 matrix_res, double number, Matrix33 matrix);
/*3.15. ��������� ����� �� ������ 
res_vector = vector*number */
Vector3 mult_num_vector3(double number, Vector3 vector);
/*3.16. ��������� ��������� (������������ ��� [vector_1,vector_1] ��� ��� vector_1 x vector_1) */
Vector3 multvect_vector3(Vector3 vector_1, Vector3 vector_2);
/*3.17. �������������� �������� ������� � �����*/
Vector3 add_num_vector3(double number, Vector3 vector);
void set_GSSK_STSK_matrix33(Matrix33 Matrix, double B0, double L0, double PSI0, double B, double L);
/*������� �����-������ �������� � ������*/
/*4.1. ������ ������� �����������*/
void print_console_vector3(Vector3 Vector);                     // ������ �� ������� (�����)
void print_file_vector3(Vector3 Vector, FILE *file_descriptor); // ������ � ���� � ������������ file_descriptor
/*4.2. ������ ������� ����������� � ���������������*/
void print_console_matrix33(Matrix33 matrix);                     // ������ �� ������� (�����)
void print_file_matrix33(Matrix33 matrix, FILE *file_descriptor); // ������ � ���� � ������������ file_descriptor

/*������ �������*/
/*5.1. ���������� ����� 10 � ����� �������*/
double pow10(int number);


#endif

