#pragma once
typedef struct Initial_Conditions
{
	double X0;
	double Y0;
	double Z0;
	double X_target;
	double Y_target;
	double Z_target;
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
	double x_cm;
	double l1_ar;
	double l2_ar;
	double l_har;
	double L_har;
	
} Object;



typedef struct _ModelParams_Struct
{
	double Step;    		//��� ������������� �����
	double PrintStep;		//��� ������ � ����
	
} ModelParams_Struct;

typedef	struct _Matrix_Struct
{
	Matrix GS_SV;		//������� ������������ ��������� ����� ���� � ����
	Matrix GS_ST;		//������� ������������ ��������� ����� ���� � ����
	Matrix ST_GS;		//������� ������������ ��������� ����� ���� � ����
	Matrix ST_SV;		//������� ������������ ��������� ����� ���� � ����

	Matrix RateSV_I;	//���oc������������� ������� ������� ��������� ����
	Matrix RateGS_I;	//���oc������������� ������� ������� ��������� ����
	Matrix RateST_I;  //������������������ ������� ������� ��������� ����

	unsigned char N;	//����� ����������� ���� double ��� �������� ������� ���������
						//Vector3 ������� �� 3-� double
						//Matrix33 ������� �� 9-�� double
} Matrix_Struct;