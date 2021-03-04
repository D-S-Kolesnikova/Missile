#include "EarthModel.h"
#include "math_lib.h"
//��������� ������ �����, ����������� ��� ������� (�������� �������� � ������� ������� ��� ���������� �����������)
static	double	A_GD, B_GD;				//������� (6378245.0 �) � ����� (6356863.0 �) ������� ������� ����������
static	double	E_GD, E2_GD;			//������ (0.08181337) � ������ ��������������� ������� ����������
static	double	EARTH_RATE;				//������� �������� �������� ����� (7.292116E-5 ���/���)
static	double	G_EQUATOR;				//�������� ��������� ���������� ������� �� �������� (9.78049)
static	double	Q_GR;					//��������� ������������ ����, ����������� ���������� �������� �����, � ���� ������� �� ��������
static	double	K1, K2, K3, K4, K5, K6;	//��������������� ������������

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
	//������� ������� �������� ����������
	//H - ������ ��� ������������ ��������-���������� �����������
	//B - ������������� ������
	Earth_Struct OutStruct;
	double  W_gd, V_gd;				//�������� �������������� �������
	double  R_phi, R_lambda;		//������� ������� �������� ����������� ���������� ��������
	double  R_gd0;					//��������������� ������-������ �� ����������� ����� � �����
	double  phi_gr0, phi_gr;		//��������������� ������
	Vector F1_GCS;					//�������� ������� ������������� ���� ��������� �� �����
	Vector F1_GS;					//�������� ������� ������������� ���� ��������� �� ����
	Vector AcclCp_GS;				//�������� ������� �������������������� ��������� �� ����
	Vector F2_GS;					//�������� ������� ������������� ���� ���� ������� �� ����

	double R_gd;					//������-������ ��������� ������ ���� �� � �����

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

	/*-------�������, ����� ������� � ��������� ������ �����-------//
	//  �������� ���� ��������� �� "�����"
		F1_GCS.X = lnau.lna.gwz*sqrt(1 - _SQR(lnau.lna.slac));
		F1_GCS.Y = naout.gr + lnau.lna.gwz*lnau.lna.slac;
		F1_GCS.Z = 0.;
	//-------------------------------------------------------------*/

	F1_GS.X = F1_GCS.X * cos(B - phi_gr) - F1_GCS.Y * sin(B - phi_gr);
	F1_GS.Y = F1_GCS.X * sin(B - phi_gr) + F1_GCS.Y * cos(B - phi_gr);
	F1_GS.Z = 0.;

	AcclCp_GS.X = _SQR(EARTH_RATE) * R_gd * (sin(B) * cos(B) * cos(B - phi_gr) + _SQR(sin(B)) * sin(B - phi_gr));
	AcclCp_GS.Y = _SQR(EARTH_RATE) * R_gd * (-_SQR(cos(B)) * cos(B - phi_gr) - sin(B) * cos(B) * sin(B - phi_gr));
	AcclCp_GS.Z = 0.;

	F2_GS = F1_GS - AcclCp_GS;
	//F2_GS = sub_vector3(F1_GS, AcclCp_GS);

	OutStruct.A_gd = A_GD;
	OutStruct.G_eq = G_EQUATOR;
	OutStruct.Rate = EARTH_RATE;
	OutStruct.E_gd = E_GD;
	OutStruct.R_phi = R_phi;
	OutStruct.R_lambda = R_lambda;
	OutStruct.R_gd = R_gd;
	OutStruct.phi_gr = phi_gr;
	OutStruct.AcclCp_GS = F1_GS;	 //�������, ����� ������� � ��������� ������ �����
	OutStruct.GField = F2_GS;

	return OutStruct;
}


/*Earth_Struct	GetEarthParameters(double H, double B) {
	//������� ������� �������� ����������
	//H - ������ ��� ������������ ��������-���������� �����������
	//B - ������������� ������
	Earth_Struct OutStruct;
	double  e_gd, e2_gd;  			//������ � ������ ��������������� ������� ����������
	double  W_gd, V_gd;				//�������� �������������� �������
	double  R_phi, R_lambda;		//������� ������� �������� ����������� ���������� ��������
	double  R_gd0;								//��������������� ������-������ �� ����������� ����� � �����
	double  phi_gr0, phi_gr;					//��������������� ������
	double  q_gr;					//��������� ������������ ����, ����������� ���������� �������� �����, � ���� ������� �� ��������
	Vector3 F1_GCS;					//�������� ������� ������������� ���� ��������� �� �����
	Vector3 F1_GS;					//�������� ������� ������������� ���� ��������� �� ����
	Vector3 AcclCp_GS;				//�������� ������� �������������������� ��������� �� ����
	Vector3 F2_GS;					//�������� ������� ������������� ���� ���� ������� �� ����
	double K1, K2, K3, K4, K5, K6;	//��������������� ������������
	double R_gd;					//������-������ ��������� ������ ���� �� � �����
	double u_gd;
	double l_gr;
	double l1_gr;
	double Dmu_gr;
	double a1, a2, a3;
	double P0_gr;

	e_gd = sqrt(_SQR(A_GD) - _SQR(B_GD)) / A_GD; //����������� � ������ �������������
	e2_gd = sqrt(_SQR(A_GD) - _SQR(B_GD)) / B_GD; //����������� � ������ �������������
	W_gd = sqrt(1 - _SQR(e_gd * sin(B)));
	V_gd = sqrt(1 + _SQR(e_gd * sin(B)));
	R_phi = A_GD * (1 - _SQR(e_gd)) / (_SQR(W_gd) * W_gd);
	R_lambda = A_GD / W_gd;
	R_gd0 = sqrt(_SQR(R_lambda * cos(B)) + _SQR(R_lambda * (1 - _SQR(e_gd)) * sin(B)));
	phi_gr0 = atan(_SQR(B_GD) / _SQR(A_GD) * tan(B));
	R_gd = sqrt(_SQR(R_gd0 * cos(phi_gr0) + H * cos(B)) + _SQR(R_gd0 * sin(phi_gr0) + H * sin(B)));
	u_gd = atan(B_GD / A_GD * tan(B));
	phi_gr = atan(B_GD / A_GD * tan(u_gd));
	l_gr = sqrt(_SQR(E_GD) / (1 - _SQR(E_GD)));
	a1 = 2 * (_SQR(A_GD) - _SQR(B_GD));
	a2 = _SQR(R_gd) - (_SQR(A_GD) - _SQR(B_GD));
	a3 = sqrt(pow(R_gd, 4) + _SQR(_SQR(A_GD) - _SQR(B_GD)) - 2 * R_gd * R_gd * (_SQR(A_GD) - _SQR(B_GD)) * cos(2 * phi_gr));
	l1_gr = sqrt((2 * (_SQR(A_GD) - _SQR(B_GD))) / (_SQR(R_gd) - (_SQR(A_GD) - _SQR(B_GD)) + sqrt(pow(R_gd, 4) + _SQR(_SQR(A_GD) - _SQR(B_GD)) - 2 * R_gd * R_gd * (_SQR(A_GD) - _SQR(B_GD)) * cos(2 * phi_gr))));
	Dmu_gr = _SQR(EARTH_RATE) * pow(l_gr, 3)/(2 * M_PI * ((3+l_gr*l_gr)*atan(l_gr)-3*l_gr));
	P0_gr = 2 * M_PI * Dmu_gr * (1 + l_gr * l_gr) / pow(l_gr, 3) * (atan(l_gr) - l_gr / (1 + l_gr * l_gr));
	//K1 = G_EQUATOR * 0.5 * (q_gr - _SQR(e_gd)) * (1.0 + _SQR(e_gd) * (7.0 * _SQR(e_gd) - 30.0 * q_gr) / (14.0 * (q_gr - _SQR(e_gd)))); //����������� � ������ �������������
	//K2 = K1 * _SQR(e_gd) * (30.0 * q_gr - 21.0 * _SQR(e_gd)) / (14.0 * (q_gr - _SQR(e_gd))); //����������� � ������ �������������
	//K3 = K1 * _SQR(e_gd) * (7.0 * _SQR(e_gd) - 10.0 * q_gr) / (2.0 * (q_gr - _SQR(e_gd))); //����������� � ������ �������������
	//K4 = -G_EQUATOR * (1.0 - 0.5 * _SQR(e_gd) * (1.0 + 0.25 * _SQR(e_gd)) + 1.5 * q_gr * (1.0 - 5.0 / 14.0 * _SQR(e_gd))); //����������� � ������ �������������
	//K5 = -G_EQUATOR * (_SQR(e_gd) * (1.0 - 0.5 * _SQR(e_gd)) - q_gr * (1.0 - 15.0 / 7.0 * _SQR(e_gd))) / 2.0; //����������� � ������ �������������
	//K6 = G_EQUATOR * (_SQR(e_gd) - q_gr - _SQR(e_gd) * (0.5 * _SQR(e_gd) - 15.0 / 7.0 * q_gr)) * 1.5; //����������� � ������ �������������

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