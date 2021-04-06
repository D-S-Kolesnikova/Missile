#pragma once 
#include "Project1.h"
#include "dig.h"


void KOEFFICIENT::Set(const char* AEROcc, const char* FileAddress)
{
	//����������������� ������ ���� ����� ��� ������� �������
	//char AEROcc[]={'C','x','_','N','O','M','I','N','A','L'};
	//������� ���������� �������� KaN_SHOW (�� �������� ����, ��� ������ ������� ����������
	//� ������� �������� � ���� ������ �������� ��� �������

	//������� ���������, ���������� ���������� ��� ��������� ����������,
	//������������ �������� ������ ������������ � ���� AEROcc
	_strnset(aN_SHOW, '\0', sizeof(aN_SHOW) / sizeof(char));
	strcpy(aN_TYPE, aN_SHOW);
	strcpy(aN_ELEMENT, aN_SHOW);

	double** Dig;//�������� ������� ������

	//telo funkcii

	//� ���� ����� ��������� ������������ ������������ ������ ������ �� ���������
	//����� ���������������� ������������� � ����������� ������������� ����������
	//aN_SHOW,aN_TYPE,aN_ELEMENT
	//��������� �������� ���������� ����������������� ������������

	//	if (strstr(AEROcc,"Cx_NOMINAL")!=NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cxp_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "CxpI_NOMINAL") != NULL)Send(AEROcc);
	//if (strstr(AEROcc, "Delta_Cx_tr_NOMINAL") != NULL)Send(AEROcc);
	//if (strstr(AEROcc, "Cy_kr_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al00_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al05_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al10_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al15_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al20_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al25_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_ar_Al30_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "f_ar_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cx_don_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cxp_H_NOMINAL") != NULL)Send(AEROcc);
	//	if (strstr(AEROcc,"Cy_Alfa_al00_NOMINAL")!=NULL)Send(AEROcc);
	if (strstr(AEROcc, "Mz_wz_Al00_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Mx_wx_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cy_M05_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cyp_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "CypI_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cd_M05_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cdp_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "CdpI_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cdp_H_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Cyp_H_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "Mz_wz_NOMINAL") != NULL)Send(AEROcc);
	//	if (strstr(AEROcc,"Cy_ar_delta_NOMINAL")!=NULL)Send(AEROcc);
	if (strstr(AEROcc, "Mx_ko_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "dCx_T_NOMINAL") != NULL)Send(AEROcc);
	if (strstr(AEROcc, "dCd_T_NOMINAL") != NULL)Send(AEROcc);

	/*	if (strstr(AEROcc,"Cn_ar_gp90_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp90_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cn_ar_gp45_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp45_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cn_ar_gp00_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gp00_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cn_ar_gm45_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm45_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cn_ar_gm90_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cn_ar_gm90_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Msh_gp90_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp90_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Msh_gp45_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp45_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Msh_gp00_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gp00_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Msh_gm45_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm45_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Msh_gm90_a00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a05_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a10_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a15_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a25_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Msh_gm90_a30_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cyp_with_OD_M20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Czp_with_OD_M20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mzp_with_OD_M20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Myp_with_OD_M20_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mxp_with_OD_M20_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cyp_with_OD_M08_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cyp_with_OD_M17_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cyp_with_OD_M23_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Czp_with_OD_M08_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Czp_with_OD_M17_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Czp_with_OD_M23_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Mzp_with_OD_M08_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mzp_with_OD_M17_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mzp_with_OD_M23_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Myp_with_OD_M08_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Myp_with_OD_M17_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Myp_with_OD_M23_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Mxp_with_OD_M08_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mxp_with_OD_M17_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mxp_with_OD_M23_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cy_int_x1_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cy_int_x1_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cy_int_x2_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cy_int_x2_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Cz_int_x1_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cz_int_x1_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cz_int_x2_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Cz_int_x2_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Mz_int_x1_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mz_int_x1_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mz_int_x2_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mz_int_x2_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"My_int_x1_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"My_int_x1_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"My_int_x2_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"My_int_x2_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);

		if (strstr(AEROcc,"Mx_int_x1_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mx_int_x1_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mx_int_x2_BettaSm02_NOMINAL")!=NULL)Send(AEROcc);
		if (strstr(AEROcc,"Mx_int_x2_BettaSm00_NOMINAL")!=NULL)Send(AEROcc);    */

		/*	if (strstr(AEROcc,"Fit_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Fik_NOMINAL")!=NULL)Send(AEROcc);
		//	if (strstr(AEROcc,"Kfi_NOMINAL")!=NULL)Send(AEROcc);

			if (strstr(AEROcc,"WW_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Temp_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Ro_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"P8_NOMINAL")!=NULL)Send(AEROcc);

		/*	if (strstr(AEROcc,"Cy_int_M08_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cy_int_M08_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cy_int_M17_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cy_int_M17_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cy_int_M23_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cy_int_M23_BettaS00_NOMINAL")!=NULL)Send(AEROcc);

			if (strstr(AEROcc,"Cz_int_M08_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cz_int_M08_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cz_int_M17_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cz_int_M17_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cz_int_M23_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Cz_int_M23_BettaS00_NOMINAL")!=NULL)Send(AEROcc);

			if (strstr(AEROcc,"Mx_int_M08_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mx_int_M08_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mx_int_M17_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mx_int_M17_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mx_int_M23_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mx_int_M23_BettaS00_NOMINAL")!=NULL)Send(AEROcc);

			if (strstr(AEROcc,"My_int_M08_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"My_int_M08_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"My_int_M17_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"My_int_M17_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"My_int_M23_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"My_int_M23_BettaS00_NOMINAL")!=NULL)Send(AEROcc);

			if (strstr(AEROcc,"Mz_int_M08_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mz_int_M08_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mz_int_M17_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mz_int_M17_BettaS00_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mz_int_M23_BettaSm2_NOMINAL")!=NULL)Send(AEROcc);
			if (strstr(AEROcc,"Mz_int_M23_BettaS00_NOMINAL")!=NULL)Send(AEROcc);  */

//���� ������
	FILE* fError;
	fError = fopen("fError.err", "w");
			//�������� ������������� ���� ������������
	FILE* AERODATE;
	if ((AERODATE = fopen(FileAddress, "r")) == NULL) return; //imya faila zadaetsya vo vhodnih parametrah funkcii
	//���� ������� ����������������� �����
	char cc[256]; char cc_temp[256];
	int num = sizeof(cc) / sizeof(char);
	char chetchik = 0;//C������, ������� �� ��� ��� ��������� ������������? ���� �������� � ����� ����������
	//�������� ������ ���������, � ���������� ��� ���� �����.
	while (fgets(cc, num, AERODATE) != NULL)
	{
		bool read_cc = true;//���������� �� ���������� ������������, ���� � ������ ���� ���������
		//������� ������� '/' � ���� "//" � ����� ������, �� ������ �� �������������
		for (int oo = 0; oo < (num - 1); oo++)
			if ((cc[oo] == '/') && (cc[oo + 1] == '/'))
			{
				read_cc = false;
				break; 
			}
		if (!read_cc) 
			continue;

		if (strstr(cc, aN_SHOW) != NULL)
		{
			++chetchik;
			//for (int oo = 0; oo < (num - 1); oo++)
			//cc_temp[oo] = cc[oo];
			SHOW = (int)StrFileToDouble(cc, fError);
			//printf("cc=%s  aN_SHOW=%s\n",cc,aN_SHOW);
			//printf("chetchik=%d  %d\n",chetchik,SHOW);
		}
		if (strstr(cc, aN_TYPE) != NULL)
		{
			++chetchik;
			TYPE = StrFileToDouble(cc, fError);
			//printf("chetchik=%d  %d\n",chetchik,TYPE);
		}
		if (strstr(cc, aN_ELEMENT) != NULL)
		{
			++chetchik;
			ELEMENT = StrFileToDouble(cc, fError);
			//printf("chetchik=%d  %d\n",chetchik,ELEMENT);
		}
	}

	fseek(AERODATE, 0, SEEK_SET);//����������� ��������� � ������ �����

	//������ �� �������
	AeroK = 0;//���-�� ��������� �� �����������
	AeroI = 0;//���-�� ��������� �� ���������

	//������� �������� ����, ���� ����� ��������� ��������������� �������, ��� ��
	//�� ������� ����-���� �� ������ � ���� �� ����� � ���������
	FILE* putevik;
	putevik = fopen("aero.tmp", "w+");
	//printf("AeroK=%d  AeroI=%d\n",AeroK,AeroI);

//����� ����������� ���� �������� ���� � ���� ���������������� ����������� ����� ������ ����������� � ���������������
//�������� OnlyDig ��� ��������
	if (chetchik == 3)//��� � ��� ������� ���������� 
		while (fgets(cc, num, AERODATE) != NULL)
		{
			bool read_cc = true;//���������� �� ���������� ������������, ���� � ������ ���� ���������
		//������� ������� '/' � ���� "//" � ����� ������, �� ������ �� �������������
			for (int oo = 0; oo < (num - 1); oo++)
				if ((cc[oo] == '/') && (cc[oo + 1] == '/')) { read_cc = false; break; }
			if (!read_cc) continue;

			if (strstr(cc, aN_SHOW) != NULL)
			{
				bool StartReadDig = false;
				while (fgets(cc, num, AERODATE) != NULL)
				{
					if (strstr(cc, "}") != NULL) break;
					if ((strstr(cc, "{") == NULL) && (!StartReadDig)) continue;
					else StartReadDig = true;

					if ((strstr(cc, "=") != NULL) || (strstr(cc, "0") != NULL) ||
						(strstr(cc, "1") != NULL) || (strstr(cc, "2") != NULL) ||
						(strstr(cc, "3") != NULL) || (strstr(cc, "4") != NULL) ||
						(strstr(cc, "5") != NULL) || (strstr(cc, "6") != NULL) ||
						(strstr(cc, "7") != NULL) || (strstr(cc, "8") != NULL) ||
						(strstr(cc, "9") != NULL) || (strstr(cc, ".") != NULL) ||
						(strstr(cc, ",") != NULL) || (strstr(cc, "-") != NULL) ||
						(strstr(cc, "+") != NULL))//������������� ������ � ������� ���������� ������ ��� �������
					{
						if (SHOW == 1) OnlyDig = StrFileToDouble(cc, fError);//���� ���������������� ����������� ����� ����� ������ ������� ��� ��������
						if (SHOW == 0)//� ��� ���� �� ����� ��������, �� ������ �� � �������� ����
						{
							//������ �������
							if (AeroK == 0)//��������� ���������� �������� � ����� �������

							{
								for (int o = 0; (cc[o] != '\n') && (o < 256); o++)
									if ((cc[o] != '0') && (cc[o] != '1') && (cc[o] != '2') && (cc[o] != '3') && (cc[o] != '4') && (cc[o] != '5') && (cc[o] != '6') &&
										(cc[o] != '7') && (cc[o] != '8') && (cc[o] != '9') && (cc[o] != '.') && (cc[o] != ',') && (cc[o] != '-') && (cc[o] != '+') &&
										(cc[o] != ' ') && (cc[o] != '\n')) {
										printf("cc[0]=%c\n", cc[o]);
									}

								//printf("cc=%s\n",cc);
								//printf("stop\n"); 

								for (int o = 0; cc[o] != '\n'; o++)
								{
									if (((cc[o] != ' ') && (cc[o + 1] == ' ')) || ((cc[o] != ' ') && (cc[o + 1] == '\n'))) ++AeroK;
								}
							}
							char sss[256]; _strnset(sss, '\0', 256);//������� ��������� ������ ��� ���������� ��������� �����
							bool ooo = false;//���� � ���, ��� ���������� ������ �����, �����, �������, ����� ��� ���� ��� ������

							for (int o = 0, oo = 0; (cc[o] != '\n') && (o < 256); o++)
							{

								//�� ������ ������ ��������� �� ��������-���������, ���� ����� � ��� �����
								if (cc[o] != ' ')
									if ((cc[o] == '0') || (cc[o] == '1') || (cc[o] == '2') || (cc[o] == '3') || (cc[o] == '4') || (cc[o] == '5') || (cc[o] == '6') ||
										(cc[o] == '7') || (cc[o] == '8') || (cc[o] == '9') || (cc[o] == '.') || (cc[o] == ',') || (cc[o] == '-') || (cc[o] == '+'))
									{
										sss[oo] = cc[o]; ++oo; ooo = true;
									}
								if (o != 0) if ((cc[o] == ' ') && ooo) { sss[oo] = ' '; ++oo; ooo = false; }
							}

							++AeroI;//��������� ���������� ����� � ����� �������

							//printf("AeroK=%d  AeroI=%d\n",AeroK,AeroI);

							fprintf(putevik, "%s\n\0", sss);//������� � �������� ����
							//printf("%s\n%s\n",cc,sss);
							//printf("AeroK=%d AeroI=%d\n",AeroK,AeroI);
						}

					}
					//                          else {printf("ccc=%s\n",cc);}//��� �������� ������ ��� ������������ ������� ����������������� �������� ���������

				}
				break;
			}

		}

	//������ ������������� �������
	if (SHOW == 0)
	{
		//������ ���������� ������������ ������ Dig  c ������� �� �������, ��� ���� ������, ���
		//1�� ������� ���� ����� � ��������, � 1�� ������ ����� ����, ������ �� �� ����������� �����
		//������� �������� ������� ������, �� ��������� ����� ��� ����� 0
		fseek(putevik, 0, SEEK_SET);//����������� ��������� � ������ �����
		Dig = (double**)malloc(AeroK * sizeof(double*));
		for (int o = 0; o < AeroK; o++)
			Dig[o] = (double*)calloc(AeroI, sizeof(double));
		//���������������� ������������� Dig -� ������

		//->AeroK
		//|
		//V
		//AeroI

		int Aerok = 0, Aeroi = 0;

		while (fgets(cc, num, putevik) != NULL)
		{
			char sss[256]; _strnset(sss, '\0', 256);
			Aerok = 0;//C ����� ������ ������� ������;
				//	printf("-|%s\n",cc);
				//	printf("-|%d %d \n",Aerok,Aeroi);
			for (int o = 0, oo = 0; Aerok < AeroK; o++)
			{
				if ((cc[o] == ' ') || (cc[o] == '\n'))
				{
					Dig[Aerok][Aeroi] = StrFileToDouble2(sss, fError);
					//		printf("%s %f\n",sss, Dig[Aerok][Aeroi]);
					_strnset(sss, '\0', 256);
					oo = 0;
					++Aerok;
				}
				else
				{
					sss[oo] = cc[o]; oo++;
				}
			}
			++Aeroi;
			if (Aeroi >= AeroI) break;
		}

		//�������������� Dig -� � ��������� �����������

		alfa = (double*)malloc((AeroI - 1) * sizeof(double));
		mah = (double*)malloc((AeroK - 1) * sizeof(double));
		TblDig = (double**)malloc((AeroK - 1) * sizeof(double*));
		for (int o = 0; o < (AeroK - 1); o++)
			TblDig[o] = (double*)calloc((AeroI - 1), sizeof(double));

		for (int oo = 0; oo < AeroI; oo++)
		{
			for (int o = 0; o < AeroK; o++)
			{
				if ((o > 0) && (oo == 0)) mah[o - 1] = Dig[o][oo];
				else
					if ((oo > 0) && (o == 0)) alfa[oo - 1] = Dig[o][oo];
					else
						if ((o != 0) && (oo != 0)) TblDig[o - 1][oo - 1] = Dig[o][oo];
				/*printf("%f ",Dig[o][oo]);*/
			}
			//printf("\n");
		}


		/*
		for(int oo=0;oo<(AeroI-1);oo++)
		{
		for(int o =0;o<(AeroK-1);o++)
		printf("%f ",TblDig[o][oo]);
		printf("\n");
		}
		*/
		/*
		for(int oo=0;oo<(AeroI-1);oo++)
		printf("%f ",alfa[oo]);
		printf("\n");
		for(int oo=0;oo<(AeroK-1);oo++)
		printf("%f ",mah[oo]);
		printf("\n");
		*/
		for (int o = 0; o < AeroK; o++)
			free(Dig[o]); free(Dig);
	}

	fclose(putevik);
	fclose(AERODATE);
	fclose(fError);
	//printf("AeroK- %f \n",(double)AeroK);	
}

KOEFFICIENT::~KOEFFICIENT()
{
	if (SHOW == 0)
	{
		for (int o = 0; o < (AeroK - 1); o++)
			free(TblDig[o]);
		free(TblDig);
		free(alfa);
		free(mah);
	}
}
//---------------------------------------------------------------------------
//�������� ������������ �� ���� ������
double alur(double x1, double y1, double x2, double y2, double x)
{
	double b, k;
	//////////////////////
	if (x2 == x1) x2 = x1 + 0.0001;
	/////////////////////
	b = (y2 - y1) / (x2 - x1);
	k = (y1 * x2 - y2 * x1) / (x2 - x1);

	return b * x + k;
}
//

//������� �������� ������o�����
double KOEFFICIENT::Get(double m, double a)
{
	if (SHOW == 1) return OnlyDig;
	int k, i;
	double c1, c2;

	for (k = 1; k < (AeroK - 1); k++)
		if (m <= mah[k]) break;
	//printf("k- %d \n",k);	

	if (m < mah[0]) { k = 1; m = mah[0]; } //k=1 ���, ��� �� ������ ������ ������ ������� ���������, ������ 
	//������� �� 1 ������, � ����� ��� ����� ����� ��������������� �������� ������� ������� ������� �����[0], 
	//� �����, �� �� �����, ������ ������� �����
	if (m > mah[AeroK - 2]) { k = AeroK - 2; m = mah[AeroK - 2]; }
	//AeroK - ���������� �������� ������� ���, � ��� - ���������� �������� �����-1, �������������
	//��� ���������� ��������� �������, ������� � 0 ������������ ������� ������ ��� ����� ������ �����-2

	//AeroI-1 - ���-�� ����� � ������� ������
	for (i = 1; i < (AeroI - 1); i++)
		if (a <= alfa[i]) break;
	//printf("i- %d \n",i);	
	//���������� � ����������
	if (a < alfa[0]) { i = 1; a = alfa[0]; }
	if (a > alfa[AeroI - 2]) { i = AeroI - 2; a = alfa[AeroI - 2]; }

	//   mah[k-1]                mah[k]
	//   ........             ............
	//TblDig[k-1][i-1]       TblDig[k][i-1]
	//TblDig[k-1][i]         TblDig[k][i]

	//double alur(double x1, double y1, double x2, double y2, double x)
	//printf("k%d i%d m%f mah[0]%f \n",k,i,m, mah[0]);	
	//printf("mah[k-1]%f,TblDig[k-1][i-1]%f,mah[k]%f,TblDig[k][i-1]%f,m%f\n",mah[k-1],TblDig[k-1][i-1],mah[k],TblDig[k][i-1],m);	
	c1 = alur(mah[k - 1], TblDig[k - 1][i - 1], mah[k], TblDig[k][i - 1], m);
	//printf("c1 %f \n",c1);	

	c2 = alur(mah[k - 1], TblDig[k - 1][i], mah[k], TblDig[k][i], m);
	//printf("c2 %f \n",c2);	

	//alfa[i-1] c1
	//alfa[i]   c2

	return alur(alfa[i - 1], c1, alfa[i], c2, a);
}

void KOEFFICIENT::Send(const char* titul)
{
	//printf("titul %s \n",titul);
	strcpy(aN_SHOW, titul);
	strcat(aN_SHOW, "_SHOW");
	//printf("aN_SHOW %s \n",aN_SHOW);

	strcpy(aN_TYPE, titul);
	strcat(aN_TYPE, "_TYPE");
	//printf("aN_TYPE %s \n",aN_TYPE);

	strcpy(aN_ELEMENT, titul);
	strcat(aN_ELEMENT, "_ELEMENT");
	//printf("aN_ELEMENT %s \n",aN_ELEMENT);
}


/*void main(void)
{
	KOEFFICIENT Cx_NOMINAL("Cx_NOMINAL");
   printf("%f \n",Cx_NOMINAL.alfa[4]);
   printf("%f \n",Cx_NOMINAL.TblDig[2][2]);
   printf("cx c3%f \n",Cx_NOMINAL.Get(0.01,0.1));
printf("AeroK%d \n",Cx_NOMINAL.AeroK);
printf("AeroI%d \n",Cx_NOMINAL.AeroI);
printf("SHOW%d \n",Cx_NOMINAL.SHOW);
printf("TYPE%d \n",Cx_NOMINAL.TYPE);
//
	KOEFFICIENT Cy_alfa_NOMINAL("Cy_alfa_NOMINAL");
	printf("SHOW%d \n",Cy_alfa_NOMINAL.SHOW);
	printf("cy c3%f \n",Cy_alfa_NOMINAL.Get(0.01,0.1));


	KOEFFICIENT Mx_delta_NOMINAL("Mx_delta_NOMINAL");
	printf("SHOW%d \n",Mx_delta_NOMINAL.SHOW);
	printf("mx c3%f \n",Mx_delta_NOMINAL.Get(0.01,0.1));


} */