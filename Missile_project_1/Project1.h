#ifndef _Project1_h_
#define _Project1_h_
//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include <conio.h>
#include "dig.h"


//---------------------------------------------------------------------------
struct KOEFFICIENT//�����, ����������� �������� ���������������� ������������� �� ����� Akfm.dt
	{
	          ~KOEFFICIENT();
	    double Get (double, double);//������� �� (���, �����) ������ �������� ����������������� ������������ ������������ �������� ������ ��������� ������
		void Send (const char *);//������� ��������� ��������� ��������� aN_SHOW[50] aN_TYPE[50] aN_ELEMENT[50];
	    char SHOW;
	    char TYPE;
	    char ELEMENT;

	 int  AeroK;
	 int  AeroI;
	 double OnlyDig;
	 double *alfa; 
	 double *mah;
	 double **TblDig;
	 //������� ���������, ���������� ��������� ��� ��������� ����������,
	 //������������ �������� ������ ������������ � ���� AEROcc
	 char aN_SHOW[50];
	 char aN_TYPE[50];
	 char aN_ELEMENT[50];
     void Set(const char *, const char* FileAddress);//������ �����������

	};
//---------------------------------------------------------------------------
//extern struct KOEFFICIENT Mx_ko_NOMINAL;




#endif