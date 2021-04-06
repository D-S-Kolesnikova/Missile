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
struct KOEFFICIENT//Класс, позволяющий получать аэродинамические коэффиициенты из файла Akfm.dt
	{
	          ~KOEFFICIENT();
	    double Get (double, double);//функция по (мах, атака) выдает значения аэродинамического коэффициента относительно которого создан экземпляр класса
		void Send (const char *);//Функция формирует строковые параметры aN_SHOW[50] aN_TYPE[50] aN_ELEMENT[50];
	    char SHOW;
	    char TYPE;
	    char ELEMENT;

	 int  AeroK;
	 int  AeroI;
	 double OnlyDig;
	 double *alfa; 
	 double *mah;
	 double **TblDig;
	 //Массивы элементов, содержащих текстовое имя требуемых параметров,
	 //определяемых вводимым именем коэффициента в виде AEROcc
	 char aN_SHOW[50];
	 char aN_TYPE[50];
	 char aN_ELEMENT[50];
     void Set(const char *, const char* FileAddress);//Бывший конструктор

	};
//---------------------------------------------------------------------------
//extern struct KOEFFICIENT Mx_ko_NOMINAL;




#endif