#ifndef _dig_h_
#define _dig_h_

#include <stdio.h>
#include <math.h>
#include <cstring>

double StrFileToDouble(char str[256], FILE*fg)
		{
		fprintf(fg,"str!!!%s\n",str);
		if ((str[0]=='/')||(str[0]=='<')||(str[0]=='>')) return 0.0;
		//Исследуем в переданной строке все символы после знака "="
		int i=0,j=0,k=0;

		for (;(str[i]!='=');i++)
		if ((str[i]=='\n')||(str[i]=='\0')) return 0.0;

			for (;(str[j]!='\n')&&(str[j]!='\0');j++);
			char *str1,*str2;

			str1= new char[j-i];
			str2= new char[j-i];

			for(k=0;(k+i+1)<=j;k++)	str1[k]=str2[k]='\0';
			for(k=0;(k+i+1)<=j;k++)	str1[k]=str[k+i+1];
	
		fprintf(fg,"str1???%s\n",str1);
		//Определяем символы, которые могут быть в числе
		char dig[]  = {'0','1','2','3','4','5','6','7','8','9','-','+','.',','};

		int kk=0;
		//Считаем сколько в исследуемом элементе символов, невходящих в
		//вышеобозначенный массив (это kk)
		for (i=0;i<k;i++)
		{
			for(j=0;j<(sizeof(dig)/sizeof(char));j++)
				if (str1[i]==dig[j]) {str2[kk]=str1[i]; kk++; break;}
		}
		
	fprintf(fg,"str2???%s\n",str2);


	
	//множитель - определяет знак числа 
	int Znak;
	//классификатор, отвечает за наличие знака - или + перед числом
	int BeZnak=0;
		double Cel=0;
		if (str2[0]=='-') {Znak=-1; BeZnak=1;} else Znak=1;
		if (str2[0]=='+') BeZnak=1;
		//Померяем целочисленную часть
		for(i=BeZnak;(str2[i]!='.')&&(i<kk);i++);
		for(j=BeZnak;j<i;j++)
			{
				switch(str2[j])
				{
			case'0':  break;
				case'1': Cel+=1*pow(10.0,(double)i-j-1); break;
					case'2': Cel+=2*pow(10.0,(double)i-j-1); break;
						case'3': Cel+=3*pow(10.0,(double)i-j-1); break;
							case'4': Cel+=4*pow(10.0,(double)i-j-1); break;
								case'5': Cel+=5*pow(10.0,(double)i-j-1); break;
									case'6': Cel+=6*pow(10.0,(double)i-j-1); break;
										case'7': Cel+=7*pow(10.0,(double)i-j-1); break;
											case'8': Cel+=8*pow(10.0,(double)i-j-1); break;
												case'9': Cel+=9*pow(10.0,(double)i-j-1); break;
				}
			}
				//Analiziruem desyatichnuyu chast chisla
			//esli konechno je ona est
	if(i!=kk)
	{
		for(j=i+1;j<kk;j++)
		{
			switch(str2[j])
			{
			case'0':  break;
				case'1': {Cel+=1*pow(10.0,(double)i-j); break;}
					case'2': {Cel+=2*pow(10.0,(double)i-j); break;}
						case'3': {Cel+=3*pow(10.0,(double)i-j); break;}
							case'4': {Cel+=4*pow(10.0,(double)i-j); break;}
								case'5': {Cel+=5*pow(10.0,(double)i-j); break;}
									case'6': {Cel+=6*pow(10.0,(double)i-j); break;}
										case'7': {Cel+=7*pow(10.0,(double)i-j); break;}
											case'8': {Cel+=8*pow(10.0,(double)i-j); break;}
												case'9': {Cel+=9*pow(10.0,(double)i-j); break;}
			}    			
		}
	}
			
    			delete [] str1;
			delete [] str2;
			Cel*=Znak;
	fprintf(fg,"str3???%f\n", Cel/*strtod(str2,e)*/);
			return Cel;
		}


double StrFileToDoubleWF(char str[256])
		{
		if ((str[0]=='/')||(str[0]=='<')||(str[0]=='>')) return 0.0;
		//Исследуем в переданной строке все символы после знака "="
		int i=0,j=0,k=0;

		for (;(str[i]!='=');i++)
		if ((str[i]=='\n')||(str[i]=='\0')) return 0.0;

			for (;(str[j]!='\n')&&(str[j]!='\0');j++);
			char *str1,*str2;

			str1= new char[j-i];
			str2= new char[j-i];

			for(k=0;(k+i+1)<=j;k++)	str1[k]=str2[k]='\0';
			for(k=0;(k+i+1)<=j;k++)	str1[k]=str[k+i+1];

		//Определяем символы, которые могут быть в числе
		char dig[]  = {'0','1','2','3','4','5','6','7','8','9','-','+','.',','};

		int kk=0;
		//Считаем сколько в исследуемом элементе символов, невходящих в
		//вышеобозначенный массив (это kk)
		for (i=0;i<k;i++)
		{
			for(j=0;j<(sizeof(dig)/sizeof(char));j++)
				if (str1[i]==dig[j]) {str2[kk]=str1[i]; kk++; break;}
		}





	//множитель - определяет знак числа
	int Znak;
	//классификатор, отвечает за наличие знака - или + перед числом
	int BeZnak=0;
		double Cel=0;
		if (str2[0]=='-') {Znak=-1; BeZnak=1;} else Znak=1;
		if (str2[0]=='+') BeZnak=1;
		//Померяем целочисленную часть
		for(i=BeZnak;(str2[i]!='.')&&(i<kk);i++);
		for(j=BeZnak;j<i;j++)
			{
				switch(str2[j])
				{
			case'0':  break;
				case'1': Cel+=1*pow(10.0,(double)i-j-1); break;
					case'2': Cel+=2*pow(10.0,(double)i-j-1); break;
						case'3': Cel+=3*pow(10.0,(double)i-j-1); break;
							case'4': Cel+=4*pow(10.0,(double)i-j-1); break;
								case'5': Cel+=5*pow(10.0,(double)i-j-1); break;
									case'6': Cel+=6*pow(10.0,(double)i-j-1); break;
										case'7': Cel+=7*pow(10.0,(double)i-j-1); break;
											case'8': Cel+=8*pow(10.0,(double)i-j-1); break;
												case'9': Cel+=9*pow(10.0,(double)i-j-1); break;
				}
			}
				//Analiziruem desyatichnuyu chast chisla
			//esli konechno je ona est
	if(i!=kk)
	{
		for(j=i+1;j<kk;j++)
		{
			switch(str2[j])
			{
			case'0':  break;
				case'1': {Cel+=1*pow(10.0,(double)i-j); break;}
					case'2': {Cel+=2*pow(10.0,(double)i-j); break;}
						case'3': {Cel+=3*pow(10.0,(double)i-j); break;}
							case'4': {Cel+=4*pow(10.0,(double)i-j); break;}
								case'5': {Cel+=5*pow(10.0,(double)i-j); break;}
									case'6': {Cel+=6*pow(10.0,(double)i-j); break;}
										case'7': {Cel+=7*pow(10.0,(double)i-j); break;}
											case'8': {Cel+=8*pow(10.0,(double)i-j); break;}
												case'9': {Cel+=9*pow(10.0,(double)i-j); break;}
			}
		}
	}

    			delete [] str1;
			delete [] str2;
			Cel*=Znak;
			return Cel;
		}

/////////////////////////////////////////////
//Эта функция позволяет избежать проблеммы если отсутствует знак "="
double StrFileToDouble2(char str[256], FILE*fg)
{
  char strr[257];
  _strnset(strr,'\0',257);
  strcat(strr,"=");
  strcat(strr,str);
  return StrFileToDouble(strr, fg);
}

double StrFileToDouble2WF(char str[256])
{
  char strr[257];
  _strnset(strr,'\0',257);
  strcat(strr,"=");
  strcat(strr,str);
  return StrFileToDoubleWF(strr);
}

/////////////////////////////////////////////


#endif
