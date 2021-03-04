#pragma once
#include "math_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>

//параметры стандартной атмосферы
struct ParamAtmosferaStruct {
	double h;		//геометрическая высота  (м)
	double H;		//геопотенциальная высота  (гп. м.)
	double M;		//молярная масса
	double T;		//Молекулярная температура (К)
	double P;		//давление в изометрическом слое (Па)
	double Ro;	//плотность в изометрическом слое (кг/м3)
	double Gam;	// удельный вес (Н/м3)
	double n;		// концентрация частиц воздуха (1/м3)
	double Nu;	//кинематическая вязкость (м2/с)
	double Mu;	//динамическая вязкость (Па*с)
	double L;		//средняя длина свободного пробега молекул (м)
	double a;		//скорость распространения звука (м/с)
	double g;		//ускорение силы тяжести (м/с2)
};


double interpol2(double x1, double y1, double x2, double y2, double x) /*интерполирование по двум точкам*/
{
	double y;
	y = ((x - x1) / (x2 - x1) * y2) + ((x - x2) / (x1 - x2) * y1);
	return y;
}

ParamAtmosferaStruct ParamAtmos(double hv) {
	ParamAtmosferaStruct Atm;
	double r;
	double G0;
	double Rz;
	double R;
	double S;
	double Mc;
	double Na;
	double x;
	double sigma;
	double bettaS;
	double Muc;
	double Tc;
	double ac;
	double Hpc;
	double lc;
	double nc;
	double Pc;
	double vc;
	double gammaC;
	double Nuc;
	double LamdaC;
	double wc;
	double Roc;

	double Roz, am;
	double Tz, Hz, Pz, Bet, logP;
	double B0, B1, B2, B3;

	int nm;
	double A0, A1, A2, A3, A4;

	r = 6356767;
	G0 = 9.80665;
	Rz = 8314.32;
	R = 287.05287;
	S = 110.4;
	Mc = 28.964420;
	Na = 602.257E+24;
	x = 1.4;
	sigma = 0.365E-9;
	bettaS = 1.458E-6;
	Muc = 17.894E-6;
	Tc = 288.15;
	ac = 340.294;
	Hpc = 8434.5;
	lc = 66.328E-9;
	nc = 25.471E+24;
	Pc = 101325;
	vc = 458.94;
	gammaC = 12.013;
	Nuc = 14.607E-6;
	LamdaC = 25.343E-3;
	wc = 6.9193E+9;
	Roc = 1.225;



	Atm.h = hv;
	Atm.g = G0 * (r / (r + Atm.h)) * (r / (r + Atm.h));
	Atm.H = (r * Atm.h) / (r + Atm.h);

	if ((Atm.h >= -2000) && (Atm.h <= 94000)) Atm.M = Mc;
	else
	{
		if ((Atm.h > 94000) && (Atm.h <= 97500)) Atm.M = 28.82 + 0.158 * sqrt(1 - 7.5E-8 * (Atm.h - 94000) * (Atm.h - 94000)) - 2.479E-4 * sqrt(97500 - Atm.h);
		else
		{
			if ((Atm.h > 97500) && (Atm.h <= 120000)) { Atm.M = 28.85 - 0.0001511 * (Atm.h - 97500); }
			else
			{
				if ((Atm.h > 120000) && (Atm.h <= 250000)) { B0 = 46.9083;  B1 = -29.71210e-5; B2 = 12.08693e-10; B3 = -1.85675e-15; };
				if ((Atm.h > 250000) && (Atm.h <= 400000)) { B0 = 40.4668;  B1 = -15.52722e-5; B2 = 3.55735e-10; B3 = -3.02340e-16; };
				if ((Atm.h > 400000) && (Atm.h <= 650000)) { B0 = 6.3770;  B1 = 6.25497e-5; B2 = -1.10144e-10; B3 = 3.36907e-17; };
				if ((Atm.h > 650000) && (Atm.h <= 900000)) { B0 = 75.6896;  B1 = -17.61243e-5; B2 = 1.33603e-10; B3 = -2.87884e-17; };
				if ((Atm.h > 900000) && (Atm.h <= 1050000)) { B0 = 112.4838; B1 = -30.68086e-5; B2 = 2.90329e-10; B3 = -9.20616e-17; };
				if ((Atm.h > 1050000) && (Atm.h <= 1200000)) { B0 = 9.8970;  B1 = -1.19732e-5; B2 = 7.78247e-12; B3 = -1.77541e-18; };
				Atm.M = B0 + B1 * Atm.h + B2 * Atm.h * Atm.h + B3 * Atm.h * Atm.h * Atm.h;
			};
		};
	};
	//ShowMessage(String(Atm.M));
	if ((Atm.H >= -2000) && (Atm.H < 0)) { Tz = 301.15;  Hz = -2000;  Bet = -0.0065;  Pz = 127774; };
	if ((Atm.H >= 0) && (Atm.H < 11000)) { Tz = 288.15;  Hz = 0;  Bet = -0.0065;  Pz = 101325; };
	if ((Atm.H >= 11000) && (Atm.H < 20000)) { Tz = 216.65;  Hz = 11000;  Bet = 0.0000;   Pz = 22632; };
	if ((Atm.H >= 20000) && (Atm.H < 32000)) { Tz = 216.65;  Hz = 20000;  Bet = 0.001;    Pz = 5474.87; };
	if ((Atm.H >= 32000) && (Atm.H < 47000)) { Tz = 228.65;  Hz = 32000;  Bet = 0.0028;   Pz = 868.014; };
	if ((Atm.H >= 47000) && (Atm.H < 51000)) { Tz = 270.65;  Hz = 47000;  Bet = 0.0000;   Pz = 110.906; };
	if ((Atm.H >= 51000) && (Atm.H < 71000)) { Tz = 270.65;  Hz = 51000;  Bet = -0.0028;  Pz = 66.9384; };
	if ((Atm.H >= 71000) && (Atm.H < 85000)) { Tz = 214.65;  Hz = 71000;  Bet = -0.002;   Pz = 3.95639; };
	if ((Atm.H >= 85000) && (Atm.H < 94000)) { Tz = 186.65;  Hz = 85000;  Bet = 0.0000;   Pz = 0.363702; };
	if ((Atm.H >= 94000) && (Atm.H < 102450)) { Tz = 186.525; Hz = 94000;  Bet = 0.00225;    Pz = 0.083634; };
	if ((Atm.H >= 102450) && (Atm.H < 117777)) { Tz = 203.81;  Hz = 102450; Bet = 0.0095;    Pz = 0.016439; };
	if ((Atm.h >= 120000) && (Atm.h < 140000)) { Tz = 334.42;  Hz = 120000; Bet = 0.011259; Pz = 0.00266618; };
	if ((Atm.h >= 140000) && (Atm.h < 160000)) { Tz = 559.60;  Hz = 140000; Bet = 0.0068;   Pz = 0.000826375; };
	if ((Atm.h >= 160000) && (Atm.h < 200000)) { Tz = 695.60;  Hz = 160000; Bet = 0.00397;  Pz = 0.00030362; };
	if ((Atm.h >= 200000) && (Atm.h < 250000)) { Tz = 834.40;  Hz = 200000; Bet = 0.00175;  Pz = 0.0000853026; };
	if ((Atm.h >= 250000) && (Atm.h < 325000)) { Tz = 941.90;  Hz = 250000; Bet = 0.00057;  Pz = 0.0000247564; };
	if ((Atm.h >= 325000) && (Atm.h < 400000)) { Tz = 984.65;  Hz = 325000; Bet = 0.00015;  Pz = 0.00000524521; };
	if ((Atm.h >= 400000) && (Atm.h < 600000)) { Tz = 995.90;  Hz = 400000; Bet = 0.00002;  Pz = 0.00000145265; };
	if ((Atm.h >= 600000) && (Atm.h < 800000)) { Tz = 999.90;  Hz = 600000; Bet = 0.0000005; Pz = 0.0000000821535; };
	if ((Atm.h >= 800000) && (Atm.h < 1200000)) { Tz = 1000.0;  Hz = 800000; Bet = 0.00000;  Pz = 0.0000000170593; };
	if ((Atm.h >= -2000) && (Atm.h <= 120000))   Atm.T = Tz + Bet * (Atm.H - Hz);
	if ((Atm.h >= 120000) && (Atm.h <= 1200000)) Atm.T = Tz + Bet * (Atm.h - Hz);

	if ((Atm.h >= -2000) && (Atm.h <= 120000))
	{
		if (Bet != 0.0)
			logP = log10(Pz) - (G0 / (Bet * R)) * (log10((Bet * (Atm.H - Hz) + Tz) / Tz));
		else
			logP = log10(Pz) - (0.434294 * G0) / (R * Atm.T) * (Atm.H - Hz);
		Atm.P = pow(10.0, logP);
		Atm.n = 7.243611 * Atm.P / Atm.T;
	};

	if ((Atm.h > 120000) && (Atm.h <= 1200000)) {

		if ((Atm.h >= 120000) && (Atm.h < 150000)) { A0 = 0.21000587e+4;  A1 = -0.5618444757e-1;  A2 = 0.5663986231e-6;  A3 = -0.2547466858e-11;  A4 = 0.4309844119e-17; nm = 17; };
		if ((Atm.h >= 150000) && (Atm.h < 200000)) { A0 = 0.10163937e+4;  A1 = -0.2119530830e-1;  A2 = 0.1671627815e-6;  A3 = -0.5894237068e-12;  A4 = 0.7826684089e-18; nm = 16; };
		if ((Atm.h >= 200000) && (Atm.h < 250000)) { A0 = 0.76315750e+3;  A1 = -0.1150600844e-1;  A2 = 0.6612598428e-7;  A3 = -0.1708736137e-12;  A4 = 0.1669823114e-18; nm = 15; };
		if ((Atm.h >= 250000) && (Atm.h < 350000)) { A0 = 0.18822030e+3;  A1 = -0.2265999519e-2;  A2 = 0.1041726141e-7;  A3 = -0.2155574922e-13;  A4 = 0.1687430962e-19; nm = 15; };
		if ((Atm.h >= 350000) && (Atm.h < 450000)) { A0 = 0.28048230e+3;  A1 = -0.2432231125e-2;  A2 = 0.8055024663e-8;  A3 = -0.1202418519e-13;  A4 = 0.6805101379e-20; nm = 14; };
		if ((Atm.h >= 450000) && (Atm.h < 600000)) { A0 = 0.55993620e+3;  A1 = -0.3714141392e-2;  A2 = 0.9358870345e-8;  A3 = -0.1058591881e-13;  A4 = 0.4525531532e-20; nm = 13; };
		if ((Atm.h >= 600000) && (Atm.h < 800000)) { A0 = 0.83587560e+3;  A1 = -0.4265393073e-2;  A2 = 0.8252842085e-8;  A3 = -0.7150127437e-14;  A4 = 0.2335744331e-20; nm = 12; };
		if ((Atm.h >= 800000) && (Atm.h < 1000000)) { A0 = 0.83699600e+2;  A1 = -0.3162492458e-3;  A2 = 0.4602064246e-9;  A3 = -0.3021858469e-15;  A4 = 0.7512304301e-22; nm = 12; };
		if ((Atm.h >= 1000000) && (Atm.h < 1200000)) { A0 = 0.38322000e+2;  A1 = -0.5098000000e-4;  A2 = 0.1810000000e-10; A3 = -0.0000000000e-00;  A4 = 0.0000000000e-00; nm = 11; };

		Atm.n = (A0 + A1 * Atm.h + A2 * Atm.h * Atm.h + A3 * Atm.h * Atm.h * Atm.h + A4 * Atm.h * Atm.h * Atm.h * Atm.h) * pow10(nm);
		Atm.P = (Atm.n * R * Atm.T) / Na;
	};
	Atm.Ro = (Atm.M * Atm.P) / (Rz * Atm.T);

	Atm.Gam = Atm.Ro * Atm.g;
	Atm.L = 2.332376 * pow10(-5) * Atm.T / Atm.P;
	Atm.a = 20.046796 * sqrt(Atm.T);
	Atm.Mu = (bettaS * pow(Atm.T, 1.5)) / (Atm.T + S);
	Atm.Nu = Atm.Mu / Atm.Ro;

	return Atm;
}