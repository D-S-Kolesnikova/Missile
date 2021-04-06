#include "AHK.h"

// ADH_OUT_PUT BB;
//-----------------------------------------------------------------------------------
double AeroDH::alur(double x1, double y1, double x2, double y2, double x)
{
	double b, k;
	//////////////////////
	if (x2 == x1) x2 = x1 + 0.0001;
	/////////////////////
	b = (y2 - y1) / (x2 - x1);
	k = (y1 * x2 - y2 * x1) / (x2 - x1);
	return b * x + k;
};
//-----------------------------------------------------------------------------------
double AeroDH::Cyp_(double M, double AlfaSpace, double H, double P)
{
	double Cyp;
	if (H < 80)
	{
		if ((M <= 0.5) && (P > 0)) { Cyp = Cy_M05_NOMINAL.Get(1, AlfaSpace); }
		else
		{
			if (AlfaSpace <= 30.0) { Cyp = Cyp_NOMINAL.Get(M, AlfaSpace); }
			else
			{
				Cyp = CypI_NOMINAL.Get(AlfaSpace, M);
			}
		}
	}
	double Cyp1, Cyp2;
	if ((H >= 80) && (H < 100))
	{
		if (M <= 0.5) { Cyp1 = Cy_M05_NOMINAL.Get(1, AlfaSpace); }
		else
		{
			if (AlfaSpace <= 30.0) { Cyp1 = Cyp_NOMINAL.Get(M, AlfaSpace); }
			else
			{
				Cyp1 = CypI_NOMINAL.Get(AlfaSpace, M);
			}
		}
		Cyp2 = Cyp_H_NOMINAL.Get(H, AlfaSpace);
		Cyp = alur(80, Cyp1, 100, Cyp2, H);
	}
	if (H >= 100) { Cyp = Cyp_H_NOMINAL.Get(H, AlfaSpace); }
	return Cyp;
};
//-----------------------------------------------------------------------------------
double AeroDH::Cdp_(double M, double AlfaSpace, double H, double P, double coFi, char B_T)
{
	double Cdp;
	double dCdp;
	if (H < 80)
	{
		if ((M <= 0.5) && (P > 0)) { Cdp = Cd_M05_NOMINAL.Get(1, AlfaSpace); }
		else
		{
			if (AlfaSpace <= 30.0) { Cdp = Cdp_NOMINAL.Get(M, AlfaSpace); }
			else
			{
				Cdp = CdpI_NOMINAL.Get(AlfaSpace, M);
			}
		}
	}
	double Cdp1, Cdp2;
	if ((H >= 80) && (H < 100))
	{
		if (M <= 0.5) { Cdp1 = Cd_M05_NOMINAL.Get(1, AlfaSpace); }
		else
		{
			if (AlfaSpace <= 30.0) { Cdp1 = Cdp_NOMINAL.Get(M, AlfaSpace); }
			else
			{
				Cdp1 = CdpI_NOMINAL.Get(AlfaSpace, M);
			}
		}
		Cdp2 = Cdp_H_NOMINAL.Get(H, AlfaSpace);
		Cdp = alur(80, Cdp1, 100, Cdp2, H);
	}
	if (H >= 100) { Cdp = Cdp_H_NOMINAL.Get(H, AlfaSpace); }
	if (B_T == 0) { dCdp = dCd_T_NOMINAL.Get(1, M) * coFi * coFi; }
	else
	{
		dCdp = 0.0;
	}
	return Cdp + dCdp;
};
//-----------------------------------------------------------------------------------
double AeroDH::Cx_(double M, double AlfaSpace, double H, double non, double coAp, char B_T)
{
	double Cxp, dCxp;
	if (H < 90)
	{
		if (AlfaSpace <= 30.0) { Cxp = Cxp_NOMINAL.Get(M, AlfaSpace); }
		else
		{
			Cxp = CxpI_NOMINAL.Get(AlfaSpace, M);
		}
	}
	double Cxp1, Cxp2;
	if ((H >= 90) && (H < 100))
	{
		if (AlfaSpace <= 30.0) { Cxp1 = Cxp_NOMINAL.Get(M, AlfaSpace); }
		else
		{
			Cxp1 = CxpI_NOMINAL.Get(AlfaSpace, M);
		}
		Cxp2 = Cxp_H_NOMINAL.Get(H, AlfaSpace);
		Cxp = alur(90, Cxp1, 100, Cxp2, H);
	}
	if (H >= 100) { Cxp = Cxp_H_NOMINAL.Get(H, AlfaSpace); }
	if (B_T == 0) { dCxp = dCx_T_NOMINAL.Get(1, M) * coAp; }
	else
	{
		dCxp = 0.0;
	}
	//Delta_Cx_tr_NOMINAL.Get(H, M) * coAp
	return Cxp + Cx_don_NOMINAL.Get(M, non) - Cx_don_NOMINAL.Get(M, 0) + dCxp;
};
//-----------------------------------------------------------------------------------
double AeroDH::Cy_ar_(double M, double Alfa, double Delta)
{
	double Cy_ar;
	int sigA = 1;
	if (Alfa < 0.0) { sigA = -1; Delta = -Delta; }
	Alfa = fabs(Alfa);
	double Cy_ar1, Cy_ar2;
	if ((Alfa >= 0.0) && (Alfa <= 5.0))
	{
		Cy_ar1 = Cy_ar_Al00_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al05_NOMINAL.Get(M, Delta);
		Cy_ar = alur(0.0, Cy_ar1, 5.0, Cy_ar2, Alfa);
	}
	if ((Alfa > 5.0) && (Alfa <= 10.0))
	{
		Cy_ar1 = Cy_ar_Al05_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al10_NOMINAL.Get(M, Delta);
		Cy_ar = alur(5.0, Cy_ar1, 10.0, Cy_ar2, Alfa);
	}
	if ((Alfa > 10.0) && (Alfa <= 15.0))
	{
		Cy_ar1 = Cy_ar_Al10_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al15_NOMINAL.Get(M, Delta);
		Cy_ar = alur(10.0, Cy_ar1, 15.0, Cy_ar2, Alfa);
	}
	if ((Alfa > 15.0) && (Alfa <= 20.0))
	{
		Cy_ar1 = Cy_ar_Al15_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al20_NOMINAL.Get(M, Delta);
		Cy_ar = alur(15.0, Cy_ar1, 20.0, Cy_ar2, Alfa);
	}
	if ((Alfa > 20.0) && (Alfa <= 25.0))
	{
		Cy_ar1 = Cy_ar_Al20_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al25_NOMINAL.Get(M, Delta);
		Cy_ar = alur(20.0, Cy_ar1, 25.0, Cy_ar2, Alfa);
	}
	if ((Alfa > 25.0) && (Alfa <= 30.0))
	{
		Cy_ar1 = Cy_ar_Al25_NOMINAL.Get(M, Delta);
		Cy_ar2 = Cy_ar_Al30_NOMINAL.Get(M, Delta);
		Cy_ar = alur(25.0, Cy_ar1, 30.0, Cy_ar2, Alfa);
	}
	return sigA * Cy_ar;
};
//-----------------------------------------------------------------------------------
double	AeroDH::Km_(double M, double Alfa, double Delta)
{
	double K_alfa, Km;
	if (M >= 3) { Km = 1.0; }
	if (M < 3)
	{
		if ((Alfa * Delta <= 0) && (Alfa >= 0)) { K_alfa = -0.0266; };
		if ((Alfa * Delta <= 0) && (Alfa < 0)) { K_alfa = 0.0266; };
		if (Alfa * Delta > 0) { K_alfa = 0; };
		Km = (0.2 - K_alfa * Alfa) * (M - 0.6) * (M - 0.6) / 5.76 + K_alfa * Alfa + 0.8;
	}
	return Km;
};
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
void AeroDH::SSet(const char* FileAddress)
{
	Cxp_NOMINAL.Set("Cxp_NOMINAL", FileAddress);
	CxpI_NOMINAL.Set("CxpI_NOMINAL", FileAddress);
	Delta_Cx_tr_NOMINAL.Set("Delta_Cx_tr_NOMINAL", FileAddress);
	Cy_kr_NOMINAL.Set("Cy_kr_NOMINAL", FileAddress);

	Cy_ar_Al00_NOMINAL.Set("Cy_ar_Al00_NOMINAL", FileAddress);
	Cy_ar_Al05_NOMINAL.Set("Cy_ar_Al05_NOMINAL", FileAddress);
	Cy_ar_Al10_NOMINAL.Set("Cy_ar_Al10_NOMINAL", FileAddress);
	Cy_ar_Al15_NOMINAL.Set("Cy_ar_Al15_NOMINAL", FileAddress);
	Cy_ar_Al20_NOMINAL.Set("Cy_ar_Al20_NOMINAL", FileAddress);
	Cy_ar_Al25_NOMINAL.Set("Cy_ar_Al25_NOMINAL", FileAddress);
	Cy_ar_Al30_NOMINAL.Set("Cy_ar_Al30_NOMINAL", FileAddress);

	f_ar_NOMINAL.Set("f_ar_NOMINAL", FileAddress);
	Cx_don_NOMINAL.Set("Cx_don_NOMINAL", FileAddress);
	Cxp_H_NOMINAL.Set("Cxp_H_NOMINAL", FileAddress);
	//		 Cy_Alfa_al00_NOMINAL.Set("Cy_Alfa_al00_NOMINAL");
	Mz_wz_Al00_NOMINAL.Set("Mz_wz_Al00_NOMINAL", FileAddress);
	Mx_wx_NOMINAL.Set("Mx_wx_NOMINAL", FileAddress);
	Cy_M05_NOMINAL.Set("Cy_M05_NOMINAL", FileAddress);
	Cyp_NOMINAL.Set("Cyp_NOMINAL", FileAddress);
	CypI_NOMINAL.Set("CypI_NOMINAL", FileAddress);
	Cd_M05_NOMINAL.Set("Cd_M05_NOMINAL", FileAddress);
	Cdp_NOMINAL.Set("Cdp_NOMINAL", FileAddress);
	CdpI_NOMINAL.Set("CdpI_NOMINAL", FileAddress);
	Cyp_H_NOMINAL.Set("Cyp_H_NOMINAL", FileAddress);
	Cdp_H_NOMINAL.Set("Cdp_H_NOMINAL", FileAddress);
	Mz_wz_NOMINAL.Set("Mz_wz_NOMINAL", FileAddress);
	//		 Cy_ar_delta_NOMINAL.Set("Cy_ar_delta_NOMINAL");
	Mx_ko_NOMINAL.Set("Mx_ko_NOMINAL", FileAddress);
	dCx_T_NOMINAL.Set("dCx_T_NOMINAL", FileAddress);
	dCd_T_NOMINAL.Set("dCd_T_NOMINAL", FileAddress);

};
//-----------------------------------------------------------------------------------

DESIGNER AeroDH::KitADH(ADH_IN_PUT AA)
{
	DESIGNER FF;  //структура выходных данных

	double coAp = cos(AA.Alfa__) * cos(AA.Betta__); //косинус пространственного угла атаки
	if (coAp >= 1) coAp = 1;
	if (coAp <= -1) coAp = -1;
	double AlfaSpace = acos(coAp) * ToDegreee;     //пространственный угол атаки в град
	AA.H__ = AA.H__ / 1000;
	double Fi_p_sin__, Fi_p_cos__;             //синус и косинус пространственного угла крена
	if ((AA.Alfa__ == 0) && (AA.Betta__ == 0)) { Fi_p_cos__ = 1;  Fi_p_sin__ = 0; }
	else
	{
		Fi_p_sin__ = (sin(AA.Betta__) / (sqrt(sin(AA.Alfa__) * sin(AA.Alfa__)
			* cos(AA.Betta__) * cos(AA.Betta__) + sin(AA.Betta__) * sin(AA.Betta__))));
		Fi_p_cos__ = ((sin(AA.Alfa__) * cos(AA.Betta__)) / (sqrt(sin(AA.Alfa__) * sin(AA.Alfa__)
			* cos(AA.Betta__) * cos(AA.Betta__) + sin(AA.Betta__) * sin(AA.Betta__))));
	}
	FF.AlfaSpace = AlfaSpace;
	FF.Fi_p = atan2(Fi_p_sin__, Fi_p_cos__) * ToDegreee;
	//printf("cosFi=%e	sinFi=%e	Alfa p=%e\n", Fi_p_cos__, Fi_p_sin__, AlfaSpace);
	//момент крена от корпуса

	//FF.Mxk_ = Mx_ko_NOMINAL.Get(AA.M__, AlfaSpace) * (4 * Fi_p_sin__ * Fi_p_cos__ * Fi_p_cos__ * Fi_p_cos__ - 4 * Fi_p_sin__
		//* Fi_p_sin__ * Fi_p_sin__ * Fi_p_cos__) * AA.S_har__ * AA.l_har__ * AA.q__;
	FF.Mxk_ = 0;
	//демпфирующий момент по крену
	FF.Mx_wx_ = Mx_wx_NOMINAL.Get(1, AA.M__) * AA.q__ * AA.l_har__ * AA.l_har__ * AA.S_har__ * AA.wx__ / AA.V__;
	//определение коэффициента демпфирующего момента по крену и тангажу
	double mz_wz;
	if (AA.M__ < 2.0) { mz_wz = Mz_wz_Al00_NOMINAL.Get(AA.x_ct__, AA.M__); }
	else
	{
		mz_wz = Mz_wz_NOMINAL.Get(AlfaSpace, AA.M__);
	}
	//демпфирующий момент по рысканью
	FF.My_wy_ = mz_wz * AA.q__ * AA.L_har__ * AA.L_har__ * AA.S_har__ * AA.wy__ / AA.V__;
	//демпфирующий момент по тангажу
	FF.Mz_wz_ = mz_wz * AA.q__ * AA.L_har__ * AA.L_har__ * AA.S_har__ * AA.wz__ / AA.V__;
	//коэффициент а/д нормальной силы в системе координат, связанной с пространственным углом атаки
	double Cyp = Cyp_(AA.M__, AlfaSpace, AA.H__, AA.Px__);
	//подьемная и боковая сила от корпуса ЛА
	FF.Fyk_ = Cyp * Fi_p_cos__ * AA.q__ * AA.S_har__;
	FF.Fzk_ = -Cyp * Fi_p_sin__ * AA.q__ * AA.S_har__;
	//коэффициент центра давления
	double Cdp = Cdp_(AA.M__, AlfaSpace, AA.H__, AA.Px__, Fi_p_cos__, AA.B_T);
	//момент тангажа и рысканья
	FF.Myk_ = Cyp * Fi_p_sin__ * AA.q__ * AA.S_har__ * (AA.x_ct__ - Cdp * AA.L_har__);
	FF.Mzk_ = Cyp * Fi_p_cos__ * AA.q__ * AA.S_har__ * (AA.x_ct__ - Cdp * AA.L_har__);
	//сила лобового сопротивления корпуса
	FF.Fxk_ = Cx_(AA.M__, AlfaSpace, AA.H__, AA.n__, coAp, AA.B_T) * AA.q__ * AA.S_har__;

	if (AA.GSS__ == 0) {
		double Alfa1 = atan(sqrt(2.0) * (tan(AA.Alfa__) - tan(AA.Betta__) / cos(AA.Alfa__)) / 2.0) * ToDegreee;
		double Alfa2 = atan(sqrt(2.0) * (tan(AA.Alfa__) + tan(AA.Betta__) / cos(AA.Alfa__)) / 2.0) * ToDegreee;
		//printf("al1 = %e	al2 = %e\n", Alfa1, Alfa2);
		double Delta1 = AA.Delta1__ * ToDegreee;
		double Delta2 = AA.Delta2__ * ToDegreee;
		double Delta3 = AA.Delta3__ * ToDegreee;
		double Delta4 = AA.Delta4__ * ToDegreee;
		double sgA1 = 1, sgA2 = 1;
		if (Alfa1 < 0.0) { sgA1 = -1; };
		if (Alfa2 < 0.0) { sgA2 = -1; };

		double Cy_kr_1 = sgA1 * Cy_kr_NOMINAL.Get(AA.M__, fabs(Alfa1));
		double Cy_kr_2 = sgA2 * Cy_kr_NOMINAL.Get(AA.M__, fabs(Alfa2));

		FF.Fx_con_ =
			0.5 * Km_(AA.M__, Alfa2, AA.Delta1__) * f_ar_NOMINAL.Get(AA.M__, Alfa1) * tan(AA.Delta1__) * Cy_kr_2
			+ 0.5 * Km_(AA.M__, Alfa2, AA.Delta1__) * f_ar_NOMINAL.Get(AA.M__, Alfa1) * tan(AA.Delta1__) * Cy_ar_(AA.M__, Alfa2, Delta1)
			+ 0.5 * Km_(AA.M__, Alfa1, AA.Delta2__) * f_ar_NOMINAL.Get(AA.M__, -Alfa2) * tan(AA.Delta2__) * Cy_kr_1
			+ 0.5 * Km_(AA.M__, Alfa1, AA.Delta2__) * f_ar_NOMINAL.Get(AA.M__, -Alfa2) * tan(AA.Delta2__) * Cy_ar_(AA.M__, Alfa1, Delta2)
			+ 0.5 * Km_(AA.M__, Alfa2, AA.Delta3__) * f_ar_NOMINAL.Get(AA.M__, -Alfa1) * tan(AA.Delta3__) * Cy_kr_2
			+ 0.5 * Km_(AA.M__, Alfa2, AA.Delta3__) * f_ar_NOMINAL.Get(AA.M__, -Alfa1) * tan(AA.Delta3__) * Cy_ar_(AA.M__, Alfa2, Delta3)
			+ 0.5 * Km_(AA.M__, Alfa1, AA.Delta4__) * f_ar_NOMINAL.Get(AA.M__, Alfa2) * tan(AA.Delta4__) * Cy_kr_1
			+ 0.5 * Km_(AA.M__, Alfa1, AA.Delta4__) * f_ar_NOMINAL.Get(AA.M__, Alfa2) * tan(AA.Delta4__) * Cy_ar_(AA.M__, Alfa1, Delta4);

		FF.Fx_con_ *= AA.q__ * AA.S_har__;


		FF.Fy_con_ = (0.5 * Cy_ar_(AA.M__, Alfa2, Delta1) * f_ar_NOMINAL.Get(AA.M__, Alfa1)
			+ 0.5 * Cy_ar_(AA.M__, Alfa1, Delta2) * f_ar_NOMINAL.Get(AA.M__, -Alfa2)
			+ 0.5 * Cy_ar_(AA.M__, Alfa2, Delta3) * f_ar_NOMINAL.Get(AA.M__, -Alfa1)
			+ 0.5 * Cy_ar_(AA.M__, Alfa1, Delta4) * f_ar_NOMINAL.Get(AA.M__, Alfa2)) * sqrt(2.0) / 2.0 * AA.q__ * AA.S_har__;

		FF.Fz_con_ = (-0.5 * Cy_ar_(AA.M__, Alfa2, Delta1) * f_ar_NOMINAL.Get(AA.M__, Alfa1)
			+ 0.5 * Cy_ar_(AA.M__, Alfa1, Delta2) * f_ar_NOMINAL.Get(AA.M__, -Alfa2)
			- 0.5 * Cy_ar_(AA.M__, Alfa2, Delta3) * f_ar_NOMINAL.Get(AA.M__, -Alfa1)
			+ 0.5 * Cy_ar_(AA.M__, Alfa1, Delta4) * f_ar_NOMINAL.Get(AA.M__, Alfa2)) * sqrt(2.0) / 2.0 * AA.q__ * AA.S_har__;

		FF.Mx_con_ = (0.5 * Cy_ar_(AA.M__, Alfa2, Delta1) * f_ar_NOMINAL.Get(AA.M__, Alfa1)
			+ 0.5 * Cy_ar_(AA.M__, Alfa1, Delta2) * f_ar_NOMINAL.Get(AA.M__, -Alfa2)
			- 0.5 * Cy_ar_(AA.M__, Alfa2, Delta3) * f_ar_NOMINAL.Get(AA.M__, -Alfa1)
			- 0.5 * Cy_ar_(AA.M__, Alfa1, Delta4) * f_ar_NOMINAL.Get(AA.M__, Alfa2)) * AA.q__ * AA.S_har__ * AA.l2_ar__;

		FF.My_con_ = FF.Fz_con_ * (AA.l1_ar__ - AA.x_ct__);
		FF.Mz_con_ = FF.Fy_con_ * (-AA.l1_ar__ + AA.x_ct__);
	}
	else
	{
		FF.Fx_con_ = 0.0;
		FF.Fy_con_ = 0.0;
		FF.Fz_con_ = 0.0;
		FF.Mx_con_ = 0.0;
		FF.My_con_ = 0.0;
		FF.Mz_con_ = 0.0;
	}
	/*printf("Fxk = %e	Fyk = %e	Fzk = %e\n", FF.Fxk_, FF.Fyk_, FF.Fzk_);
	printf("Fxy = %e	Fyy = %e	Fzy = %e\n", FF.Fx_con_, FF.Fy_con_, FF.Fz_con_);
	printf("Mxk = %e	Myk = %e	Mzk = %e\n", FF.Mxk_, FF.Myk_, FF.Mzk_);
	printf("Mxy = %e	Myy = %e	Mzy = %e\n", FF.Mx_con_ ,FF.My_con_, FF.Mz_con_);
	printf("Mxwx = %e	Mywy = %e	Mzwx = %e\n", FF.Mx_wx_, FF.My_wy_, FF.Mz_wz_);*/
	//getch();
	return FF;
};