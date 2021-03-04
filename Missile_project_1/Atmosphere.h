#pragma once
class Atm
{
public:

	const double p0N = 101325,
		R_earth = 6371000,
		omega_earth = 7.2921E-5,
		Mu = 3.98600E14;
	const double H_i[9] = { -2,0,11,20,32,47,51,71,85 };
	const double T_i[9] = { 301.15, 288.15, 216.65, 216.65, 288.65, 270.65, 270.65, 214.65, 186.65 };
	const double beta[9] = { -6.5, -6.5, 0, 1, 2.8, 0, -2.8, -2, 0 };; //градиент температуры
	const double p_i[9] = { 127774, 101325, 22632, 5474.87, 868.014, 110.906, 66.9384, 3.95639, 0.44571 };
	const double p_0 = 101325;
	const double T_0 = 288.15;
	const double g_0 = 9.80665;
	const double hi = 1.4;
	const double R = 287.05287;
	const double R_mu = 8314.32;

	

	double p(double h)
	{
		double T;
		double p;
		double H = (R_earth * h) / (1000 * (R_earth + h));
		if (H > 85)

			return(0);
		for (int i = 0; i <= 8; i++)
		{
			if ((H >= H_i[i]) && (H < H_i[i + 1]))
			{
				T = T_i[i] + beta[i] * (H - H_i[i]);

				if ((beta[i]) < 0.0000001)
				{
					return p = p_i[i] * exp(-g_0 * (H - H_i[i]) * 1000 / (R * T_i[i]));
				}
				else return p = p_i[i] * pow((1 + beta[i] * (H - H_i[i]) / T_i[i]), -(g_0 * 1000) / (R * beta[i]));
			}
		}
	}

	double ro(double h)
	{
		double T;
		double p;
		double H = (R_earth * h) / (1000 * (R_earth + h)); //геопотенциальная высота
		if (H > 85)
			return(0);

		for (int i = 0; i <= 8; i++)
		{
			if ((H >= H_i[i]) && (H < H_i[i + 1]))
			{
				T = T_i[i] + beta[i] * (H - H_i[i]);

				if ((beta[i]) < 0.0000001)
				{
					p = p_i[i] * exp(-g_0 * (H - H_i[i]) * 1000 / (R * T_i[i]));
				}
				else p = p_i[i] * pow((1 + beta[i] * (H - H_i[i]) / T_i[i]), -(g_0 * 1000) / (R * beta[i]));
				return(p / (R * T));

			}
		}
	}
	double a(double h)
	{
		double T;
		double p;
		double H = (R_earth * h) / (1000 * (R_earth + h)); //геопотенциальная высота
		if (H > 85)
		{
			T = T_i[8];
			return 20.046796 * sqrt(T);

		}


		for (int i = 0; i <= 8; i++)
		{
			if ((H >= H_i[i]) && (H < H_i[i + 1]))
			{
				T = T_i[i] + beta[i] * (H - H_i[i]);

				if ((beta[i]) < 0.0000001)
				{
					p = p_i[i] * exp(-g_0 * (H - H_i[i]) * 1000 / (R * T_i[i]));
				}
				else p = p_i[i] * pow((1 + beta[i] * (H - H_i[i]) / T_i[i]), -(g_0 * 1000) / (R * beta[i]));
				return 20.046796 * sqrt(T);

			}
		}
	}

	void Calculation(double H, double &P, double &Ro, double &A)
	{
		P = p(H);
		Ro = ro(H);
		A = a(H);
	}

};
