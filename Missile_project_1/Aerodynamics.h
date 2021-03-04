#pragma once
#include <stdio.h>
#include <fstream>
#include <xlocale>

class Aerodynamic
{	public:
	Aerodynamic(double Mah, double alfaSpace, double alfa, double betta)
	{
		double Calculation(double Mah, double alfaSpace, double alfa, double betta);
	}

	/*void Calculation(double Mah, double alfaSpace, double alpha, double betta, double *Cx, double *Cy) {
		*Cx = Cx(alfaSpace, Mah);
		*Cy = Cy(alpha, Mah);
	}*/


	double Cx(double alphap, double Mah) {

		double Cxa14_1 = 0, Cxa5_1 = 0;
		const int s = 18, q = 8;
		double Cya[q] = { 0.0391, 0.0345, 0.0306, 0.0276, 0.0254, 0.0236, 0.0224, 0.0216 };
		double alpha[s] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 };

		double Cxa14[s] = { 0.0462, 0.0478, 0.05, 0.0525, 0.0573, 0.062,
							0.0697, 0.077, 0.0867, 0.0956, 0.1069, 0.1183,
							0.1268, 0.1397, 0.1533, 0.1678, 0.1830, 0.1989 };

		double Cxa5[s] = { 0.0310, 0.0320, 0.0330, 0.0350, 0.0380, 0.0410,
							0.0440, 0.0490, 0.0540, 0.0600, 0.0660, 0.0745,
							0.0819, 0.0913, 0.1018, 0.1131, 0.1251, 0.1379 };

		int j = 0;
		for (j; j <= 18; j += 1) {
			if (alphap >= alpha[j] && alphap < alpha[j + 1]) {
				Cxa14_1 = (alphap * (Cxa14[j + 1] - Cxa14[j]) /
					(alpha[j + 1] - alpha[j])) + Cxa14[j] - alpha[j] *
					((Cxa14[j + 1] - Cxa14[j]) / (alpha[j + 1] - alpha[j]));
			}

			if (alphap >= alpha[17]) {
				Cxa14_1 = Cxa14[17];
			}

			if (alphap <= alpha[0]) {
				Cxa14_1 = Cxa14[0];
			}

		}

		j = 0;
		for (j; j <= 18; j += 1) {
			if (alphap >= alpha[j] && alphap < alpha[j + 1]) {
				Cxa5_1 = (alphap * (Cxa5[j + 1] - Cxa5[j]) / (alpha[j + 1] - alpha[j])) + Cxa5[j] - alpha[j] * ((Cxa5[j + 1] - Cxa5[j]) / (alpha[j + 1] - alpha[j]));
			}

			if (alphap >= alpha[17]) {
				Cxa5_1 = Cxa5[17];
			}

			if (alphap <= alpha[0]) {
				Cxa5_1 = Cxa5[0];
			}

		}
		return (Cxa5_1 * (Mah - 1.4) + Cxa14_1 * (5 - Mah)) / (5 - 1.4);
	}

	 double  Cy(double alpha, double Mah) 
	 {
		 const int q = 8;
		 double Mahcy[q] = { 1.4, 2, 2.5, 3, 3.5, 4, 4.5, 5 };
		 double Cya[q] = { 0.0391, 0.0345, 0.0306, 0.0276, 0.0254, 0.0236, 0.0224, 0.0216 };
		int j = 0;
		 for (j; j <= 8; j += 1) {
			 if (Mah >= Mahcy[j] && Mah < Mahcy[j + 1]) {
				 return (Mah * (Cya[j + 1] - Cya[j]) /
					 (Mahcy[j + 1] - Mahcy[j])) + Cya[j] -
					 Mahcy[j] * ((Cya[j + 1] - Cya[j]) /
					 (Mahcy[j + 1] - Mahcy[j]));
			 }

			 if (Mah >= Mahcy[7]) {
				 return  Cya[7];
			 }

			 if (Mah <= Mahcy[0]) {
				 return  Cya[0];
			 }
		 }
	 }

	 double Cz(double betta, double Mah) \
	 {

		 const int q = 8;
		 double Mahcz[q] = { 1.4, 2, 2.5, 3, 3.5, 4, 4.5, 5 };
		 double Cza[q] = { 0.0391, 0.0345, 0.0306, 0.0276, 0.0254, 0.0236, 0.0224, 0.0216 };
		 int j = 0;
		 for (j; j <= 8; j += 1) {
			 if (Mah >= Mahcz[j] && Mah < Mahcz[j + 1]) {
				 return (Mah * (Cza[j + 1] - Cza[j]) /
					 (Mahcz[j + 1] - Mahcz[j])) + Cza[j] -
					 Mahcz[j] * ((Cza[j + 1] - Cza[j]) /
					 (Mahcz[j + 1] - Mahcz[j]));
			 }

			 if (Mah >= Mahcz[7]) {
				 return  Cza[7];
			 }

			 if (Mah <= Mahcz[0]) {
				 return  Cza[0];
			 }
		 }
	 }
};