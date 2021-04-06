#pragma once
/* математические макросы */
#define _SQR(w) ((w)*(w)) // квадрат числа

#define	pr1Speed_X Left[0]
#define	pr1Speed_Y Left[1]
#define	pr1Speed_Z Left[2]
#define	pr1Point_X Left[3]
#define	pr1Point_Y Left[4]
#define	pr1Point_Z Left[5]
#define	pr1omegaX   Left[6]
#define	pr1omegaY   Left[7]
#define	pr1omegaZ   Left[8]
#define	pr1ro_rg    Left[9]
#define	pr1lamda_rg Left[10]
#define	pr1mu_rg    Left[11]
#define	pr1nu_rg    Left[12]
#define pr1B		Left[13]
#define pr1L		Left[14]
#define pr1H		Left[15]

#define target_r target[0]
#define target_fi target[1]
#define target_hi target[2]
#define target_pr1r target[3]
#define target_pr1fi target[4]
#define target_pr1hi target[5]
#define target_X target[6]
#define target_Y target[7]
#define target_Z target[8]
#define MISS target[9]



#define	Speed_X el[0]
#define	Speed_Y el[1]
#define	Speed_Z el[2]
#define	Point_X el[3]
#define	Point_Y el[4]
#define	Point_Z el[5]
#define	omegaX   el[6]
#define	omegaY   el[7]
#define	omegaZ   el[8]
#define	ro_rg    el[9]
#define	lamda_rg el[10]
#define	mu_rg    el[11]
#define	nu_rg    el[12]
#define	B_    el[13]
#define	L_    el[14]
#define	H_    el[15]


#define TIME Time[0]

#define ALFA angle[0]
#define BETTA angle[1]
#define PITCH angle[2]
#define YAW angle[3]
#define ROLL angle[4]
#define TETA angle[5]
#define PSI angle[6]
#define L42 angle[7]
#define B42 angle[8]
#define PSI0 angle[9]

#define pr1Pitch Eiler[0]
#define pr1Yaw Eiler[1]
#define pr1Roll Eiler[2]
#define DeltaPitch Eiler[3]
#define DeltaYaw Eiler[4]
#define DeltaRoll Eiler[5]

#define Gx G[0]
#define Gy G[1]
#define Gz G[2]

#define delta_1 delta[0]
#define delta_2 delta[1]
#define delta_3 delta[2]
#define delta_4 delta[3]

/*
#define A11 cos(PITCH)*cos(YAW);
#define A12 sin(PITCH);
#define A13 -cos(PITCH) * sin(YAW);
#define A21 -sin(PITCH)*cos(YAW) *cos(ROLL) + sin(YAW) *sin(ROLL);
#define A22 cos(PITCH) * cos(ROLL);
#define A23 cos(YAW) * sin(ROLL) + sin(PITCH) * sin(YAW) * cos(ROLL);
#define A31 sin(PITCH) * cos(YAW) * sin(ROLL) + sin(YAW) * cos(ROLL) ;
#define A32 -cos(PITCH) * sin(ROLL);
#define A33 cos(YAW) * cos(ROLL) - sin(YAW) * sin(PITCH) * sin(ROLL);
#define pr1A11 omegaZ*A21 - omegaY*A31;
#define pr1A12 omegaZ*A22 - omegaY*A32;
#define pr1A13 omegaZ*A23 - omegaY*A33;
#define pr1A21 -omegaZ*A11 + omegaX*A31;
#define pr1A22 -omegaZ*A12 + omegaX*A32;
#define pr1A23 -omegaZ*A13 + omegaX*A33;
#define pr1A31 omegaY*A11 - omegaX*A21;
#define pr1A32 omegaY*A12 - omegaX*A22;
#define pr1A33 omegaY*A13 - omegaX*A23;
*/
