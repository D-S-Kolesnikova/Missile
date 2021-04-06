#pragma once

#include <conio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "structs.h"
using namespace std;

void Print_()
{
	ofstream file;
	if (TIME < 1E-10)
	{
		
		file.open("file_result.txt");

		file << setw(20) << " Time,ñ" << " "
			<< setw(20) << "õ,m" << " "
			<< setw(20) << "y,m" << " "
			<< setw(20) << "z,m" << " "
			<< setw(20) << "Speed_X,m/s" << " "
			<< setw(20) << "Speed_Y,m/s" << " "
			<< setw(20) << "Speed_Z,m/s" << " "
			<< setw(20) << "Omegax,rad/c" << " "
			<< setw(20) << "Omegay,rad/c" << " "
			<< setw(20) << "Omegaz,rad/c" << " "
			<< setw(20) << "Alfa,degree" << " "
			<< setw(20) << "Betta,degree" << " "
			<< setw(20) << "Pitch,degree" << " "
			<< setw(20) << "Yaw,degree" << " "
			<< setw(20) << "Roll,degree" << " "
			<< setw(20) << "TETA,degree" << " "
			<< setw(20) << "PSI,degree" << " "
			<< setw(20) << "B,degree" << " "
			<< setw(20) << "L,degree" << " "
			<< setw(20) << "H,m" << " "
			<< setw(20) << "deltaPitch" << " "
			<< setw(20) << "deltaYaw" << " "
			<< setw(20) << "deltaRoll" << " "
			<< setw(20) << "Gx" << " "
			<< setw(20) << "Gy" << " "
			<< setw(20) << "Gz" << " "

			<< endl;
		
		file << setw(20) << TIME << " "
			<< setw(20) << Point_X << " "
			<< setw(20) << Point_Y << " "
			<< setw(20) << Point_Z << " "
			<< setw(20) << Speed_X << " "
			<< setw(20) << Speed_Y << " "
			<< setw(20) << Speed_Z << " "
			<< setw(20) << omegaX << " "
			<< setw(20) << omegaY << " "
			<< setw(20) << omegaZ << " "
			<< setw(20) << ALFA * TO_DEG << " "
			<< setw(20) << BETTA * TO_DEG<< " "
			<< setw(20) << PITCH * TO_DEG << " "
			<< setw(20) << YAW * TO_DEG << " "
			<< setw(20) << ROLL * TO_DEG << " "
			<< setw(20) << TETA * TO_DEG << " "
			<< setw(20) << PSI * TO_DEG << " "
			<< setw(20) << B_ * TO_DEG << " "
			<< setw(20) << L_ * TO_DEG << " "
			<< setw(20) << H_ << " "
			<< setw(20) << DeltaPitch * TO_DEG << " "
			<< setw(20) << DeltaYaw * TO_DEG << " "
			<< setw(20) << DeltaRoll * TO_DEG << " "
			<< setw(20) << Gx << " "
			<< setw(20) << Gy << " "
			<< setw(20) << Gz << " "

			<< endl;
		
		file.close();
	}
	else
	{
		ofstream file("file_result.txt", ios::app);
		if (file.is_open())
		{
					file << setw(20) << TIME << " "
						<< setw(20) << Point_X << " "
						<< setw(20) << Point_Y << " "
						<< setw(20) << Point_Z << " "
						<< setw(20) << Speed_X << " "
						<< setw(20) << Speed_Y << " "
						<< setw(20) << Speed_Z << " "
						<< setw(20) << omegaX << " "
						<< setw(20) << omegaY << " "
						<< setw(20) << omegaZ << " "
						<< setw(20) << ALFA * TO_DEG << " "
						<< setw(20) << BETTA * TO_DEG << " "
						<< setw(20) << PITCH * TO_DEG << " "
						<< setw(20) << YAW * TO_DEG << " "
						<< setw(20) << ROLL * TO_DEG << " "
						<< setw(20) << TETA * TO_DEG << " "
						<< setw(20) << PSI * TO_DEG << " "
						<< setw(20) << B_ * TO_DEG << " "
						<< setw(20) << L_ * TO_DEG << " "
						<< setw(20) << H_ << " "
						<< setw(20) << DeltaPitch * TO_DEG << " "
						<< setw(20) << DeltaYaw * TO_DEG << " "
						<< setw(20) << DeltaRoll * TO_DEG << " "
						<< setw(20) << Gx << " "
						<< setw(20) << Gy << " "
						<< setw(20) << Gz << " "
						<< endl;
		}
		
		file.close();
	}
	
}

void Print_zone()
{
	ofstream file("print_zone.txt", ios::app);
	if (file.is_open())
	{
		file 
			<< setw(20) << target_X << " "
			<< setw(20) << target_Z << " "
			<< setw(20) << Point_X << " "
			<< setw(20) << Point_Z << " "
			<< endl;
	}

	file.close();
}

