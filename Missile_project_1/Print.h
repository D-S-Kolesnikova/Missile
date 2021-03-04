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

		file << setw(20) << " Time,ρ" << " "
			<< setw(20) << "υ,κμ" << " "
			<< setw(20) << "y,κμ" << " "
			<< setw(20) << "z,κμ" << " "
			<< setw(20) << "Speed_X" << " "
			<< setw(20) << "Speed_Y" << " "
			<< setw(20) << "Speed_Z" << " "
			<< setw(20) << "Omegax" << " "
			<< setw(20) << "Omegay" << " "
			<< setw(20) << "Omegaz" << " "
			<< setw(20) << "Alfa" << " "
			<< setw(20) << "Betta" << " "
			<< setw(20) << "Pitch" << " "
			<< setw(20) << "Yaw" << " "
			<< setw(20) << "Roll" << " "
			<< setw(20) << "TETA" << " "
			<< setw(20) << "PSI" << " "
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
						<< endl;
		}
		
		file.close();
	}
	
}

