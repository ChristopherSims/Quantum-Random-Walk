// Quantum Random Walk.cpp 
//

#include "stdafx.h"


using namespace std;
typedef complex<double> dcomp;

int main(){
	double const PI = 3.14;
	int const steps = 3;
	int const repetitions = 20;
	int const probability = 10;
	int const max = 103 * 103;
	int const itr = 103;
	int sx[probability+1][repetitions];
	int sy[probability+1][repetitions];
	vector<int> x;
	vector<int> y;
	int sumx = 0;
	int sumy = 0;
	complex<double> img = sqrt(complex<double>(-1));
	dcomp phase[103][103];
	/*start of operation*/
	for (int p = 0; p < probability; p++) { // p
		for (int n = 1; n < repetitions; n++){ // n
			for (int i = 0; i < probability + 1; i++){
				for (int j = 0; j < repetitions; j++){
					sx[i][j] = 0;
					sy[i][j] = 0;
				}
			}

			/*
			Initialize quantum state variables
			*/
			dcomp phase[103][103];
			vector<dcomp> C1(max);
			vector<dcomp> C2(max);
			vector<dcomp> C3(max);
			vector<dcomp> C4(max);
			vector<dcomp> C1f(max);
			vector<dcomp> C2f(max);
			vector<dcomp> C3f(max);
			vector<dcomp> C4f(max);
			vector<double> mv(max);
			vector<double> w2(max);
			C1[(103 * 52) + 52] = 1.0;

			/*
			creates phase state
			*/
			for (int i = 0; i < 103; i++){
				for (int j = 0; j < 103; j++){
					phase[i][j] = 2 * PI * (rand()/(double)(RAND_MAX + 1));
				}
			}
			
			/*
			quantum step routine 
			*/
		 for (int k = 0; k < steps; k++){
				/*
				// Performs the hadamard and flip operations
				*/
			for (int i = 2; i < itr - 2 ; i++){
				for (int j = 2; j < itr - 2; j++){
					C1f[(103 * i) + j] = operator*(exp(operator*(img, phase[i][j])),
						0.5 * operator+(operator+(std::operator+(C1[(103 * (i - 1)) + (j - 1)], C2[(103 * (i - 1)) + (j - 1)]), C3[(103 * (i - 1)) + (j - 1)]), C4[(103 * (i - 1)) + (j - 1)]));
					C2f[(103 * i) + j] = operator*(exp(operator*(img, phase[i][j])),
						 0.5 * operator+(operator+(operator+(C1[(103 * (i - 1)) + (j + 1)], C2[(103 * (i - 1)) + (j + 1)]), C3[(103 * (i - 1)) + (j + 1)]), C4[(103 * (i - 1)) + (j + 1)]));
					C3f[(103 * i) + j] = operator*(exp(operator*(img, phase[i][j])),
					 0.5 *  operator+(operator+(operator+(C1[(103 * (i + 1)) + (j - 1)], C2[(103 * (i + 1)) + (j - 1)]), C3[(103 * (i + 1)) + (j - 1)]), C4[(103 * (i + 1)) + (j - 1)]));
					C4f[(103 * i) + j] = operator*(exp(operator*(img, phase[i][j])),
						0.5 * operator+(operator+(operator+(C1[(103 * (i + 1)) + (j + 1)], C2[(103 * (i + 1)) + (j + 1)]), C3[(103 * (i + 1)) + (j + 1)]), C4[(103 * (i + 1)) + (j + 1)]));
				}
			}
			/* Adds final values of origional state array*/
			for (int i = 0; i < 103; i++){
				for (int j = 0; j <  103; j++){
					C1[(103 * i) + j] = C1f[(103 * i) + j];
					C2[(103 * i) + j] = C2f[(103 * i) + j];
					C3[(103 * i) + j] = C3f[(103 * i) + j];
					C4[(103 * i) + j] = C4f[(103 * i) + j];
				}
			}
			
			for (int i = 1; i < itr - 1; i++){
				for (int j = 1; j < itr - 1; j++){
					w2[(103 * (i - 1)) + j] = 
						(pow(abs(C1[(103 * i) + j]), 2)
						+ pow(std::abs(C2[(103 * i) + j]), 2)
						+ pow(std::abs(C3[(103 * i) + j]), 2) 
						+ pow(std::abs(C4[(103 * i) + j]), 2));
				}
			}
			

			/*simulation of the probability of a measurment\
			
			*/
			double measurement_vairable1 = (rand() / (double)(RAND_MAX + 1));
			if (measurement_vairable1 < (p / probability)){
				mv[0] = 0;
				mv[10608] = 0;
				mv[2] = w2[2];
				for (int i = 1; i < max - 1; i++){
					mv[i] = mv[i - 1] + w2[i];

				}
				double measurment_variable2 = (rand()/(double)(RAND_MAX + 1));
				int measurement_pos = 1;
					while (measurment_variable2 >= mv[measurement_pos] && measurement_pos < (max - 1)){
						measurement_pos++;
					}
				double measurement_variable3 = (rand()/(double)(RAND_MAX + 1));
				double c1_norm = pow(std::abs(C1f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)]), 2);
					double c2_norm = pow(std::abs(C2f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)]), 2);
						double c3_norm = pow(std::abs(C3f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)]), 2);
							double c4_norm = pow(std::abs(C4f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)]), 2);
				double normal_square = sqrt(c1_norm + c2_norm + c3_norm + c4_norm);
				double inv_norm_square = pow((1 / normal_square), 2);
				if ((inv_norm_square * c1_norm) < measurement_variable3 && measurement_variable3 <= (inv_norm_square * (c1_norm + c2_norm))){
					C2f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)] = 1;
				}
				if (inv_norm_square * (c1_norm + c2_norm) < measurement_variable3 && measurement_variable3 <= inv_norm_square * (c1_norm + c2_norm + c3_norm)){
					C2f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)] = 1;
				}
				if (measurement_variable3 >= (1 - (inv_norm_square * c4_norm))){
					C4f[(itr * div(measurement_pos, 103).quot + 1) + ((measurement_pos - 1 % 103) + 1)] = 1;
				}
				for (int i = 1; i < itr; i++){
					for (int j = 1; j < itr; j++){
						w2[(itr * (i - 1)) + j] =
							pow(abs(C1[(103 * i) + j]), 2) +
							pow(abs(C2[(103 * i) + j]), 2) +
							pow(abs(C3[(103 * i) + j]), 2) +
							pow(abs(C4[(103 * i) + j]), 2);
					}
				}

				/* end of p/prob if statement*/
			}
			/* end of steps loop*/
		}
		mv[1] = 0;
		mv[10608] = 0;
		mv[2] = w2[2];
		for (int i = 2; i < max - 2; i++){
			mv[i] = mv[i - 1] + w2[i];
		}
		double measurment_variable2 = (rand() / (double)(RAND_MAX + 1));;
		
		int measurment_pos = 1;
		while (measurment_variable2 >= mv[measurment_pos] && measurment_pos < (max - 1)){
			measurment_pos++;
		}
		int a = 52;
		
		 sx[p+1][n] = (int) pow(abs(a - (measurment_pos / 103) + 1), 2);
		 sy[p+1][n] = (int) pow(abs(a - ((measurment_pos - 1) % 103) + 1), 2);
		sumy += sy[p+1][n];
		sumx += sx[p+1][n];
		//cout << sumx << endl;
	cout << "P:" << p << " N: " << n << endl;
	mv.clear();
	w2.clear();
	C1.clear();
	 C2.clear();
	 C3.clear();
	 C4.clear();
	 C1f.clear();
	C2f.clear();
	C3f.clear();
	 C4f.clear();
	
		}
		y.push_back(sumy);
		x.push_back(sumx);
		sumy = 0;
		sumx = 0;
	}

	vector<double> rmsx;
	vector<double> rmsy;
	for (int i = 0; i < probability; i++){
		rmsx.push_back(i / probability);
	}
	double sumsq = 0;
	float RMS = 0;
	int rep = repetitions;
	double inv = (1.0/rep);
	for (int h = 0; h < probability; h++){
		sumsq = y[h] + x[h];
		RMS = sqrt(inv * sumsq);
			rmsy.push_back(RMS);
	}
	for (int i = 0; i < probability; i++){
		cout << "X: " << rmsx[i] << " " << "Y: " << rmsy[i] << endl;
	}
	system("pause");

	return 0;
}

