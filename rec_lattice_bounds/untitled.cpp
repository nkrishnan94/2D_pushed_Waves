#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


const xdemes = 5000;
const ydemes = 5000;
const zdemes = 5000;
unsigned int n_gens = 10000;


long double deme[xdemes][ydemes][zdemes] = {{0}};

int main(){




	srand (time(NULL));

	/* generate secret number between 1 and 10: */
	iSecret = rand() % 10 + 1;

	ofstream fprof;
	time_t time_start;
	clock_t c_init = clock();
	struct tm * timeinfo;
	char buffer [80];
	float ht = 0.5;
	int prof_count = 4;


	time (&time_start);
	timeinfo = localtime (&time_start);
	strftime (buffer,80,"%F-%H-%M-%S",timeinfo);

	ostringstream date_time, Kstr,  Mstr, Bstr, Gstr;
	date_time << buffer;

	fprof.open("prof_"+date_time.str() + ".txt")

	for(int i = int(xdemes*.5)-10; i < int(xdemes*.5)+10; i++){
		for(int j = int(ydemes*.5)-10; j < int(ydemes*.5)+10; j++){


			deme[i][j][0] = 1;
			deme[i][j][0] = 1;


		}

	}

	for (int dt = 0 ; dt < n_gens; dt++ ){

		for(int i = 0; i < xdemes; i++){
			for(int j = 0; j < ydemes; j++){
				for(int k = 0; k < zdemes; j++){
					pick = rand() % 6 ;


				}
			}
		}
	}






}