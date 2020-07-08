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


unsigned long K  = 1000;
unsigned int n_gens = 1000;
const int n_demes = 300; 
const unsigned int n_spec = 2;
float M = 0.25;
float B = 2;
float g0 = 0.01;
unsigned long prof_hist = 0;


double sumDeme(long double arr[][n_demes][n_spec], int arrSize){
			// sum of population
	double sum = 0.0;
	for (int i=0; i<n_demes; i++){
		for (int j=0; j<n_demes; j++){

			for (int k = 0; k<n_spec; k++){
				sum+= arr[i][j][k];

			};
		};
	};

	return sum;




}




float calcHet(long double arr[][n_demes][n_spec], const int arrSize){




	int cnt =  0 ;
	long double H = 0.0;

	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){

			double deme_pop = arr[i][j][0]+arr[i][j][1];


			//std::cout << i << "\n";
			if (deme_pop > 0.0){
				H += (2*arr[i][j][0]*(deme_pop - arr[i][0][j]))/(deme_pop*deme_pop);
				//std::cout << arr[i][0] << "\n";
				//std::cout << (2*arr[i][0]*(deme_pop - arr[i][0]))/(deme_pop*deme_pop)<< "\n";
				cnt+=1;


			}

		}

	}





	return  H/cnt;

}





int main (int argc, char * argv[]){
	using namespace std;

	int c;
    while ((c = getopt (argc, argv, "K:Z:B:M:G")) != -1)
    {
        if (c == 'K')
            K  = atoi(optarg); // carrying capacity
        else if (c == 'Z')
            prof_hist = atoi(optarg); //keep track of profile through time 
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        else if (c == 'M')
            M = atof(optarg); // migration probability
        else if (c == 'G')
            g0 = atof(optarg); // growth rate

    }


	const gsl_rng_type * T;
	gsl_rng * r;
	long double deme[n_demes][n_demes][n_spec] = {{0}};
	long double deme_aux[n_demes][n_demes][n_spec] = {{0}};


	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	//int sysRandom;
	gsl_rng_set(r, time(NULL));

	double new_prob[n_spec + 1];
	unsigned int new_cnt[n_spec + 1];
	//int n_data = 1000;
	//int record_time = int(n_gens/n_data);
	int record_time = 10;
	int n_data = int(n_gens/record_time);

	double pop_shift = 0.0;
	double w_s = 1.0;
	double w_avg;
	double w_v;
	vector <double> pop_hist;
	vector <double> het_hist;

	//data files
	ofstream flog, fpop, fhet, fprof;
	time_t time_start;
	struct tm * timeinfo;
	char buffer [80];


    time (&time_start);
	timeinfo = localtime (&time_start);

	strftime (buffer,80,"%F-%H-%M-%S",timeinfo);

	ostringstream date_time, Kstr,  Mstr, Bstr, Gstr;
	date_time << buffer;
	Kstr << K;
	Mstr << M;
	Bstr << B;
	Gstr << g0;
	string param_string =  "K"+Kstr.str()+"_M" + Mstr.str() + "_B" +Bstr.str() + "_G" +Gstr.str() +"_";



	string logName = "log_" + param_string + date_time.str() + ".txt";
	string hetName = "het_" + param_string +  date_time.str() + ".txt";
	string popName = "pop_"+ param_string +  date_time.str() + ".txt";
	string profName = "prof_" + param_string + date_time.str() + ".txt";

    flog.open(logName);
    fhet.open(hetName);
    fpop.open(popName);
    fprof.open(profName);







    //initial population in one side
	/*for(int i = 0; i < int(n_demes*.5); i++){
		for(int j = 0; j < int(n_demes); j++){

		deme[i][j][0] = .5*K;
		deme[i][j][1] = .5*K;


		}


	}*/
	//initial population in middle
	deme[int(n_demes/2)][int(n_demes/2)][0] = .5*K;
	deme[int(n_demes/2)][int(n_demes/2)][1] = .5*K;





	for (int dt = 0 ; dt < n_gens; dt++ ){
		for(int ii = 0; ii < int(n_demes); ii++){
			for(int jj = 0; jj < int(n_demes); jj++){


			deme_aux[ii][jj][0] = deme[ii][jj][0];
			deme_aux[ii][jj][1] = deme[ii][jj][1];

			}


		}
		


		for(int i = 0; i < n_demes ; i++){
			for(int j = 0; j < n_demes; j++){
				int arr[2] = {i, j};
				int neighb_vec[4][2] = {{0,1},{0,-1},{1,0},{-1,0}};
				std::vector<std::vector <int> > neighbs;

	
				for(int ne=0; ne <4; ne++){

					if( ((arr[0] + neighb_vec[ne][0])<n_demes && (arr[0] + neighb_vec[ne][0])>= 0 )){
						std::vector<int> temp;
						
						temp.push_back(int(arr[0] + neighb_vec[ne][0]));
						temp.push_back( int(arr[1] + neighb_vec[ne][1]) % n_demes);
						

						neighbs.push_back(temp);


					}


				}
				
				

				

			}
		}

		

		




		if (sumDeme(deme,n_demes)/K > .5*n_demes*n_demes){
			int shift =  int(sumDeme(deme,n_demes)/K - .5*n_demes*n_demes)+1;
			if ((shift< 0) == true){

	            cout << "Negative shift of array!" << endl;
	            exit(EXIT_FAILURE);

	        }


			for (int i = 0; i < n_demes - shift; i++){
				for (int j = 0; j < n_demes; j++){ 

					for (int k = 0; k <n_spec;k++){

						deme[i][j][k] = deme[i+shift][j][k];


					}
				}
			}


	        for (int i = n_demes - shift; i < n_demes; i++){
	        	for (int j = 0; j < n_demes; j++){ 

				    for (unsigned int k = 0; k < n_spec; k++){
				        deme[i][j][k] = 0;
				    }
			    }
	        }


			

	        pop_shift += shift*K*n_demes;





		}


		
		//cout << record_time << endl;
		if (dt % record_time == 0){

			het_hist.push_back(calcHet(deme, n_demes)); // Store heterozygosity
	        pop_hist.push_back(pop_shift+sumDeme(deme,n_demes));

	        if (prof_hist !=0){
	        	ostringstream strT;
	        	strT << dt;
	        	string proftName = "prof_T"+ strT.str() + "_" + date_time.str() + ".txt";
	        	ofstream fproft;
	            fproft.open(proftName);
	            for(int i = 0; i <n_demes; i++){
	            	for(int j = 0; j <n_demes; j++){
	            		fproft << i << ", " << j << ", " << deme[i][j][0] << ", " << deme[i][j][1] <<endl;


	            	}
            	}

	        }


		}
	



    }

   for(int i=0;i < n_data;i++){

    	fhet << i*record_time << ", "  << het_hist[i] << endl;
    	fpop << i*record_time << ", "  << pop_hist[i] << endl;

    }

    for(int i=0;i < n_demes; i++){
    	for(int j =0;j < n_demes; j++){

    		fprof << i << ", " << j<< ", " << deme[i][j][0] << ", "  << deme[i][j][1] << endl;

		}

    }

    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demes << time_start<< endl;

    fhet.close();
    fpop.close();
    flog.close();
    fprof.close();

    cout << "Finished!" << "\n";
    cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);









	return 0;
	




}