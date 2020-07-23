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


unsigned long K  = 50000;
unsigned int n_gens = 1000000;
unsigned int n_gens = 1000000;
const int n_demes = 80; 
const unsigned int n_spec = 2;
float M = 0.25;
float B = 8;
float g0 = 0.01;
unsigned long prof_hist = 0;
unsigned long fast_samp_flag = 1;



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

	return sum/K;




}




float calcHet(long double arr[][n_demes][n_spec], const int arrSize){




	int cnt =  0 ;
	long double H = 0.0;

	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){

			double deme_pop = arr[i][j][0]+arr[i][j][1];


			//std::cout << i << "\n";
			if (deme_pop > 0.0){
				H += (2*arr[i][j][0]*arr[i][j][1])/(deme_pop*deme_pop);
				//std::cout << arr[i][0] << "\n";
				//std::cout << (2*arr[i][0]*(deme_pop - arr[i][0]))/(deme_pop*deme_pop)<< "\n";
				cnt+=1;


			}

		}

	}





	return  H/cnt;

}

float calcVarHet(long double arr[][n_demes][n_spec], const int arrSize){

	long double hets[n_demes][n_demes];
	long double H = 0.0;
	long double varH = 0.0;
	int cnt=0;
	float average;


	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){

			double deme_pop = arr[i][j][0]+arr[i][j][1];
			//std::cout << i << "\n";
			if (deme_pop > 0.0){
				hets[i][j]= (2*arr[i][j][0]*arr[i][j][1])/(deme_pop*deme_pop);
				H+=hets[i][j];
				cnt+=1;
				//std::cout << arr[i][0] << "\n";
				//std::cout << (2*arr[i][0]*(deme_pop - arr[i][0]))/(deme_pop*deme_pop)<< "\n";
			}

		}
	}

	average =  H/cnt;
	cnt=0;
	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){
			double deme_pop = arr[i][j][0]+arr[i][j][1];
			if (deme_pop > 0.0){
				varH+= (hets[i][j]-average)*(hets[i][j]-average);
				cnt+=1;



			}

		}
	}

	

	return  varH/cnt;

}

long double deme[n_demes][n_demes][n_spec] = {{0}};
long double deme_aux[n_demes][n_demes][n_spec] = {{0}};






int main (int argc, char * argv[]){
	using namespace std;

	int c;
    while ((c = getopt (argc, argv, "K:Z:B:M:G:F")) != -1)
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
        else if (c == 'F')
            fast_samp_flag = atoi(optarg); // growth rate

    }
    //n_gens = 50*K;
    /*if (B >= 2.0)
    	n_gens = 50*K;
    else
        n_gens = 25*K;*/


	const gsl_rng_type * T;
	gsl_rng * r;



	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	//int sysRandom;
	gsl_rng_set(r, time(NULL));

	double new_prob[n_spec + 1];
	unsigned int new_cnt[n_spec + 1];
	int n_data = 1000;
	int record_time = int(n_gens/n_data);
	//int record_time = 10;
	//int n_data = int(n_gens/record_time);

	double pop_shift = 0.0;
	double w_s = 1.0;
	double w_avg;
	double w_v;
	vector <double> pop_hist;
	vector <double> het_hist;

	//data files
	ofstream flog, fpop, fhet, fprof;
	time_t time_start;
	clock_t c_init = clock();
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

    flog.open( logName);
    fhet.open(hetName);
    fpop.open(popName);
    fprof.open(profName);








	for(int i = 0; i < int(n_demes*.5); i++){
		for(int j = 0; j < n_demes; j++){

		deme[i][j][0] = .5*K;
		deme[i][j][1] = .5*K;


		}


	}
	//initial population in middle
	//deme[int(n_demes/2)][int(n_demes/2)][0] = K;
	//deme[int(n_demes/2)][int(n_demes/2)][1] = K;




    /*if (prof_hist !=0){
		ostringstream strT;

		string proftName = "prof_T_start_" + date_time.str() + ".txt";
		ofstream fproft;
	    fproft.open(proftName);
	    for(int i = 0; i <n_demes; i++){
	    	for(int j = 0; j <n_demes; j++){
	    		fproft << i << ", " << j << ", " << deme[i][j][0] << ", " << deme[i][j][1] <<endl;


	    	}
		}

    }*/

	for (int dt = 0 ; dt < n_gens; dt++ ){

		
		for(int ii = 0; ii < int(n_demes); ii++){
			for(int jj = 0; jj < int(n_demes); jj++){


				deme_aux[ii][jj][0] = deme[ii][jj][0];
				deme_aux[ii][jj][1] = deme[ii][jj][1];

			

			}

		}
		//int i=0;
		for( int j=0; j<n_demes; j++){
			int arr[2] = {0, j};
			int neighb_vec[3][2] = {{0,1},{1,0},{0,-1}};
			int neighbs[3][2];
			int pop_sum = fast_samp_flag;
			
			for(int ne=0; ne <3; ne++){
				neighbs[ne][0] = (arr[0] + n_demes+neighb_vec[ne][0]) % n_demes;
				neighbs[ne][1] = (arr[1] + n_demes+neighb_vec[ne][1]) % n_demes;
				//cout << neighbs[ne][0] << ", " << neighbs[ne][1] << endl;

			}


			long double f1 = deme[0][j][0]/int(K);
			long double f2 = deme[0][j][1]/int(K);
			f1 = (1 - M*(3/4))*f1;
			f2 = (1 - M*(3/4))*f2; 


			for(int ne = 0; ne <3; ne++){

				f1+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][0]/int(K);
				f2+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][1]/int(K);


			}

			w_v = 1 - g0*(1+ B*(f1+f2));
			w_avg = w_v + (w_s - w_v)*(f1 + f2);

			f1 *= w_s/w_avg;
			f2 *= w_s/w_avg;
			new_prob[0] = 1-f1-f2;
			new_prob[1] = f1;
			new_prob[2] = f2;
			gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
			deme[0][j][0] = new_cnt[1];
			deme[0][j][1] = new_cnt[2];
		}

		for(int i = 1; i < n_demes-2 ; i++){
			for(int j = 0; j < n_demes; j++){
				int arr[2] = {i, j};
				int neighb_vec[4][2] = {{0,1},{0,-1},{1,0},{-1,0}};
				int neighbs[4][2];
				int pop_sum = fast_samp_flag;

				for(int ne=0; ne <4; ne++){
					neighbs[ne][0] = (arr[0] + n_demes+neighb_vec[ne][0]) % n_demes;
					neighbs[ne][1] = (arr[1] + n_demes+neighb_vec[ne][1]) % n_demes;
					//pop_sum+=deme_aux[neighbs[ne][0]][neighbs[ne][1]][0] +deme_aux[neighbs[ne][0]][neighbs[ne][1]][1];
					//pop_sum+=deme[neighbs[ne][0]][neighbs[ne][1]][0] +deme[neighbs[ne][0]][neighbs[ne][1]][1];

					/*int n_neighbs[4][2];
					int ne_arr[2] = {neighbs[ne][0], neighbs[ne][1]};

					
					for(int ne_=0; ne_ <4; ne_++){
						n_neighbs[ne_][0] = (ne_arr[0] + n_demes+neighb_vec[ne_][0]) % n_demes;
						n_neighbs[ne_][1] = (ne_arr[1] + n_demes+neighb_vec[ne_][1]) % n_demes;
						pop_sum+=deme_aux[n_neighbs[ne_][0]][n_neighbs[ne_][1]][0] + deme_aux[n_neighbs[ne_][0]][n_neighbs[ne_][1]][1];
						pop_sum+=deme[n_neighbs[ne_][0]][n_neighbs[ne_][1]][0] + deme[n_neighbs[ne_][0]][n_neighbs[ne_][1]][1];

					}*/




				}


				if (((deme[i][j][0] + deme[i][j][1]+deme_aux[i][j][1]+deme_aux[i][j][0]) !=0) || (pop_sum != 0)){

					long double f1 = deme[i][j][0]/int(K);
					long double f2 = deme[i][j][1]/int(K);
					f1 = (1 - M)*f1;
					f2 = (1 - M)*f2; 
					for(int ne = 0; ne <4; ne++){

						f1+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][0]/int(K);
						f2+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][1]/int(K);


					}

					w_v = 1 - g0*(1+ B*(f1+f2));
					w_avg = w_v + (w_s - w_v)*(f1 + f2);

					f1 *= w_s/w_avg;
					f2 *= w_s/w_avg;
					new_prob[0] = 1-f1-f2;
					new_prob[1] = f1;
					new_prob[2] = f2;
					gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
					deme[i][j][0] = new_cnt[1];
					deme[i][j][1] = new_cnt[2];




				}



				

			}
		}

		int i=n_demes-1;
		for( int j=0; j<n_demes; j++){
			int arr[2] = {i, j};
			int neighb_vec[3][2] = {{0,1},{-1,0},{0,-1}};
			int neighbs[3][2];
			int pop_sum = fast_samp_flag;
			//cout << i << ", " << j << endl;

			for(int ne=0; ne <3; ne++){
				neighbs[ne][0] = (arr[0] + n_demes+neighb_vec[ne][0]) % n_demes;
				neighbs[ne][1] = (arr[1] + n_demes+neighb_vec[ne][1]) % n_demes;
				//cout << neighbs[ne][0] << ", " << neighbs[ne][1] << endl;

			}


			long double f1 = deme[i][j][0]/int(K);
			long double f2 = deme[i][j][1]/int(K);
			f1 = (1 - M*(3/4))*f1;
			f2 = (1 - M*(3/4))*f2; 


			for(int ne = 0; ne <3; ne++){

				f1+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][0]/int(K);
				f2+= (M/4)*deme_aux[neighbs[ne][0]][neighbs[ne][1]][1]/int(K);


			}

			w_v = 1 - g0*(1+ B*(f1+f2));
			w_avg = w_v + (w_s - w_v)*(f1 + f2);

			f1 *= w_s/w_avg;
			f2 *= w_s/w_avg;
			new_prob[0] = 1-f1-f2;
			new_prob[1] = f1;
			new_prob[2] = f2;
			gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
			deme[i][j][0] = new_cnt[1];
			deme[i][j][1] = new_cnt[2];
		}

		

		




		if (sumDeme(deme,n_demes) > .5*n_demes*n_demes){
			int shift =  int(sumDeme(deme,n_demes) - .5*n_demes*n_demes)+1;
			//cout << shift << endl;
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
			
			
	        for (int i = int(n_demes - shift); i < n_demes; i++){
	        	for (int j = 0; j < n_demes; j++){ 

				    for (int k = 0; k < n_spec; k++){
				    	//cout <<i << endl;
				        deme[i][j][k] = 0;
				    }
			    }
	        }


			

	        pop_shift += shift;





		}


		//cout << dt << endl;
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

    	fhet << int(i*record_time) << ", "  << het_hist[i] << endl;
    	fpop << int(i*record_time) << ", "  << pop_hist[i] << endl;

    }

    for(int i=0;i < n_demes; i++){
    	for(int j =0;j < n_demes; j++){

    		fprof << i << ", " << j<< ", " << deme[i][j][0] << ", "  << deme[i][j][1] << endl;

		}

    }

    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demes << time_start<< run_time<< endl;

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