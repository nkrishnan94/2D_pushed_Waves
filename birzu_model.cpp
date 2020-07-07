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
unsigned int n_gens = 10000;
const int n_demes = 300; 
const unsigned int n_spec = 2;
float M = 0.25;
float B = 2;
float g0 = 0.01;
int prof_hist = 0;


double sumDeme(long double arr[][n_spec], int arrSize){
			// sum of population
	double sum = 0.0;
	for (int i=0; i<n_demes; i++){
		for (int j = 0; j<n_spec; j++){
			sum+= arr[i][j];

		};
	};

	return sum;




}




float calcHet(long double arr[][n_spec], const int arrSize){




	int cnt =  0 ;
	long double H = 0.0;

	for(int i = 0; i < arrSize; i++){

		float deme_pop = arr[i][0]+arr[i][1];


		//std::cout << i << "\n";
		if (deme_pop > 0.0){
			H += (2*arr[i][0]*(deme_pop - arr[i][0]))/(deme_pop*deme_pop);
			//std::cout << arr[i][0] << "\n";
			//std::cout << (2*arr[i][0]*(deme_pop - arr[i][0]))/(deme_pop*deme_pop)<< "\n";
			cnt+=1;


		}



	}





	return  H/cnt;

}



int main (int argc, char * argv[]){
	using namespace std;

	int c;
    while ((c = getopt (argc, argv, "K:B:G:M:H")) != -1)
    {
        if (c == 'K')
            K  = atoi(optarg); // carrying capacity
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        else if (c == 'M')
            M = atof(optarg); // migration probability
        else if (c == 'G')
            g0 = atof(optarg); // growth rate
        else if (c == 'H')
            prof_hist = atoi(optarg); //keep track of profile through time 
    }


	const gsl_rng_type * T;
	gsl_rng * r;
	long double deme[n_demes][n_spec] = {{0}};


	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	//int sysRandom;
	cout << "Random seed: " << time(NULL)<< endl;
	gsl_rng_set(r, time(NULL));
	

	double new_prob[n_spec + 1];
	unsigned int new_cnt[n_spec + 1];
	int n_data = 1000;
	int record_time = int(n_gens/n_data);
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








	for(int i = 0; i < int(n_demes*.5); i++){

		deme[i][0] = gsl_ran_binomial(r, 0.5, K);
		deme[i][1] = K - deme[i][0];


	}




	for (int dt = 0 ; dt < n_gens; dt++ ){
		long double f1 = deme[0][0]/int(K);
		long double f2 = deme[0][1]/int(K);
		long double f1_ = f1;
		long double f2_ = f2;
		f1 = (1 - M/2)*f1 + (M/2) *deme[1][0]/int(K);
		f2 = (1 - M/2)*f2 + (M/2) *deme[1][1]/int(K);
		long double f1_new = f1_;
		long double f2_new = f2_;
		w_v = 1 - g0*(1+ B*(f1+f2));
		w_avg = w_v + (w_s - w_v)*(f1 + f2);

		f1 *= w_s/w_avg;
		f2 *= w_s/w_avg;
		new_prob[0] = 1-f1-f2;
		new_prob[1] = f1;
		new_prob[2] = f2;
		gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
		deme[0][0] = new_cnt[1];
		deme[0][1] = new_cnt[2];





		for(int i = 1; i < n_demes - 1; i++){

			f1 = deme[i][0]/int(K);
			f2 = deme[i][1]/int(K);
			f1_ = f1;
			f2_ = f2;
			f1 = (1 - M)*f1 + (M/2) *deme[i+1][0]/int(K)+ (M/2)*f1_new;
			f2 = (1 - M)*f2 + (M/2) *deme[i+1][1]/int(K)+ (M/2)*f2_new;
			f1_new = f1_;
			f2_new = f2_;
			w_v = 1 - g0*(1+ B*(f1+f2));
			w_avg = w_v + (w_s - w_v)*(f1 + f2);

			f1 *= w_s/w_avg;
			f2 *= w_s/w_avg;
			new_prob[0] = 1-f1-f2;
			new_prob[1] = f1;
			new_prob[2] = f2;
			gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
			deme[i][0] = new_cnt[1];
			deme[i][1] = new_cnt[2];
		};

		int i = n_demes-1;

		f1 = deme[i][0]/int(K);
		f2 = deme[i][1]/int(K);
		f1 = (1 - M/1)*f1 +  (M/2)*f1_new;
		f2 = (1 - M/1)*f2 + (M/2)*f2_new;
		w_v = 1 - g0*(1+ B*(f1+f2));
		w_avg = w_v + (w_s - w_v)*(f1 + f2);

		f1 *= w_s/w_avg;
		f2 *= w_s/w_avg;
		new_prob[0] = 1-f1-f2;
		new_prob[1] = f1;
		new_prob[2] = f2;
		gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
		deme[i][0] = new_cnt[1];
		deme[i][1] = new_cnt[2];




		if (sumDeme(deme,n_demes)/K > .5*n_demes){
			int shift =  int(sumDeme(deme,n_demes)/K - .5*n_demes)+1;
			if ((shift< 0) == true){

	            cout << "Negative shift of array!" << endl;
	            exit(EXIT_FAILURE);

	        }


			for (int i = 0; i < n_demes - shift; i++){
				for (int j = 0; j <n_spec;j++){

					deme[i][j] = deme[i+shift][j];


				}
			}


	        for (int i = n_demes - shift; i < n_demes; i++){
			    for (unsigned int j = 0; j < n_spec; j++){
			        deme[i][j] = 0;
			    }
	        }


			

	        pop_shift += shift*K;





		}



		if (dt % record_time == 0){

			het_hist.push_back(calcHet(deme, n_demes)); // Store heterozygosity
	        pop_hist.push_back(pop_shift+sumDeme(deme,n_demes));

	        if (prof_hist !=0){
	        	ostringstream strT;
	        	strT << dt;
	        	string proftName = "prof_T"+ strT.str() + "_" + date_time.str() + ".txt";
	        	ofstream fproft;
	            fproft.open(proftName);
	            for(int i = 0; i< n_demes;i++){
	            	fproft << i << ", " << deme[i][0] << ", " << deme[i][0] << endl;

	            }
	            


	        }


		}
	



    }

    for(int i=0;i < n_data;i++){

    	fhet << i*record_time << ", "  << het_hist[i] << endl;
    	fpop << i*record_time << ", "  << pop_hist[i] << endl;

    }

    for(int i=0;i < n_demes; i++){

    	fprof << i << ", "  << deme[i][0] << ", "  << deme[i][1] << endl;



    }
    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demes << time_start<< run_time <<endl;

    fhet.close();
    fpop.close();
    flog.close();
    fprof.close();

    cout << "Finished!" << "\n";
    cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	//puts (buffer);









	return 0;
	




}