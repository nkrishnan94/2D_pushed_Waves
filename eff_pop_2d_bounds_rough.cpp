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


unsigned long K  = 1000000;
unsigned int n_gens = 10;
float h_thresh = .3;
const int n_demesh = 120; 
const int n_demesw = 80; 
const unsigned int n_spec = 2;
float M = 0.25;
float B = 0;
float g0 = 0.01;
unsigned long prof_hist = 0;
unsigned long fast_samp_flag = 1;
unsigned long freeze_flag = 0;



double sumDeme(long double arr[n_demesh][n_demesw][n_spec], int arrSize){
			// sum of population
	double sum = 0.0;
	for (int i=0; i<n_demesh; i++){
		for (int j=0; j<n_demesw; j++){

			for (int k = 0; k<n_spec; k++){
				sum+= arr[i][j][k];

			};
		};
	};

	return sum/K;




}




float calcHet(long double arr[n_demesh][n_demesw][n_spec], const int arrSize){




	int cnt =  0 ;
	long double H = 0.0;

	for(int i = 0; i < n_demesh; i++){
		for(int j=0; j<n_demesw;j++){

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

float calcRoughness(long double arr[n_demesw][n_spec]){
	long double freqs[n_demesw];
	float sumFreq=0;
	float sumDiffs=0;

	for(int i=0;i<n_demesw;i++){

		freqs[i] = arr[i][0]/ (arr[i][0]+arr[i][1]);
		sumFreq+=arr[i][0]/ (arr[i][0]+arr[i][1]);
	}

	for(int i=0;i<n_demesw;i++){
		sumDiffs+= abs(freqs[i] - sumFreq/n_demesw);

	}

	return sumDiffs/n_demesw;

}
int calcBoundPos(long double arr[n_demesw][n_spec]){
	long double freqs[n_demesw];
	float sumFreq=0;
	int boundDeme=0;
	float sumDiffs=0;

	for(int i=0;i<n_demesw;i++){

		freqs[i] = arr[i][0]/ (arr[i][0]+arr[i][1]);
		
	}
	int i =1;
	while((boundDeme==0)&&(i<n_demesw)){
	//for(int i=1;i<n_demesw;i++){
		if ((freqs[i]>.5 &&freqs[i-1]<.5) || (freqs[i]<.5 &&freqs[i-1]>.5)){
			boundDeme = i;


		}
		i+=1;

	}

	return boundDeme;

}



int calcLastRow(long double arr[n_demesh][n_demesw][n_spec]){
	int lastrow = 0;
	for(int i = 0; i < n_demesh; i++){
		int row_tot=0;

		for(int j=0; j < n_demesw; j++){
			if ((arr[i][j][0]+arr[i][j][1])>0){
				row_tot+=1;



			}


		}
		if (row_tot ==n_demesw){
			//std::cout<< row_tot<<std::endl;
			lastrow=i;



		}
	}

	return lastrow;


}

float calcVarHet(long double arr[n_demesh][n_demesw][n_spec], const int arrSize){

	long double hets[n_demesh][n_demesw];
	long double H = 0.0;
	long double varH = 0.0;
	int cnt=0;
	float average;


	for(int i = 0; i < n_demesh; i++){
		for(int j=0; j<n_demesw; j++){

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
	for(int i = 0; i < n_demesh; i++){
		for(int j=0; j<n_demesw;j++){
			double deme_pop = arr[i][j][0]+arr[i][j][1];
			if (deme_pop > 0.0){
				varH+= (hets[i][j]-average)*(hets[i][j]-average);
				cnt+=1;



			}

		}
	}



	return  varH/cnt;

}

long double deme[n_demesh][n_demesw][n_spec] = {{0}};
long double deme_aux[n_demesh][n_demesw][n_spec] = {{0}};






int main (int argc, char * argv[]){
	using namespace std;

	int c;
    while ((c = getopt (argc, argv, "K:Z:B:T:M:G:R:F")) != -1)
    {
        if (c == 'K')
            K  = atoi(optarg); // carrying capacity
        else if (c == 'Z')
            prof_hist = atoi(optarg); //keep track of profile through time 
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        else if (c == 'T')
            h_thresh = atof(optarg); // cooperativity
        else if (c == 'M')
            M = atof(optarg); // migration probability
        else if (c == 'G')
            g0 = atof(optarg); // growth rate
        else if (c == 'F')
            freeze_flag = atoi(optarg); 
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
	//int n_data = 1000;
	//int record_time = int(n_gens/n_data);
	int record_time = 10;
	int n_data = int(n_gens/record_time);
	//int dt = 0;

	double pop_shift = 0.0;
	double w_s = 1.0;
	double w_avg;
	double w_v;

	vector <double> pop_hist;
	vector <double> het_hist;
	vector <double> rough_0;
	vector <double> rough_10;
	vector <double> rough_20;
	vector <double> rough_30;
	vector <double> rough_40;
	vector <double> rough_50;
	vector <double> pos_0;
	vector <double> pos_10;
	vector <double> pos_20;
	vector <double> pos_30;
	vector <double> pos_40;
	vector <double> pos_50;


	//vector <double> varhet_hist;

	//data files
	ofstream flog, fpop, fhet, fprof, frough,fpos;
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
	Kstr << K;
	Mstr << M;
	Bstr << B;
	Gstr << g0;
	string param_string =  "K"+Kstr.str()+"_M" + Mstr.str() + "_B" +Bstr.str() + "_G" +Gstr.str() +"_";



	string logName = "log_" + param_string + date_time.str() + ".txt";
	string hetName = "het_" + param_string +  date_time.str() + ".txt";
	//string varhetName = "varhet_" + param_string +  date_time.str() + ".txt";
	string popName = "pop_"+ param_string +  date_time.str() + ".txt";
	string profName = "prof_" + param_string + date_time.str() + ".txt";
	string roughName = "rough_" + param_string + date_time.str() + ".txt";
	string posName = "pos_" + param_string + date_time.str() + ".txt";
	string folder = "";

    flog.open(folder+logName);
    fhet.open(folder+hetName);
    fpop.open(folder+popName);
    fprof.open(folder+profName);
    frough.open(folder+roughName);
    fpos.open(folder+posName);
    //fvarhet.open(varhetName);








	for(int i = 0; i < int(n_demesh*.5); i++){
		for(int j = 0; j < int(n_demesw/2); j++){

		deme[i][j][0] = .5*K;

		}

	}

	for(int i = 0; i < int(n_demesh*.5); i++){
		for(int j = int(n_demesw/2); j < n_demesw; j++){

		deme[i][j][1] = K;

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
	//while(ht>h_thresh){
		int d_start = 0;
		int d_end= n_demesh ;
		int empty_found=0;


		for(int ii = 0; ii < int(n_demesh); ii++){
			int full = 0;
			int empty = 0;
			for(int jj = 0; jj < int(n_demesw); jj++){

				deme_aux[ii][jj][0] = deme[ii][jj][0];
				deme_aux[ii][jj][1] = deme[ii][jj][1];
				if(deme[ii][jj][0] + deme[ii][jj][1] == K){
					full+=1;
				}	

				if(deme[ii][jj][0] + deme[ii][jj][1] == 0){
					empty+=1;
				}

			}
			if ((full ==n_demesw)&&(freeze_flag == 1)){
				d_start = ii;
			}
			if ((empty ==n_demesw)&&(empty_found==0)&&(fast_samp_flag == 1)){
				d_end = ii;
				empty_found=1;
			}

		}

		if (d_end <n_demesh-10){
			d_end =d_end+10;
		}
		//cout<<d_end<<endl;
		//int i=0;
		for( int j=0; j<n_demesw; j++){
			int arr[2] = {d_start, j};
			int neighb_vec[3][2] = {{0,1},{1,0},{0,-1}};
			int neighbs[3][2];
			//int pop_sum = fast_samp_flag;

			for(int ne=0; ne <3; ne++){
				neighbs[ne][0] = (arr[0] + n_demesh+neighb_vec[ne][0]) % n_demesh;
				neighbs[ne][1] = (arr[1] + n_demesw+neighb_vec[ne][1]) % n_demesw;
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
			deme[d_start][j][0] = new_cnt[1];
			deme[d_start][j][1] = new_cnt[2];
		}

		for(int i = d_start+1; i < d_end-1 ; i++){
			for(int j = 0; j < n_demesw; j++){
				int arr[2] = {i, j};
				int neighb_vec[4][2] = {{0,1},{0,-1},{1,0},{-1,0}};
				int neighbs[4][2];
				int pop_sum = 1;

				for(int ne=0; ne <4; ne++){
					neighbs[ne][0] = (arr[0] + n_demesh+neighb_vec[ne][0]) % n_demesh;
					neighbs[ne][1] = (arr[1] + n_demesw+neighb_vec[ne][1]) % n_demesw;
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

		int i=n_demesh-1;
		for( int j=0; j<d_end; j++){
			int arr[2] = {i, j};
			int neighb_vec[3][2] = {{0,1},{-1,0},{0,-1}};
			int neighbs[3][2];
			//int pop_sum = fast_samp_flag;
			//cout << i << ", " << j << endl;

			for(int ne=0; ne <3; ne++){
				neighbs[ne][0] = (arr[0] + n_demesh+neighb_vec[ne][0]) % n_demesh;
				neighbs[ne][1] = (arr[1] + n_demesw+neighb_vec[ne][1]) % n_demesw;
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








		if (sumDeme(deme,n_demesh) > .5*n_demesh*n_demesw){
			int shift =  int(sumDeme(deme,n_demesh) - .5*n_demesh*n_demesw)+1;
			//cout << shift << endl;
			if ((shift< 0) == true){

	            cout << "Negative shift of array!" << endl;
	            exit(EXIT_FAILURE);

	        }


			for (int i = 0; i < n_demesh - shift; i++){
				for (int j = 0; j < n_demesw; j++){ 

					for (int k = 0; k <n_spec;k++){

						deme[i][j][k] = deme[i+shift][j][k];


					}
				}
			}


	        for (int i = int(n_demesh - shift); i < n_demesh; i++){
	        	for (int j = 0; j < n_demesw; j++){ 

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
			//varhet_hist.push_back(calcVarHet(deme, n_demesh)); 
			het_hist.push_back(calcHet(deme, n_demesh)); // Store heterozygosity
	        pop_hist.push_back(pop_shift+sumDeme(deme,n_demesh));
	        rough_0.push_back(calcRoughness(deme[0]));
	        rough_10.push_back(calcRoughness(deme[10]));
	        rough_20.push_back(calcRoughness(deme[20]));
	        rough_30.push_back(calcRoughness(deme[30]));
	        rough_40.push_back(calcRoughness(deme[40]));
	        rough_50.push_back(calcRoughness(deme[50]));
	        pos_0.push_back(calcBoundPos(deme[0]));
	        pos_10.push_back(calcBoundPos(deme[10]));
	        pos_20.push_back(calcBoundPos(deme[20]));
	        pos_30.push_back(calcBoundPos(deme[30]));
	        pos_40.push_back(calcBoundPos(deme[40]));
	        pos_50.push_back(calcBoundPos(deme[50]));
	        int lastRow = calcLastRow(deme);



	        //cout <<calcLastRow(deme)<<endl;


	        //ht= calcHet(deme, n_demesh);

	        /*if ((het_hist[het_hist.size()-2]> prof_count*.1) && (het_hist[het_hist.size()-1]< prof_count*.1)  ){
	        	ostringstream strh;
	        	strh << prof_count*.1 ;

	        	string proftName = "prof_H"+ strh.str() + "_" +param_string + date_time.str() + ".txt";
	        	ofstream fproft;
	            fproft.open("sim_data/"+proftName);
	            for(int i = 0; i <n_demesh; i++){
	            	for(int j = 0; j <n_demesw; j++){
	            		fproft << i << ", " << j << ", " << deme[i][j][0] << ", " << deme[i][j][1] <<endl;


	            	}
            	}
            	prof_count-=1;

	        }*/

	        if (prof_hist !=0){
	        	ostringstream strT;
	        	strT << dt;
	        	string proftName = "prof_T"+ strT.str() + "_" + date_time.str() + ".txt";
	        	ofstream fproft;
	            fproft.open(proftName);
	            for(int i = 0; i <n_demesh; i++){
	            	for(int j = 0; j <n_demesw; j++){
	            		fproft << i << ", " << j << ", " << deme[i][j][0] << ", " << deme[i][j][1] <<endl;


	            	}
            	}

	        }


		}

		//dt+=1;


    }
   //int n_data = int(dt/record_time);
   //frough <<'time' <<  "0"  << ", " <<"10" << ", "<< "20" << ", " << "30"<< ", " << "40"<< ", " << "50" << ", "<<"front_8"<< ", "<<"front_6"<< ", "<<"front_4"<< ", "<<"front_2"<< ", "<<"front"<< ", "<<"front_deme" <<endl;

   for(int i=0;i < n_data;i++){
   		//fvarhet << int(i*record_time) << ", "  << varhet_hist[i] << endl;
    	fhet << int(i*record_time) << ", "  << het_hist[i] << endl;
    	fpop << int(i*record_time) << ", "  << pop_hist[i] << endl;
    	frough << int(i*record_time) << ", " << rough_0[i] << ", "<< rough_10[i] << ", "<< rough_20[i]<< ", " << rough_30[i]<< ", " << rough_40[i]<< ", " << rough_50[i]<<endl;
    	fpos << int(i*record_time) << ", " << pos_0[i] << ", "<< pos_10[i] << ", "<< pos_20[i]<< ", " << pos_30[i]<< ", " << pos_40[i]<< ", " << pos_50[i]<<endl;
    }

    for(int i=0;i < n_demesh; i++){
    	for(int j =0;j < n_demesw; j++){

    		fprof << i << ", " << j<< ", " << deme[i][j][0] << ", "  << deme[i][j][1] << endl;

		}

    }

    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes height,Number of demes width, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demesh<<n_demesw << time_start<< run_time<< endl;

    fhet.close();
    fpop.close();
    flog.close();
    fprof.close();
    frough.close();
    fpos.close();
    //fvarhet.close();

    cout << "Finished!" << "\n";
    cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);









	return 0;





}


