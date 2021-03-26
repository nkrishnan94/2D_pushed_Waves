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


unsigned long K  = 500;
unsigned int n_gens = 2000;
const int n_demes = 500;
const unsigned int n_spec = 2;
float M = 0.25;
float B = 5;
float g0 = 0.01;
int initMut = 250;
unsigned long prof_hist = 0;
unsigned long fast_samp_flag = 1;
double j_mu = 20;
int d_migs = 5;



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

int mutBubFlag = 0;


int neighborMatch(long double arrNum[][n_demes][n_spec], long double arrGrid[][n_demes], int x, int y){
	int mutSecFlag = 0;


	//std::cout<<"hi"<<std::endl;
	if ((arrGrid[x][y]==-1) && ((arrNum[x][y][1]/ (arrNum[x][y][0]+arrNum[x][y][1]))  > .8 ) ){
		arrGrid[x][y]=0;

		if (arrNum[x-1][y][1]/ (arrNum[x-1][y][0]+arrNum[x-1][y][1]) >.8){
			mutSecFlag+=neighborMatch(arrNum, arrGrid, x-1,y);
			arrGrid[x-1][y]=0;
			mutBubFlag +=1;
		}
		if (arrNum[x+1][y][1]/ (arrNum[x+1][y][0]+arrNum[x+1][y][1]) >.8){
			mutSecFlag+=neighborMatch(arrNum, arrGrid, x+1,y);
			arrGrid[x+1][y]=0;
			mutBubFlag +=1;
		}
		if (arrNum[x][y-1][1]/ (arrNum[x][y-1][0]+arrNum[x][y-1][1]) >.8){
			mutSecFlag+=neighborMatch(arrNum, arrGrid, x,y-1);
			arrGrid[x][y-1]=0;
			mutBubFlag +=1;
		}
		if (arrNum[x][y+1][1]/ (arrNum[x][y+1][0]+arrNum[x][y+1][1]) >.8){
			mutSecFlag+=neighborMatch(arrNum, arrGrid, x,y+1);
			arrGrid[x][y+1]=0;
			mutBubFlag +=1;
		}
		if ((arrNum[x][y+1][0]+arrNum[x][y+1][1] ==0) || (arrNum[x][y-1][0]+arrNum[x][y-1][1] ==0) || (arrNum[x+1][y][0]+arrNum[x+1][y][1] ==0) ||(arrNum[x-1][y][0]+arrNum[x-1][y][1] ==0)){
			mutSecFlag+= 1;

		}

	}

	if ((mutSecFlag>0)){
		mutSecFlag=1;

	}
	else
		mutSecFlag=0;


	return mutSecFlag;


}







float calcHet(long double arr[][n_demes][n_spec], const int arrSize){



	int cnt =  0 ;
	long double H = 0.0;

	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){

			double deme_pop = arr[i][j][0]+arr[i][j][1];


			//std::cout << i << "\n";
			if (deme_pop > 0.0){
				H += (2*arr[i][j][0]*(deme_pop - arr[i][j][0]))/(deme_pop*deme_pop);
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
    while ((c = getopt (argc, argv, "K:Z:B:M:U:G:F")) != -1)
    {
        if (c == 'K')
            K  = atoi(optarg); // carrying capacity
        else if (c == 'Z')
            prof_hist = atoi(optarg); //keep track of profile through time 
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        else if (c == 'M')
            M = atof(optarg); // migration probability
        else if (c == 'U')
            initMut = atoi(optarg); // migration probability
        else if (c == 'G')
            g0 = atof(optarg); // growth rate
        else if (c == 'F')
            fast_samp_flag = atoi(optarg); // growth rate

    }
    /*if ((B >= 2))
    n_gens = 1*K;


	if ((B < 2))
	    n_gens = 15*int(sqrt(K));*/


	const gsl_rng_type * T;
	gsl_rng * r;



	gsl_rng_env_setup();
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	//int sysRandom;
	gsl_rng_set(r, time(NULL));

	double new_prob[n_spec + 1];
	unsigned int new_cnt[n_spec + 1];
	int n_data = 30;
	int record_time = int(n_gens/n_data);
	//int record_time = 10;
	//int n_data = int(n_gens/record_time);

	double pop_shift = 0.0;
	double w_s1 = 1.0;
	double w_s2 = 1.0;
	double w_avg;
	double w_v;
	vector <double> pop_hist;
	vector <double> het_hist;
	vector <double> sect_hist;
	vector <double> mut_hist;
	vector <double> full_hist;
	//vector <double> varhet_hist;

	//data files
	ofstream flog, fpop, fhet, fprof, fsect,fmut,ffull;
	time_t time_start;
	clock_t c_init = clock();
	struct tm * timeinfo;
	char buffer [80];


    time (&time_start);
	timeinfo = localtime (&time_start);


	strftime (buffer,80,"%F-%H-%M-%S",timeinfo);

	ostringstream date_time, Kstr,  Mstr, Bstr, Gstr, Istr, Jstr, JNstr;
	date_time << buffer;
	Kstr << K;
	Mstr << M;
	Bstr << B;
	Istr << initMut;
	Gstr << g0;
	Jstr << j_mu;
	JNstr << d_migs;
	string param_string =  "K"+Kstr.str()+"_M" + Mstr.str() + "_B" +Bstr.str() + "_G" +Gstr.str() +"_I"+Istr.str()+"_J"+Jstr.str()+"_JN"+JNstr.str()+"_";



	string logName = "log_" + param_string + date_time.str() + ".txt";
	string hetName = "het_" + param_string +  date_time.str() + ".txt";
	string varhetName = "varhet_" + param_string +  date_time.str() + ".txt";
	string popName = "pop_"+ param_string +  date_time.str() + ".txt";
	string profName = "prof_" + param_string + date_time.str() + ".txt";
	string sectName = "sect_" + param_string + date_time.str() + ".txt";
	string mutName = "mut_" + param_string + date_time.str() + ".txt";
	string fullName = "full_" + param_string + date_time.str() + ".txt";
	//string folder = "sim_data_neutral/";
	string folder = "";


    flog.open(folder+logName);
    fhet.open(folder+hetName);
    fpop.open(folder+popName);
    fprof.open(folder + profName);
    //fvarhet.open(varhetName);
    //fsect.open(folder + sectName);
    //fmut.open(folder + mutName);
    //ffull.open(folder + fullName);








	for(int i = int(n_demes*.5)-50; i < int(n_demes*.5)+50; i++){
		for(int j = int(n_demes*.5)-50; j < int(n_demes*.5)+50; j++){
			if (sqrt(abs(i-int(n_demes*.5))*abs(i-int(n_demes*.5)) + abs(j-int(n_demes*.5))*abs(j-int(n_demes*.5))) < 50)
			{
				deme[i][j][1] = initMut;
				deme[i][j][0] = K - initMut;


			}

		}


	}
	//initial population in middle
	//deme[int(n_demes/2)][int(n_demes/2)][0] = K;
	//deme[int(n_demes/2)][int(n_demes/2)][1] = K;






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
				int neighbs[4][2];
				int pop_sum = fast_samp_flag;

				for(int ne=0; ne <4; ne++){
					neighbs[ne][0] = (arr[0] + n_demes+neighb_vec[ne][0]) % n_demes;
					neighbs[ne][1] = (arr[1] + n_demes+neighb_vec[ne][1]) % n_demes;
					pop_sum += deme[i][j][neighbs[ne][0]] + deme[i][j][neighbs[ne][1]]+deme_aux[i][j][neighbs[ne][0]]+deme_aux[i][j][neighbs[ne][1]];




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
					w_avg = w_v + (w_s1 - w_v)*f1 +(w_s2 - w_v)*f2;
					if ((deme[i][j][0]+deme[i][j][1]) < K){

						f1 *= w_s1/w_avg;
						f2 *= w_s2/w_avg;
						new_prob[0] = 1-f1-f2;
						new_prob[1] = f1;
						new_prob[2] = f2;
						gsl_ran_multinomial(r,n_spec+1,K,new_prob,new_cnt);
						deme[i][j][0] = new_cnt[1];
						deme[i][j][1] = new_cnt[2];
					}


				}


			}
		}

		for(int dm=0; dm<d_migs; dm++ ){
			int sign1;
			double rd1 = gsl_ran_flat(r, 0, 1);
			int sign2;
			double rd2 = gsl_ran_flat(r, 0, 1);
			if(rd1<.5){
				sign1 =-1;

			}else{
				sign1 =1;
			}

			if(rd2<.5){
				sign2 =-1;

			}else{
				sign2 =1;
			}


			unsigned int j1 = gsl_ran_poisson(r, j_mu);
			unsigned int j2 = gsl_ran_poisson(r, j_mu);
			int d1 =  gsl_ran_flat(r, 0, n_demes);
			int d2 =  gsl_ran_flat(r, 0, n_demes);


			if ((d2+sign2*j2)<n_demes && (d2+sign2*j2)>0 && (d1+sign1*j1)<n_demes && (d1+sign1*j1)>0){

				long double d1a1 = deme[d1][d2][0];
				long double d1a2 = deme[d1][d2][1];
				long double d2a1 = deme[d1+sign1*j1][d2+sign2*j2][0];
				long double d2a2 = deme[d1+sign1*j1][d2+sign2*j2][1];
				deme[d1][d2][0] =d2a1;
				deme[d1][d2][1] =d2a2;
				deme[d1+sign1*j1][d2+sign2*j2][0] = d1a1;
				deme[d1+sign1*j1][d2+sign2*j2][1] = d1a2;




			}
		}







		

		




		/*if (sumDeme(deme,n_demes)/(K*n_demes) > .5*n_demes*n_demes){
			int shift =  int(sumDeme(deme,n_demes)/K - .5*n_demes)+1;
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





		}*/


		//cout << dt << endl;
		//cout << record_time << endl;
		if (dt % record_time == 0){

			long double arrGrid[n_demes][n_demes] = {{0}};
			for (int i=0; i<n_demes; i++){
				for (int j=0; j<n_demes; j++){
					arrGrid[i][j] = -1;

				}

			}

			/*int sectCounts = 0;
			for(int i =0; i<n_demes;i++){
				for(int j = 0; j< n_demes;j++){
					if (arrGrid[i][j] == -1){
						int mutBubFlag = 0;
						sectCounts+= neighborMatch(deme, arrGrid, i,j );
						//cout <<sectCounts<<endl;

					}

				}

			}

			
			int fullDeme=0;
			int mutDeme=0;
			for(int i =0; i<n_demes;i++){
				for(int j = 0; j< n_demes;j++){
					if((deme[i][j][0]+deme[i][j][1]) >0){
						fullDeme+=1;
					}
					if(deme[i][j][1] >.5){
						mutDeme+=1;
					}
				}
			}*/
			//varhet_hist.push_back(calcVarHet(deme, n_demes)); 
			het_hist.push_back(calcHet(deme, n_demes)); 
	        pop_hist.push_back(pop_shift+sumDeme(deme,n_demes));
	        //sect_hist.push_back(sectCounts);
	        //mut_hist.push_back(mutDeme);
	        //full_hist.push_back(fullDeme);

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

    	//fvarhet << int(i*record_time) << ", "  << varhet_hist[i] << endl;
    	fhet << int(i*record_time) << ", "  << het_hist[i] << endl;
    	fpop << int(i*record_time) << ", "  << pop_hist[i] << endl;
    	//fsect << int(i*record_time) << ", "  << sect_hist[i] << endl;
    	//fmut << int(i*record_time) << ", "  << mut_hist[i] << endl;
    	//ffull << int(i*record_time) << ", "  << full_hist[i] << endl;

    }

    for(int i=0;i < n_demes; i++){
    	for(int j =0;j < n_demes; j++){

    		fprof << i << ", " << j<< ", " << deme[i][j][0] << ", "  << deme[i][j][1] << endl;

		}

    }

   

	

	//float mutProp = mutDeme/fullDeme;
	//fsect << sectCounts << ", " << mutDeme << ", " <<fullDeme  << endl;

    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demes << time_start<< run_time<< endl;

    fhet.close();
    fpop.close();
    flog.close();
    fprof.close();
    //fsect.close();
    //fmut.close();
    //ffull.close();
    //fvarhet.close();

    cout << "Finished!" << "\n";
    //cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);

	








	return 0;
	


}
