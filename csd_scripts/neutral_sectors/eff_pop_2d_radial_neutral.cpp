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


unsigned long K  = 50 ;
unsigned int n_gens = 650;
const int n_demes = 1000;
const unsigned int n_spec = 2;
float M = 0.1;
float B =  3;
float g0 = 0.5;
int initMut = 25;
int initRad = 10;
unsigned long prof_hist = 0;
unsigned long fast_samp_flag = 1;
unsigned int ID_seed = 2;



double sumDeme(int arr[][n_demes][n_spec], int arrSize){
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
int mutSecFlag;



int checkMax(std::vector<std::vector<std::vector<long double> > > arr, const int arrSize){
	float distMax = 0;

	for(int i = 0; i < arrSize; i++){
		for(int j=0; j<arrSize;j++){
			if ((arr[i][j][0]+arr[i][j][1]) > 0){
				float dist =  pow(pow(abs(arrSize/2 - i),2) + pow(abs(arrSize/2 - j),2),.5);

				if (dist>distMax){

					distMax = dist;
				}


			}



		}
	}






	return distMax;
}


int checkEmpty(std::vector<std::vector<std::vector<long double> > > arr, const int arrSize) {
	int buff = 2;
	int check_pop = 0;
	int emptyBounds = 1;

	for(int i = 0; i < arrSize; i++){


		check_pop+= arr[i][buff][0]+arr[i][buff][1];
		check_pop+=  arr[buff][i][0]+arr[buff][i][1];
		check_pop+=  arr[i][arrSize-1-buff][0]+arr[i][arrSize-1-buff][1];
		check_pop+=  arr[arrSize-1-buff][i][0]+arr[arrSize-1-buff][i][1];

	}
	if (check_pop>0){
		emptyBounds = 0;
	}



	return emptyBounds;


}

int neighborMatch(std::vector<std::vector<std::vector<long double> > > arrNum, int (&grid)[n_demes][n_demes], int x, int y){
	int mutSecFlag = 0;
	//float thresh =.95;
	grid[x][y] =1;


	//std::cout<<"hi"<<std::endl;


	if (grid[x-1][y] == 0){
		mutSecFlag+=neighborMatch(arrNum, grid, x-1,y);
		grid[x-1][y]=1;
		//mutBubFlag +=1;
	}
	if (grid[x+1][y] ==0){
		mutSecFlag+=neighborMatch(arrNum, grid, x+1,y);
		grid[x+1][y]=1;
		//utBubFlag +=1;
	}
	if (grid[x][y-1]==0){
		mutSecFlag+=neighborMatch(arrNum, grid, x,y-1);
		grid[x][y-1]=1;
		//mutBubFlag +=1;
	}
	if (grid[x][y+1] == 0){
		mutSecFlag+=neighborMatch(arrNum, grid, x,y+1);
		grid[x][y+1]=1;
		//mutBubFlag +=1;
	}

	if ( (grid[x][y+1] ==-2) || (grid[x][y-1] ==-2) || (grid[x+1][y] ==-2) || (grid[x-1][y] ==-2) ){
		mutSecFlag+= 1;


	}

	//if ((mutSecFlag>0)){
	//	mutSecFlag=1;

	//}
	//else
	//	mutSecFlag=0;


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

//long double deme[n_demes][n_demes][n_spec] = {{0}};
//long double deme_aux[n_demes][n_demes][n_spec] = {{0}};

//int * deme = new int[n_demes][n_demes][n_spec];
//int * deme_aux = new int[n_demes][n_demes][n_spec];




std::vector<std::vector<std::vector<long double> > > deme(n_demes, std::vector<std::vector<long double> >(n_demes, std::vector<long double>(n_spec)));
std::vector<std::vector<std::vector<long double> > > deme_aux(n_demes, std::vector<std::vector<long double> >(n_demes, std::vector<long double>(n_spec)));




int main (int argc, char * argv[]){
	using namespace std;
	//cout <<"Hi"<< endl;

	int c;
    while ((c = getopt (argc, argv, "T:B:I:U:G:M")) != -1)
    {
        if (c == 'T')
            n_gens  = atoi(optarg); // carrying capacity
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        else if (c == 'I')
            initMut = atoi(optarg); // migration probability
        else if (c == 'U')
            initRad = atoi(optarg); // migration probability
        else if (c == 'G')
            ID_seed  = atoi(optarg); // growth rate
        else if (c == 'M')
            M = atof(optarg); // growth rate


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
	gsl_rng_set(r, ID_seed);

	double new_prob[n_spec + 1];
	unsigned int new_cnt[n_spec + 1];
	//int n_data = 10;
	int record_time = 5;
	//int record_time = 10;
	//int n_data = int(n_gens/record_time);

	double pop_shift = 0.0;
	double w_s1 = 1.0;
	double w_s2 = 1.0;
	double w_avg;
	double w_v;
	double mutThresh =0.9;
	int minSize =10;
	//vector <double> pop_hist;
	//vector <double> het_hist;
	vector <double> sect_hist;
	vector <double> max_hist;
	//vector <double> mut_hist;
	//vector <double> full_hist;
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

	ostringstream date_time, Kstr,  Mstr, Bstr, Gstr, Ustr, Istr;
	date_time << buffer;
	Kstr << K;
	Mstr << M;
	Bstr << B;
	Ustr << initRad;
	Istr << ID_seed;
	Gstr << g0;
	string param_string =  "K"+Kstr.str()+"_M" + Mstr.str() + "_B" +Bstr.str() + "_G" +Gstr.str() +"_U"+Ustr.str()+"_I"+Istr.str()+"_";



	string logName = "log_" + param_string + date_time.str() + ".txt";
	string hetName = "het_" + param_string +  date_time.str() + ".txt";
	string varhetName = "varhet_" + param_string +  date_time.str() + ".txt";
	string popName = "pop_"+ param_string +  date_time.str() + ".txt";
	string profName = "prof_" + param_string + date_time.str() + ".txt";
	string sectName = "sect_" + param_string + date_time.str() + ".txt";
	//string mutName = "mut_" + param_string + date_time.str() + ".txt";
	//string fullName = "full_" + param_string + date_time.str() + ".txt";
	//string folder = "sim_data/";
	string folder = "";


    flog.open(folder+logName);
    //fhet.open(folder+hetName);
    //fpop.open(folder+popName);
    fprof.open(folder + profName);
    //fvarhet.open(varhetName);
   	fsect.open(folder + sectName);
    //fsect.open(folder + sectName)
    //fsect.open(folder + "sector_results.txt" , ios_base::app);
    //fmut.open(folder + mutName);
    //ffull.open(folder + fullName);
    //cout <<"Hi"<< endl;








	for(int i = 0; i < int(n_demes); i++){
		for(int j = 0; j < int(n_demes); j++){
			if ( round(sqrt( abs(i-int(n_demes*.5))*abs(i-int(n_demes*.5)) + abs(j-int(n_demes*.5))*abs(j-int(n_demes*.5))) ) < initRad)
			{
				deme[i][j][1] = initMut;
				deme[i][j][0] = K - initMut;


			}

		}


	}
	//initial population in middle
	//deme[int(n_demes/2)][int(n_demes/2)][0] = K;
	//deme[int(n_demes/2)][int(n_demes/2)][1] = K;
	int dt =0;





	while(checkEmpty(deme, n_demes) == 1){
	//for (int dt = 0 ; dt < n_gens; dt++ ){

		
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
			//cout << checkMax(deme,n_demes) << endl;

			int arrGrid[n_demes][n_demes] = {{0}};
			for (int i=0; i<n_demes; i++){
				for (int j=0; j<n_demes; j++){
					if (deme[i][j][0] +deme[i][j][1] ==0){
						arrGrid[i][j] = -2;

					}
					if ((deme[i][j][1] / (deme[i][j][0] + deme[i][j][1])) <mutThresh){
						arrGrid[i][j] = -1;

					}

				}

			}

			int sectCounts = 0;
			int rawSec;
			for(int i =1; i<n_demes-1;i++){
				for(int j = 1; j< n_demes-1;j++){

					if (arrGrid[i][j] == 0){
						rawSec = neighborMatch(deme, arrGrid, i,j );

						int mutBubFlag = 0;
						if( rawSec> minSize){
							sectCounts+= 1;

						}

					}

				}

			}

			
			/*int fullDeme=0;
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
			//het_hist.push_back(calcHet(deme, n_demes)); 
	        //pop_hist.push_back(pop_shift+sumDeme(deme,n_demes));
	        sect_hist.push_back(sectCounts);
	        max_hist.push_back( checkMax(deme,n_demes) );

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
	

		dt+=1;

    }
	int arrGrid[n_demes][n_demes] = {{0}};
	for (int i=0; i<n_demes; i++){
		for (int j=0; j<n_demes; j++){
			if (deme[i][j][0] +deme[i][j][1] ==0){
				arrGrid[i][j] = -2;

			}
			if ((deme[i][j][1] / (deme[i][j][0] + deme[i][j][1])) <mutThresh){
				arrGrid[i][j] = -1;

			}

		}

	}

	//ofstream fgrid;
    //fgrid.open("finalgrid.txt");
    //for(int i = 0; i <n_demes; i++){
    //	for(int j = 0; j <n_demes; j++){
   // 		fgrid << i << " " << j << " " << arrGrid[i][j] <<endl;


    //	}
	//}

	int sectCounts = 0;
	int rawSec;
	for(int i =1; i<n_demes-1;i++){
		for(int j = 1; j< n_demes-1;j++){

			if (arrGrid[i][j] == 0){
				rawSec = neighborMatch(deme, arrGrid, i,j );

				int mutBubFlag = 0;
				if (rawSec>0){
					//cout <<rawSec<<endl;

				}
				
				if( rawSec> minSize){
					//cout <<rawSec<<endl;
					sectCounts+= 1;

				}

			}

		}

	}
	//cout <<sectCounts<<endl;
	//sect_hist.push_back(sectCounts);
	//fsect << B << " " << initMut << " " << ID_seed << " "<< sectCounts<<endl; 
	//cout <<"Hi"<<endl;
	for(int i=0;i < sect_hist.size();i++){

    	//fvarhet << int(i*record_time) << ", "  << varhet_hist[i] << endl;
    	//fhet << int(i*record_time) << " "  << het_hist[i] << endl;
    	//fpop << int(i*record_time) << " "  << pop_hist[i] << endl;
    	fsect << max_hist[i] << " "  << sect_hist[i] << endl;
    	//fmut << int(i*record_time) << ", "  << mut_hist[i] << endl;
    	//ffull << int(i*record_time) << ", "  << full_hist[i] << endl;

    }

    for(int i=0;i < n_demes; i++){
    	for(int j =0;j < n_demes; j++){

    		fprof << i << " " << j<< " " << deme[i][j][0] << " "  << deme[i][j][1] << endl;

		}

    }

   

	

	//float mutProp = mutDeme/fullDeme;
	//fsect << sectCounts << ", " << mutDeme << ", " <<fullDeme  << endl;*/

    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;



    flog << "Number of generations, Number of species, Growth rate, Migration rate, B, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << n_gens << ", " <<  n_spec << ", " << g0 << ", " << M << ", " << n_demes << time_start<< run_time<< endl;

    //fhet.close();
    //fpop.close();
    flog.close();
    fprof.close();
    fsect.close();
    //fmut.close();
    //ffull.close();
    //fvarhet.close();



    cout << "Finished!" << "\n";
    cout << "B: " <<B <<" init mut. "<< initMut <<  " seed: "<< ID_seed << endl;
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);

	








	return 0;
	


}