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
unsigned int n_gens = 100000;
//float h_thresh = .1;
const int n_demesh = 120; 
const int n_demesw =1100; 
const unsigned int n_spec = 2;
float M = 0.25;
float B =0;
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

	return sum/(K*n_demesw);




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


int getFrontDeme(long double arr[n_demesh][n_demesw][n_spec], int L){

	long double N = n_demesw;
	long double frontSum = 0; 
	long double frontDiffs=0;
	long double fronts[n_demesw];

	for(int j =int(n_demesw/2) -int(L/2) ; j< int(n_demesw/2)+int(L/2); j++){

		int frontFound=0;
		for(int i = 1; i < n_demesh; i++){
			if( (arr[i][j][0]+arr[i][j][1]==0) && (arr[i-1][j][0]+arr[i-1][j][1]>0) &&(frontFound==0)){
				fronts[j]=(i-1);
				frontSum+=(i-1);
				frontFound=1;

			}


		}
	}

	for(int j =int(n_demesw/2) -int(L/2) ; j< int(n_demesw/2)+int(L/2); j++){
		frontDiffs+=  pow((fronts[j] - (frontSum/L)  ),2);

	}


	//long double front_rough = frontDiffs/N;
	
	return frontDiffs/L;

}

int getEdgeDeme(long double arr[n_demesh][n_demesw][n_spec], int L){

	long double N = n_demesw;
	 
	long double frontDiffs=0;
	long double fronts[L];

	for(int j =0 ; j< L; j++){

		int frontFound=0;
		for(int i = 1; i < n_demesh; i++){
			if( (arr[i][int(n_demesw/2)- int(L/2)+j][0]+arr[i][int(n_demesw/2)- int(L/2)+j][1]>0) &&(frontFound==0)){
				fronts[j]=i;


			}


		}
	}
	long double frontSum = 0; 
	for(int j =0 ; j< L; j++){
		frontSum+= fronts[j];

	}

	for(int j =0 ; j< L; j++){
		frontDiffs+=  pow((fronts[j] - (frontSum/L)  ),2);

	}



	//long double front_rough = frontDiffs/N;
	
	return frontDiffs;

}

int countSectors(long double arr[n_demesh][n_demesw][n_spec]){
	long double front[n_demesw];

	for(int i = calcLastRow(arr); i < n_demesh; i++){

		for(int j=0;j<n_demesw;j++){

			if ((arr[i][j][0]+arr[i][j][1])>0){
				front[j] = arr[i][j][0]/( arr[i][j][0]+ arr[i][j][1]);

			}
		}


	}

	int sectors = 0;

	for(int i=1; i<n_demesw-1;i++){

		if (((front[i]==0)|| (front[i]==1)) && (front[i-1]!=front[i])&& (front[i+1]!=front[i])){

			sectors+=1;


		}


	}

	return sectors+1;
	
}




float calcVarHet(long double arr[n_demesh][n_demesw][n_spec], const int arrSize){

	long double hets[n_demesh][n_demesw];
	long double H = 0.0;
	long double varH = 0.0;
	int cnt=0;
	float average;


	for(int i = 0; i < n_demesh; i++){
		for(int j=0; j<n_demesw;j++){

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
    while ((c = getopt (argc, argv, "K:Z:B:T:M:G:F")) != -1)
    {
        if (c == 'K')
            K  = atoi(optarg); // carrying capacity
        else if (c == 'Z')
            prof_hist = atoi(optarg); //keep track of profile through time 
        else if (c == 'B')
            B = atof(optarg); // cooperativity
        //else if (c == 'T')
        //    h_thresh = atof(optarg); // cooperativity
        else if (c == 'M')
            M = atof(optarg); // migration probability
        else if (c == 'G')
            g0 = atof(optarg); // growth rate
        else if (c == 'F')
            fast_samp_flag = atoi(optarg); // growth rate

    }
    

    //n_gens = (B+1)*K;


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
	//int dt = 0;

	double pop_shift = 0.0;
	double w_s = 1.0;
	double w_avg;
	double w_v;

	vector <double> pop_hist;
	vector <double> het_hist;
	vector <double> sects_hist;
	vector <double> rough_hist_16;
	vector <double> rough_hist_32;
	vector <double> rough_hist_64;
	vector <double> rough_hist_128;
	vector <double> rough_hist_256;
	vector <double> rough_hist_512;
	vector <double> rough_hist_1024;
	
	//vector <double> varhet_hist;

	//data files
	ofstream flog, fpop, fhet, fprof,fsects,frough_16,frough_32,frough_64,frough_128,frough_256,frough_512,frough_1024;
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
	string sectsName = "sects_" + param_string + date_time.str() + ".txt";
	string rough10Name = "rough_16_" + param_string + date_time.str() + ".txt";
	string rough30Name = "rough_32_" + param_string + date_time.str() + ".txt";
	string rough50Name = "rough_64_" + param_string + date_time.str() + ".txt";
	string rough70Name = "rough_128_" + param_string + date_time.str() + ".txt";
	string rough90Name = "rough_256_" + param_string + date_time.str() + ".txt";
	string rough110Name = "rough_512_" + param_string + date_time.str() + ".txt";
	string rough130Name = "rough_1024_" + param_string + date_time.str() + ".txt";

	//string folder = "KPZ/sim_data_KPZ/";
	string folder = "";

    flog.open(folder+logName);
    fhet.open(folder+hetName);
    fpop.open(folder+popName);
    fprof.open(folder+profName);
    //fsects.open(folder+sectsName);
    frough_16.open(folder+rough10Name);
    frough_32.open(folder+rough30Name);
    frough_64.open(folder+rough50Name);
    frough_128.open(folder+rough70Name);
    frough_256.open(folder+rough90Name);
    frough_512.open(folder+rough110Name);
    frough_1024.open(folder+rough130Name);

    //fvarhet.open("sim_data/"+varhetName);





    cout<< "hi";
    string line;
	ifstream myfile ("KPZ/fisher_wave_files/K"+Kstr.str()+"_"+"B"+Bstr.str()+"_"+"fisher.txt");
	//cout<< "KPZ/fisher_wave_files/K"+Kstr.str()+"_"+"B"+Bstr.str()+"_"+"fisher.txt";
	

	int j = 0;
	if (myfile.is_open())
  	{


		while ( getline (myfile,line) )
	    {
	      string::iterator it;
	      int index = 0;
	      //for ( it = line.begin() ; it < line.end(); it++ ,index++)
	      //{
	       // cout << *it;
	        //cout << line << '\n';
	      //}
	      deme[j][0][0] = stof(line);
	      //cout << stoi(line);
	      deme[j][0][1] = 0;
	      j+=1;
	      //cout << line[1,2,3];

	   
	    }
	    myfile.close();


	}
	


	for(int i = 1; i < int(n_demesh); i++){
		for(int j = 0; j < n_demesw; j++){

			deme[i][j][0] = deme[i][0][0];
			//deme[i][j][1] = .5*K;


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
			int full =0;
			int empty=0;
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

		for(int i = d_start; i < d_end -1; i++){
			for(int j = 0; j < n_demesw; j++){
				int arr[2] = {i, j};
				int neighb_vec[4][2] = {{0,1},{0,-1},{1,0},{-1,0}};
				int neighbs[4][2];
				//int pop_sum = fast_samp_flag;

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

		int i=d_end;
		for( int j=0; j<n_demesw; j++){
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

		

		




		if (sumDeme(deme,n_demesh) > .5*n_demesh){
			int shift =  int(sumDeme(deme,n_demesh) - .5*n_demesh)+1;
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
	        //sects_hist.push_back(countSectors(deme));
	        rough_hist_16.push_back(getEdgeDeme(deme,16));
	        rough_hist_32.push_back(getEdgeDeme(deme,32));
	        rough_hist_64.push_back(getEdgeDeme(deme,64));
	        rough_hist_128.push_back(getEdgeDeme(deme,128));
	        rough_hist_256.push_back(getEdgeDeme(deme,256));
	        rough_hist_512.push_back(getEdgeDeme(deme,512));
	        rough_hist_1024.push_back(getEdgeDeme(deme,1024));



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
   for(int i=0;i < n_data;i++){
   		//fvarhet << int(i*record_time) << ", "  << varhet_hist[i] << endl;
    	fhet << int(i*record_time) << ", "  << het_hist[i] << endl;
    	fpop << int(i*record_time) << ", "  << pop_hist[i] << endl;
    	//fsects<< int(i*record_time) << ", "  << sects_hist[i] << endl;
    	frough_16<< int(i*record_time) << ", "  << rough_hist_16[i] << endl;
    	frough_32<< int(i*record_time) << ", "  << rough_hist_32[i] << endl;
    	frough_64<< int(i*record_time) << ", "  << rough_hist_64[i] << endl;
    	frough_128<< int(i*record_time) << ", "  << rough_hist_128[i] << endl;
    	frough_256<< int(i*record_time) << ", "  << rough_hist_256[i] << endl;
    	frough_512<< int(i*record_time) << ", "  << rough_hist_512[i] << endl;
    	frough_1024<< int(i*record_time) << ", "  << rough_hist_1024[i] << endl;



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
    frough_16.close();
    frough_32.close();
    frough_64.close();
    frough_128.close();
    frough_256.close();
    frough_512.close();
    frough_1024.close();



    //fvarhet.close();

    cout << "Finished!" << "\n";
    cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);







	return 0;
	




}