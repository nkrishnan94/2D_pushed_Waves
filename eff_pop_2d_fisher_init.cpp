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


unsigned long K  = 10000;
unsigned int n_gens = 5000;
//float h_thresh = .1;
const int n_demesh = 120; 
const int n_demesw =900; 
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

	for(int j=0 ;j< n_demesw; j++){

		int frontFound=0;
		for(int i = 1; i < n_demesh; i++){
			if( (arr[i][j][0]+arr[i][j][0]==0) && (arr[i-1][j][0]+arr[i-1][j][0]>0) &&(frontFound==0)){
				fronts[j]=i-1;
				frontSum+=i-1;
				frontFound=1;

			}


		}
	}

	for(int j =int(n_demesw/2 -L/2) ; j< int(n_demesw/2+L/2); j++){
		frontDiffs+=  pow(((frontSum/N) - fronts[j]),2);

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
	vector <double> rough_hist_10;
	vector <double> rough_hist_50;
	vector <double> rough_hist_100;
	vector <double> rough_hist_150;
	vector <double> rough_hist_200;
	vector <double> rough_hist_250;
	vector <double> rough_hist_300;
	vector <double> rough_hist_350;
	vector <double> rough_hist_400;
	vector <double> rough_hist_450;
	vector <double> rough_hist_500;
	vector <double> rough_hist_550;
	vector <double> rough_hist_600;
	vector <double> rough_hist_650;
	vector <double> rough_hist_700;
	vector <double> rough_hist_750;
	vector <double> rough_hist_800;
	vector <double> rough_hist_850;
	//vector <double> varhet_hist;

	//data files
	ofstream flog, fpop, fhet, fprof,fsects,frough_10,frough_50,frough_100,frough_150,frough_200,frough_250,frough_300,frough_350,frough_400,frough_450,frough_500,frough_550,frough_600,frough_650,frough_700,frough_750,frough_800,frough_850;
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
	string rough10Name = "rough_10_" + param_string + date_time.str() + ".txt";
	string rough50Name = "rough_50_" + param_string + date_time.str() + ".txt";
	string rough100Name = "rough_100_" + param_string + date_time.str() + ".txt";
	string rough150Name = "rough_150_" + param_string + date_time.str() + ".txt";
	string rough200Name = "rough_200_" + param_string + date_time.str() + ".txt";
	string rough250Name = "rough_250_" + param_string + date_time.str() + ".txt";
	string rough300Name = "rough_300_" + param_string + date_time.str() + ".txt";
	string rough350Name = "rough_350_" + param_string + date_time.str() + ".txt";
	string rough400Name = "rough_400_" + param_string + date_time.str() + ".txt";
	string rough450Name = "rough_450_" + param_string + date_time.str() + ".txt";
	string rough500Name = "rough_500_" + param_string + date_time.str() + ".txt";
	string rough550Name = "rough_550_" + param_string + date_time.str() + ".txt";
	string rough600Name = "rough_600_" + param_string + date_time.str() + ".txt";
	string rough650Name = "rough_650_" + param_string + date_time.str() + ".txt";
	string rough700Name = "rough_700_" + param_string + date_time.str() + ".txt";
	string rough750Name = "rough_750_" + param_string + date_time.str() + ".txt";
	string rough800Name = "rough_800_" + param_string + date_time.str() + ".txt";
	string rough850Name = "rough_850_" + param_string + date_time.str() + ".txt";
	string folder = "KPZ/sim_data_KPZ/";
	//string folder = "";

    flog.open(folder+logName);
    fhet.open(folder+hetName);
    fpop.open(folder+popName);
    fprof.open(folder+profName);
    //fsects.open(folder+sectsName);
    frough_10.open(folder+rough10Name);
    frough_50.open(folder+rough50Name);
    frough_100.open(folder+rough100Name);
    frough_150.open(folder+rough150Name);
    frough_200.open(folder+rough200Name);
    frough_250.open(folder+rough250Name);
    frough_300.open(folder+rough300Name);
    frough_350.open(folder+rough350Name);
    frough_400.open(folder+rough400Name);
    frough_450.open(folder+rough450Name);
    frough_500.open(folder+rough500Name);
    frough_550.open(folder+rough550Name);
    frough_600.open(folder+rough600Name);
    frough_650.open(folder+rough650Name);
    frough_700.open(folder+rough700Name);
    frough_750.open(folder+rough750Name);
    frough_800.open(folder+rough800Name);
    frough_850.open(folder+rough850Name);
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
	      deme[0][j][0] = stof(line);
	      //cout << stoi(line);
	      deme[0][j][1] = 0;
	      j+=1;
	      //cout << line[1,2,3];

	   
	    }
	    myfile.close();


	}
	


	for(int i = 1; i < int(n_demesh); i++){
		for(int j = 0; j < n_demesw; j++){

			deme[i][j][0] = deme[0][j][0];
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
	        //sects_hist.push_back(countSectors(deme));
	        rough_hist_10.push_back(getFrontDeme(deme,10));
	        rough_hist_50.push_back(getFrontDeme(deme,50));
	        rough_hist_100.push_back(getFrontDeme(deme,100));
	        rough_hist_150.push_back(getFrontDeme(deme,150));
	        rough_hist_200.push_back(getFrontDeme(deme,200));
	        rough_hist_250.push_back(getFrontDeme(deme,250));
	        rough_hist_300.push_back(getFrontDeme(deme,300));
	        rough_hist_350.push_back(getFrontDeme(deme,350));
	        rough_hist_400.push_back(getFrontDeme(deme,400));
	        rough_hist_450.push_back(getFrontDeme(deme,450));
	        rough_hist_500.push_back(getFrontDeme(deme,500));
	        rough_hist_550.push_back(getFrontDeme(deme,550));
	        rough_hist_600.push_back(getFrontDeme(deme,600));
	        rough_hist_650.push_back(getFrontDeme(deme,650));
	        rough_hist_700.push_back(getFrontDeme(deme,700));
	        rough_hist_750.push_back(getFrontDeme(deme,750));
	        rough_hist_800.push_back(getFrontDeme(deme,800));
	        rough_hist_850.push_back(getFrontDeme(deme,850));

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
    	frough_10<< int(i*record_time) << ", "  << rough_hist_10[i] << endl;
    	frough_50<< int(i*record_time) << ", "  << rough_hist_50[i] << endl;
    	frough_100<< int(i*record_time) << ", "  << rough_hist_100[i] << endl;
    	frough_150<< int(i*record_time) << ", "  << rough_hist_150[i] << endl;
    	frough_200<< int(i*record_time) << ", "  << rough_hist_200[i] << endl;
    	frough_250<< int(i*record_time) << ", "  << rough_hist_250[i] << endl;
    	frough_300<< int(i*record_time) << ", "  << rough_hist_300[i] << endl;
    	frough_350<< int(i*record_time) << ", "  << rough_hist_350[i] << endl;
    	frough_400<< int(i*record_time) << ", "  << rough_hist_400[i] << endl;
    	frough_450<< int(i*record_time) << ", "  << rough_hist_450[i] << endl;
    	frough_500<< int(i*record_time) << ", "  << rough_hist_500[i] << endl;
    	frough_550<< int(i*record_time) << ", "  << rough_hist_550[i] << endl;
    	frough_600<< int(i*record_time) << ", "  << rough_hist_600[i] << endl;
    	frough_650<< int(i*record_time) << ", "  << rough_hist_650[i] << endl;
    	frough_700<< int(i*record_time) << ", "  << rough_hist_700[i] << endl;
    	frough_750<< int(i*record_time) << ", "  << rough_hist_750[i] << endl;
    	frough_800<< int(i*record_time) << ", "  << rough_hist_800[i] << endl;
    	frough_850<< int(i*record_time) << ", "  << rough_hist_850[i] << endl;

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
    frough_10.close();
    frough_50.close();
    frough_100.close();
    frough_150.close();
    frough_200.close();
    frough_250.close();
    frough_300.close();
    frough_350.close();
    frough_400.close();
    frough_450.close();
    frough_500.close();
    frough_550.close();
    frough_600.close();
    frough_650.close();
    frough_700.close();
    frough_750.close();
    frough_800.close();
    frough_850.close();
    //fvarhet.close();

    cout << "Finished!" << "\n";
    cout << "Final Heterozygosity: " <<het_hist[n_data-1] << "\n";
    cout << "Finished in " << run_time << " seconds \n";

	puts (buffer);







	return 0;
	




}