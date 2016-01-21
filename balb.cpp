/*
 * Copyright 1/16/2016 Quan Zhou, Baylor College of Medicine 
 * Joint work with Philip Ernst and L. C. G. Rogers 
 *
 * This code is published under GPL-3 licenses. Please refer to License.txt for details
*/

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sstream>
#include <algorithm>

using namespace std;

const double RERR = 1e-9; 
// Introduce noise on the first day so that bins differ. Small enough compared to h
// But not useful for change of measure since return will often be 1. (unless everytime you introduce such an error)

const double ZERO = 1e-12;
const double MAXV = 999999.00;
typedef pair<int, double> pid;
const gsl_rng_type* gsl_type;
gsl_rng* gsl_r;

// double r = 0; // FORCE r = 0 
double h = 0.1 / 250.0;
double sigma = 1.0;  // no need to change
int term = 250;  // expiration day
int ahead = 5;   // lookback day
int bin = 200;  // number of bins
int b = 200;   // number of samples per bin
int n_low = 50000;   
int n_high = 10000;
int n_sub = 50;   // used in high bound simulations
double alpha = 0.5 * sigma * sigma; // here is the change of measure. Otherwise it is -0.5
double sd = sigma * sqrt(h);  
double p_unit; 

vector < vector < vector <double>* >* > transition; // Day 0 to Day term-1
vector < vector <double>* > bin_med; // Day 0 to Day term 
vector < vector <double>* > bin_bound; // 
vector < vector <int> > bin_one; // 1 is special. many bins have median/bound = 1  
vector < vector <double>* > rule_bound;  // merged from bin_bound
vector < vector <int>* > stop_rule;   // low bound stopping rule
vector < vector <double>* > max_return;  // value function 

ofstream LOG; 
ofstream LOW; 
ofstream HIGH; 

bool cmp_pid (pid a, pid b){return(a.second < b.second);}

// Max of a vector 
double vec_max (vector<double>* vec){
	double m = vec->at(0);
	for (int i=1; i<vec->size(); i++){
		if (vec->at(i) > m){ 
			m = vec->at(i); 
		}
	}
	return m;
}

// Discretization. Calculate bin median and bounds. Return bin labels (ordered).
void discretize( vector<double>& v, vector<int>& lab, int day){
	if (v.size() != b*bin){
		cerr << "The vector size is not equal to bin_size*n_bin!" << endl;	
	}
	
	vector < pid > vc;
	for (int i=0; i<v.size(); i++){
		pid ele = make_pair(i, v[i]);
		vc.push_back(ele);
	}
	sort(vc.begin(), vc.end(), cmp_pid);
	for (int i=0; i<v.size(); i++){
		int bn = i/b;
		lab[ vc[i].first ] = bn;
	}

	for (int i=0; i < bin; i++){
		double med =( vc[(i+0.5)*b - 1].second + vc[(i+0.5)*b].second) / 2.0;
		bin_med[day]->at(i) = med;
		if (i < bin-1){
			bin_bound[day]->at(i) = (vc[(i+1)*b - 1].second + vc[(i+1)*b].second)/2.0; 
		}else{
			bin_bound[day]->at(i) = vc[(i+1)*b-1].second;	
//			bin_bound[day]->at(i) = MAXV;	
		}
	}
	
	LOG << "Bin_median\tDay" << day << "\t"; 
	for (int i=0; i < bin; i++){
		LOG << "\t" << bin_med[day]->at(i);
	}
	LOG << "\nBin_bounds\tDay" << day << "\t";
	for (int i=0; i < bin; i++){
		LOG << "\t" << bin_bound[day]->at(i);
	}
	LOG << endl;
	
	return;
}

void calculate_transition(vector <int>& w1, vector <int>& w0, int day ){
	for (int i=0; i< w1.size(); i++){	
		transition[day-1]->at(w0[i])->at(w1[i]) += p_unit;
	}
	return;
}

void initialize (){
	for (int t=0; t<=term; t++){
		vector <double>* t1 = new vector <double>;
		vector <double>* t2 = new vector <double>;
		vector <double>* t3 = new vector <double>;
		for (int i=1; i <= bin; i++){
			t1->push_back(0.0);
			t2->push_back(0.0);
			t3->push_back(0.0);
		}
		bin_med.push_back(t1);
		bin_bound.push_back(t2);
		max_return.push_back(t3);
		vector <double>* tmp_b = new vector <double>;
		vector <int>* tmp_s = new vector <int>;
		rule_bound.push_back(tmp_b);
		stop_rule.push_back(tmp_s);
	}

	for (int t=1; t<=term; t++){
		vector < vector<double>* >* t1 = new vector< vector <double>* >;
		for (int i=1; i<=bin; i++){
			vector < double >* t2 = new vector <double>;
			for (int j=1; j<=bin; j++){
				t2->push_back(0.0);	
			}
			t1->push_back(t2);
		}
		transition.push_back(t1);
	}	
	return;
}

void transition_simulation (int ns){
	cout << "Calcuating transition matrix   ";
	
	vector<double> v (ns, 1.0);
	vector<double> z (ns, 1.0);
	vector<int> label (ns, 0);
	vector< vector <double> > lookback;

/********************** Day 0 ******************/
	for (int i=0; i<ns; i++){	
//		double e0 = exp ( gsl_rng_uniform(gsl_r) *RERR );  // Give up
		v[i] = 1.0;  z[i] = 1.0;
		vector<double> tmp (ahead + 1, 0.0);
		tmp[ahead] = 1.0;   // always no noise here
		lookback.push_back(tmp);
	}
	discretize(z, label, 0);
	vector<int> label_last = label;

/**************Day 1 to Day Expiration ******************/
	int print_unit = term / 20;
	for (int t=1; t <= term; t++){
		if (t % print_unit == 0){
			cout << "=";
			fflush(stdout);
		}
		for (int i=0; i < ns; i++){		
			double x = gsl_ran_gaussian(gsl_r, sd) + alpha*h;		
			double s = lookback[i][ahead] * exp(x);
			if (v[i] == lookback[i][0]){
				lookback[i].erase(lookback[i].begin());
				lookback[i].push_back(s);
				v[i] = vec_max(&lookback[i]);
			}else{
				lookback[i].erase(lookback[i].begin());
				lookback[i].push_back(s);
				if (s > v[i]){v[i] = s;}
			}
			z[i] = v[i] / s;
		}
		discretize(z, label, t);
		calculate_transition(label, label_last, t);
		label_last = label;
	}
	cout << endl;
}

void low_bound (void){
	cout << "Calculating low-bound stopping rule" << endl ;
	vector < vector <int> > stop; // Temporary stopping rule
	for (int t=0; t<=term; t++){
		vector <int> t2;
		for (int i=1; i <= bin; i++){ t2.push_back(0); }
		stop.push_back(t2);
	}

	LOG << "Day" << term << " ";
	for (int i=0; i<bin; i++){
		max_return[term]->at(i) = bin_med[term]->at(i);	
		stop[term][i] = 1;
		LOG << "\t" << max_return[term]->at(i) ;
	}
	LOG << endl;

	for (int t = term-1; t >= 0; t--){	
		LOG << "Day" << t << " ";
		for (int i=0; i<bin; i++){			
			double sr = 0.0;
			for (int j=0; j<bin; j++){
				sr += ( transition[t]->at(i)->at(j) * max_return[t+1]->at(j) );
			}
			if (sr > bin_med[t]->at(i)){
				max_return[t]->at(i) = sr;
			}else{
				max_return[t]->at(i) = bin_med[t]->at(i);
				stop[t][i] = 1;
			}
			LOG << "\t" << sr;
		}
		LOG << endl;
	}
	
	double sim_max = 0.0;
	for (int i=0; i<bin; i++){
		sim_max += (max_return[0]->at(i) / (double) bin);	
	}	
	LOG << "Max expected return over the sampled paths = " << sim_max << endl;

	int merge = 0;
	int warns = 0;
	for (int t = 0; t <= term; t++){
		double bound = bin_bound[t]->at(0);
		int if_s = stop[t][0];
		int count_a = 1;
		int count_b = 0;
		int mark = 0;
		LOG << "Rule at Day " << t << "\t";
		for (int i=1; i<bin; i++){
			if (abs(bin_bound[t]->at(i) - bound)<ZERO){
				mark = 1;
				if (if_s == stop[t][i]){count_a ++;}
				else {count_b ++;}
				continue;
			}else{
				if (mark == 1){
					if (count_b == 0){
						merge ++;
					//	cout << "Successfully merged bins on day " << t << endl;
					}else{
						warns ++;
					//	cerr << "Warning! at Day " << t << " and bound = " << bound << ": " << count_a << " vs " << count_b << "; Stop = " << if_s << endl;
					}
					if (count_a < count_b){
						if_s = 1 - if_s;	
					}
					count_a = 1;
					count_b = 0;
					mark = 0;
				}	
			}

			if (if_s == stop[t][i]){
				bound = bin_bound[t]->at(i);
			}else{
				rule_bound[t]->push_back(bound);
				stop_rule[t]->push_back(if_s);
				LOG << "  " << bound << ":" << if_s << ";";
				bound = bin_bound[t]->at(i);
				if_s = 1 - if_s;
			}
		}	
		rule_bound[t]->push_back(MAXV);
		stop_rule[t]->push_back(if_s);
		LOG << "  " << bound << ":" << if_s << ";" << endl;
	}	


	if (merge > 0){
		cout << "Successfully merged bins on " << merge << " days" << endl;	
	}
	if (warns > 0){
		cout << warns << " warnings!" << endl;	
	}
	
	return;	
}

int If_Stop(double v, int day){
	int max = rule_bound[day]->size() - 1;
	int min = 0;
	int now = (max+min)/2;
	while (1){
		if (v <= rule_bound[day]->at(now)){
			if (now == 0){break;}
			if (v > rule_bound[day]->at(now-1)){
				break;	
			}else{
				max = now - 1;	
			}
		}else{
			min = now+1;
		}
		if (min == max){
			now = min;
			break;
		}else{
			now = (max+min)/2;
		}
	}
	return stop_rule[day]->at(now);
}

void return_simulation (int ns){
	cout << "Running low bound simulations  ";
	double ss = 0.0;

	int print_unit = ns / 20;
	for (int i=1; i<=ns; i++){	
		if (i % print_unit == 0){
			cout << '=' ;
			fflush(stdout);
		}
		double v = 1.0;
		double z = 1.0;
		vector <double> lookback (ahead + 1, 0.0);
		lookback[ahead] = 1.0;  
		if ( If_Stop(z, 0) == 1){
			LOW << z << "\t" << 0 << endl;
			ss += z;
			continue;
		}
		
		for (int t=1; t <= term; t++){
			double x = gsl_ran_gaussian(gsl_r, sd) + alpha*h;		
			double s = lookback[ahead] * exp(x);		
			if (v == lookback[0]){
				lookback.erase(lookback.begin());
				lookback.push_back(s);
				v = vec_max(&lookback);
			}else{
				lookback.erase(lookback.begin());
				lookback.push_back(s);
				if (s > v){v = s;}
			}
			z = v/s;
			if ( If_Stop(z, t) == 1){
				LOW << z << "\t" << t << endl;
				ss += z;
				break;
			}
		}
	}
	double mean = ss / (double) ns;
	cout << "\nEstimated lower bound = " << mean << endl;
	LOG << "Estimated lower bound = " << mean << endl;
	return;
}

void find_bin_one(void){
	for (int t=0; t<=term; t++){
		vector <int> lab (2, 0);
		int mark = 0;
		for (int i=0; i<bin; i++){
			if ( abs(bin_bound[t]->at(i) - 1.0)<ZERO ){
				if (mark == 0){
					mark = 1;
					lab[0] = i;
				}
			}else{
				if (mark == 1){
					lab[1] = i - 1; 		
					break;
				}
			}
		}
		bin_one.push_back(lab);
		LOG << "Day " << t << "; one-start: " << lab[0] << "; one-end: " << lab[1] << endl;
	}
}

int search_bin(double z, int day){
	if (abs(z-1.0)<ZERO){
		int k = bin_one[day][0] + gsl_rng_uniform_int(gsl_r, bin_one[day][1] - bin_one[day][0] + 1);
	//	LOG << "z = " << z << "; day = " << day << "; start = " << bin_one[day][0] << "; end = " << bin_one[day][1] << "; k = " << k<< endl;
		return k;
	}

	int max = bin_bound[day]->size() - 1;
	int min = 0;
	int now = (max+min)/2;
	while (1){
		if ( abs(z - bin_bound[day]->at(now)) < ZERO ){
		// this happens when z == 1 (mostly). Return a random bin by uniform sampling. 
			int start = min;
			int end = max;
			for (int i = (now + 1); i <= max; i ++ ){
				if ( abs(z - bin_bound[day]->at(i)) > ZERO ){
					end = i - 1;
					break;
				}
			}
			for (int i = now-1; i >= min; i -- ){
				if ( abs(z - bin_bound[day]->at(i)) > ZERO ){
					start = i + 1;
					break;
				}
			}
			if (start == end){
				return start;	
			}else{
				int k = start + gsl_rng_uniform_int(gsl_r, end - start + 1);
				LOG << "z = " << z << "; day = " << day << "; start = " << start << "; end = " << end << "; k = " << k<< endl;
				return k;
			}
		}else if (z < bin_bound[day]->at(now)){
			if (now == 0){
				break;
			}
			if (z > bin_bound[day]->at(now-1)){
				break;	
			}else{
				max = now - 1;	
			}
		}else{
			min = now+1;
		}
		if (min == max){
			now = min;
			break;
		}else{
			now = (max+min)/2;
		}
	}
	return now;
}


void high_bound (int ns){
	cout << "Running high bound simulations  ";
		
	int print_unit = ns/20;
	double ss = 0.0;
	for (int i=1; i<=ns; i++){
		if (i % print_unit == 0){
			cout << '=';
			fflush(stdout);
		}
//		double e0 = exp( gsl_rng_uniform(gsl_r)*RERR );	
		double v = 1.0;
		vector <double> lookback(ahead+1, 0.0);
		vector <double> z (term+1, 1.0);
		vector <double> m (term+1, 0.0);
		lookback[ahead] = 1.0;
		double snell = 1.0;
		int sday = 0;

//		LOG << " " << m[0];
		for (int t=1; t <= term; t++){
			double x = gsl_ran_gaussian(gsl_r, sd) + alpha*h;		
			double s0 = lookback[ahead];
			double v0 = v; // min reward at next day
			double s = s0 * exp(x);		
			if (v == lookback[0]){
				lookback.erase(lookback.begin());
				v0 = vec_max(&lookback); // current 2nd max
				lookback.push_back(s);
				v = vec_max(&lookback);
			}else{
				lookback.erase(lookback.begin());
				lookback.push_back(s);
				if (s > v){v = s;}
			}
			z[t] = v/s;
			int z_b = search_bin(z[t], t);
			double mr_now = max_return[t]->at(z_b);

			double emr = 0.0;
			for (int j=0; j<n_sub; j++){	
				double x0 = gsl_ran_gaussian(gsl_r, sd) + alpha*h;		
				double s1 = s0 * exp(x0);
				double z1 = v0 / s1;
				if (s1 >= v0){ 
					z1 = 1.0;
				}
				int zb1 = search_bin(z1, t);
				emr += max_return[t]->at(zb1);
			}
			emr = emr /( (double) n_sub);
			m[t] = m[t-1] + mr_now - emr;
//			LOG << ", " << m[t];
			if (z[t] - m[t] > snell){
				snell = z[t] - m[t];	
				sday = t;
			}
		}
//		LOG << endl;
		HIGH << snell << "\t" << sday << endl;
		ss += snell;
	}	
	double mean = ss /( (double) n_high );
	cout << "\nEstimated upper bound = " << mean << endl;
	LOG << "Estimated upper bound = " << mean << endl;
}

void Help (void){
	cout << "Example command:\n    ./blb -a 5 -o test\n";
	cout << "This will run simulation for a=5(days) and generate output files: test.high.out and test.low.out. They contain the reward and stopping day of each upper/lower bound simulation, from which you can calculate the mean and standard error.\n";
	cout << "\nList of arguments (in the square brackets is the argument type and default value)" << endl;
	cout << "-o [string='bermuda']: the output file prefix" << endl;
	cout << "-a [int=5]: the number of look-back days" << endl;
	cout << "-t [int=250]: the expiration day" << endl;
	cout << "-h [double=4e-4]: time unit per day" << endl;
	cout << "-u [int=10000]: the number of simulations for calculating upper bound" << endl;
	cout << "-l [int=50000]: the number of simulations for calculating lower bound" << endl;
	cout << "-b [int=200]: the number of bins" << endl;
	cout << "-d [int=200]: the number of samples per bin" << endl;
	cout << "-s [int=50]: the number of sub-simulations in upper bound simulation" << endl;	
	exit(EXIT_FAILURE);	
}


int main (int argc, char** argv){
	gsl_rng_env_setup();
	gsl_type = gsl_rng_default;
	gsl_r = gsl_rng_alloc(gsl_type);

	if (argc > 1){
		if (argv[1][0] == 'h'){
			Help();
		}
	}

	string out = "bermuda";
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			char arg = argv[i][1];
			switch (arg){
				case 'o':
					out = argv[i+1]; break;
				case 'a':
					ahead = atoi(argv[i+1]); break;
				case 'b':
					bin = atoi(argv[i+1]); break;
				case 'd':
					b = atoi(argv[i+1]); break;
				case 'h':
					h = atof(argv[i+1]); break;	
				case 'l':
					n_low = atoi(argv[i+1]); break;
				case 's':
					n_sub = atoi(argv[i+1]); break;
				case 't':
					term = atoi(argv[i+1]); break;
				case 'u':
					n_high = atoi(argv[i+1]); break;
				default:
					cout << "Wrong argument!" << endl;
					exit(EXIT_FAILURE);
			}
			i ++ ;
		}	
	}
	

	string log_file = out;
	log_file.append(".log");	
	string low_file = out;
	low_file.append(".low.out");	
	string high_file = out;
	high_file.append(".high.out");
	
	p_unit = 1.0/((double) b);

	LOG.open(log_file.c_str());	
	LOW.open(low_file.c_str());	
	HIGH.open(high_file.c_str());	
		
	initialize();
	transition_simulation(b * bin);
	low_bound();
	return_simulation(n_low);
	find_bin_one();
	high_bound(n_high);

	LOG.close();
	LOW.close();
	HIGH.close();

	return 0;
}



