#include "qstest.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
using namespace std;

/* Private function */
double normcdf(double value)
{
    return 0.5 + 0.5 * erf(value * M_SQRT1_2);
}

/*
 * Chung, F. & Lu, L. 
 * Connected Components in Random Graphs with Given Expected Degree Sequences. 
 * Ann. Comb. 6, 125–145 (2002). 
 *
 * Miller, J. C. & Hagberg, A. 
 * Efficient Generation of Networks with Given Expected Degrees. 
 * in Algorithms and Models for the Web Graph (eds. Frieze, A., Horn, P. & Prałat, P.) 
 * 6732 LNCS, 115–126 (Springer Berlin Heidelberg, 2011). 
 *
 * */
void generate_randomised_net(
    const vector<double>& deg,
    const vector<int>& nodes, // id of node in descending order of degree
    vector<vector<int>>& A,
    vector<vector<double>>& W,
    bool isunweighted,
    mt19937_64 mtrnd // random number generator
    )
{
    bool noSelfloop = true;
    uniform_real_distribution<double> udist(0.0, 1.0);
    int N = deg.size();
    double M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
	
	
    for(int u = 0; u < N-1; u++){
	int v = u;

	if(noSelfloop){
		v = v + 1;
	}	
	
	double p = MIN(1, deg[ nodes[u] ] * deg[ nodes[v] ] / (2.0*M) );
	while(v < N && p > 0){
		if(p!=1){
			double r = MIN( MAX(udist(mtrnd), 1e-30), 1-1e-30);
			v = v + MAX(floor( log(r) / log(1-p) ), 0.0);
		}
		if(v < N ){
			double q = MIN(deg[ nodes[u] ] * deg[ nodes[v] ] / (2.0*M), 1);
			double w = 1;
			bool addEdge = false;
			if(isunweighted){
				double r = udist(mtrnd);
				addEdge = r < q / p;	
			}else{
	    			poisson_distribution<int> distribution(q / p);
	    			w = distribution(mtrnd);
				addEdge = w>0; 	
			}
			if(addEdge){
				if(u!=v){
	            			A[nodes[u]].push_back(nodes[v]);
	            			A[nodes[v]].push_back(nodes[u]);
	            			W[nodes[u]].push_back(w);
	            			W[nodes[v]].push_back(w);
				}else{
	            			A[nodes[u]].push_back(nodes[u]);
	            			W[nodes[u]].push_back(2*w);
				}
			}
			p = q;
			v = v + 1;
		}
	}
    }
}



/* ---- Estimating statistical significance of a community ----*/
void estimate_statistical_significance(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<vector<bool>>& xlist,
    double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
    double (*calc_qind)(const vector<vector<int>>&, const vector<vector<double>>&, const  vector<vector<bool>>&, int k),
    double (*calc_s)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<bool>&),
    void (*find_communities)(const vector<vector<int> >&, const vector<vector<double>>&, vector<vector<bool>>&, mt19937_64& ),
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values,
    vector<int>& nhat,
    vector<double>& qhat, 
    vector<int>& rgindex // index of randomised networks. 
	){

    /* Initialise variables */
    int K = xlist.size();
    int N = A.size();
    nhat.clear();
    qhat.clear();
    rgindex.clear();
    vector<double> deg(N);	
    bool isunweighted = true;
    for (int i = 0; i < N; i++){
	deg[i] = accumulate(W[i].begin(), W[i].end(), 0.0);
	if((pow(deg[i]-W[i].size(),2))>1){
		isunweighted = false;
	}
    }
    
    vector<double> q(K);
    vector<int> n(K);
    for (int k = 0; k < K; k++) {
	n[k] = calc_s(A, W, xlist[k]);
	q[k] = calc_qind(A, W, xlist, k);
    };
    
    vector<int> nodes(N);
    iota(nodes.begin(), nodes.end(), 0);
    sort(
        nodes.begin(),
        nodes.end(),
        [&](int x, int y){return deg[x] > deg[y];}
    );
   
    // create random number generator per each thread
    int numthread = 1;
    #ifdef _OPENMP
    	# pragma omp parallel
    	{
    		numthread = omp_get_num_threads();
    	}
    #endif
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mt19937_64 mtrnd = init_random_number_generator();
	mtrnd_list[i] = mtrnd;
    }
    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    
    #ifdef _OPENMP
    #pragma omp parallel for shared(nhat, qhat, deg, nodes, mtrnd_list)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {

        // Generate a randomised network using the configuration model.
        vector<vector<int>> A_rand(N, vector<int>(0));
        vector<vector<double>> W_rand(N, vector<double>(0));
        
        int tid = 0;
    	#ifdef _OPENMP
        	tid = omp_get_thread_num();
    	#endif
	
	
        mt19937_64 mtrnd = mtrnd_list[tid];
       	generate_randomised_net(deg, nodes, A_rand, W_rand, isunweighted, mtrnd);
        // Detect core-periphery pairs using the KM algorithm
	bool valid; 
	vector<vector<bool>> x_rand;
        do{
	
		// repeat com detect num_of_runs times 
		double Qr = -numeric_limits<double>::max();
		x_rand.clear();
		vector<int> cbest;
		for (int r = 0; r < num_of_runs; r++){
	            vector<vector<bool>> x_rand_tmp(K, vector<bool>(N, false) );
			
		    /* added by kojaku
         	    for(int i = 0; i < N; i++){
      			int newcid = std::uniform_int_distribution<>(0, K-1)(mtrnd);
      			x_rand_tmp[newcid][i] = true;
      		    }
                    */
	            
		    find_communities(A_rand, W_rand, x_rand_tmp, mtrnd);

		    double Qi = calc_q(A_rand, W_rand, x_rand_tmp);
		    if(Qi != Qi) {continue;}
		    
		    if (Qi > Qr) {
		        Qr = Qi;
		        x_rand = x_rand_tmp;
		    }
		}

	        
		// validate the detected community
		valid = true;
	        double tmp = calc_q(A_rand, W_rand, x_rand );
		if(tmp != tmp) {valid = false;}
	}while(valid==false);

        // Save the quality and size of core-periphery pairs in the randomised network.
        int K = x_rand.size();
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        for(int k = 0; k < K; k++) {
            nhat.push_back( calc_s(A_rand, W_rand, x_rand[k]) );
            qhat.push_back( calc_qind(A_rand, W_rand, x_rand, k) );
            rgindex.push_back( it );
        }
    }
    /* Compute mean and covariance */
    int S = nhat.size();
    double mu_n = (double)accumulate(nhat.begin(), nhat.end(), 0.0) / (double)S;
    double mu_q = (double)accumulate(qhat.begin(), qhat.end(), 0.0) / (double)S;
    double sig_nn = 0;
    double sig_qq = 0;
    double sig_nq = 0;
    for (int s = 0; s < S; s++) {
        sig_nn += pow((double)nhat[s] - mu_n, 2) / (double)(S - 1);
        sig_qq += pow(qhat[s] - mu_q, 2) / (double)(S - 1);
        sig_nq += ((double)nhat[s] - mu_n) * (qhat[s] - mu_q) / (double)(S - 1);
    }

    /* Compute p-values using the Gaussian kernel density estimator */
    double h = MAX(pow((double)S, -1.0 / 6.0), 1e-32);
    p_values.clear();
    p_values.assign(K, 1.0);
    for (int k = 0; k < K; k++) {
        double numer = 0.0;
        double denom = 0.0;
        for (int s = 0; s < S; s++) {
            double qbar = qhat[s] + sig_nq / sig_nn * (double)(n[k] - nhat[s]);
	
            double t = sig_nn * (q[k] - qbar) / (sqrt(sig_nn * sig_qq - sig_nq * sig_nq ) * h);
            double cum = normcdf(t);
	    
            double w = exp(-(double)pow(n[k] - nhat[s], 2) / (2.0 * h * h * sig_nn )) + 1e-33;
            numer += cum * w;
            denom += w;
        }
        p_values[k] = 1.0 - numer / denom;
    }
}

/* ---- Size function for a community ----*/
double calc_n(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<bool>& x){

    double cnt = 0;
    unsigned int N = x.size();
    for(unsigned int i = 0;i<N;i++) cnt+=!!(x[i]);
    return cnt;
    //return count(x.begin(),x.end(),true); 
}

double calc_e(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<bool>& x){

    double cnt = 0;
    int N = x.size();
    for(int i = 0;i< N ;i++){
	if(x[i]) cnt+=accumulate(W[i].begin(), W[i].end(), 0.0);
    }
    return cnt;
}

typedef double (*size_func)(const vector<vector<int> >&, const vector<vector<double>>&, const vector<bool>&);
map<string, size_func> size_functions = {
        {"nodes", calc_n},
        {"edges", calc_e}
};

