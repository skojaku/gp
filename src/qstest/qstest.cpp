#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>
#ifdef _OEPNMP
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

/* Private function */
void generate_randomised_net(
    vector<double>& deg,
    vector<vector<int>>& A,
    vector<vector<double>>& W)
{
    int N = deg.size();
    double M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
    A.clear();
    W.clear();
    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
        vector<double> tmp2;
        W.push_back(tmp2);
    }

    double rnd = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            poisson_distribution<int> distribution(deg[i] * deg[j] / (2.0 * (double)M));
            rnd = distribution(mtrnd);
            if (rnd <1 ) continue;
            A[i].push_back(j);
            W[i].push_back((double)rnd);
            
            if(i == j) continue;
                
            A[j].push_back(i);
            W[j].push_back((double)rnd);
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
    void (*find_communities)(const vector<vector<int> >&, const vector<vector<double>>&, vector<vector<bool>>&),
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values,
    vector<int>& nhat,
    vector<double>& qhat){

    /* Initialise variables */
    int K = xlist.size();
    int N = A.size();
    double Q;
    nhat.clear();
    qhat.clear();
    vector<double> deg(N);
    for (int i = 0; i < N; i++) deg[i] = accumulate(W[i].begin(), W[i].end(), 0.0);
    
    vector<double> q(K);
    vector<int> n(K);
    for (int k = 0; k < K; k++) {
	n[k] = calc_s(A, W, xlist[k]);
	q[k] = calc_qind(A, W, xlist, k);
    };

    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    #ifdef _OPENMP
    #pragma omp parallel for shared(nhat, qhat, deg)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {

        // Generate a randomised network using the configuration model.
        vector<vector<int>> A_rand;
        vector<vector<double>> W_rand;
        generate_randomised_net(deg, A_rand, W_rand);

        // Detect core-periphery pairs using the KM algorithm
	bool valid; 
	vector<vector<bool>> x_rand;
        int K_rand = x_rand.size();
        do{
	
		// repeat com detect num_of_runs times 
		double Qr = -numeric_limits<double>::max();
		x_rand.clear();
		vector<int> cbest;
		for (int r = 0; r < num_of_runs; r++) {
	            vector<vector<bool>> x_rand_tmp(K, vector<bool>(N, false) );
		    
	            find_communities(A_rand, W_rand, x_rand_tmp);

		    double Qi = calc_q(A_rand, W_rand, x_rand_tmp);
		    if(Qi != Qi) {continue;}
		    
		    if (Qi > Qr) {
		        Qr = Qi;
		        x_rand = x_rand_tmp;
		    }
		}

	        
		// validate the detected community
		valid = true;
        	int K_rand = x_rand.size();
	        double tmp = calc_q(A_rand, W_rand, x_rand );
		if(tmp != tmp) {valid = false;}
	}while(valid==false);

        // Save the quality and size of core-periphery pairs in the randomised network.

        #ifdef _OPENMP
        #pragma omp critical
        #endif
        for (int k = 0; k < x_rand.size(); k++) {
            nhat.push_back( calc_s(A_rand, W_rand, x_rand[k]) );
            qhat.push_back( calc_qind(A_rand, W_rand, x_rand, k) );
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
    for(int i = 0;i<x.size();i++) cnt+=!!(x[i]);
    return cnt;
    //return count(x.begin(),x.end(),true); 
}

double calc_e(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<bool>& x){

    double cnt = 0;
    for(int i = 0;i<x.size();i++){
	if(x[i]) cnt+=accumulate(W[i].begin(), W[i].end(), 0.0);
    }
    return cnt;
}

typedef double (*size_func)(const vector<vector<int> >&, const vector<vector<double>>&, const vector<bool>&);
map<string, size_func> size_functions = {
        {"nodes", calc_n},
        {"edges", calc_e}
};

