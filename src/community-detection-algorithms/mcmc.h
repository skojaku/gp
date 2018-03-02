/*
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>     // cout, end
*/
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

using namespace std;

void run_mcmc(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
	mt19937_64& mtrnd
	){
	
	// initialise variables
        uniform_real_distribution<double> udist(0.0, 1.0);	
	int N = A.size();
	int K = xlist.size();
	vector<double>Nk(K,0.0);
	vector<vector<double>>Wrs(K,vector<double>(K, 0.0));
	vector<double>Dk(K,0.0);
	vector<double>SelfLoop(N,0.0);
	vector<vector<double>> toU(N, vector<double>(K, 0.0));
	vector<int>ord(N,0);
	vector<int>C(N,0); 
	double M = 0;

	// mcmc params	
	double beta = 1; 
	int maxStableLoopNum = 100;
	int maxItNum = 1000;
	int itnum=0;int lastupdate = 0;
	
	// randomise 
	for(int i = 0; i < N; i++){
		int cid = -1;
		ord[i] = i;
		for(int k = 0; k < K; k++){
			if(!xlist[k][i]) continue;
			cid = k;
			break;
		}
		if(cid<0){
			int newcid = std::uniform_int_distribution<>(0, K-1)(mtrnd);
			xlist[newcid][i] = true;
		}
	}
	init_com_param(A,  W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);
	
        vector<vector<bool>> xbest=xlist; 
	
	// best quality
	double Q = calc_q(A, W, xlist);
	double Qbest = Q;
	double Qmin = Q;
	
	// mcmc
	while( (itnum <=maxItNum) & ((itnum - lastupdate)<=maxStableLoopNum) ){
	
		shuffle(ord.begin(), ord.end(), mtrnd);
		
		for(int it = 0;it < N; it++){
			int nid = ord[it];	
		
			// propose a new label
			//int nid = std::uniform_int_distribution<>(0, N-1)(mtrnd);
			int newcid = std::uniform_int_distribution<>(0, K-1)(mtrnd);
	
			int cid = C[nid]; 
			if(cid == newcid){
				continue;
			}
				
			// compute the increment 	
			double dQ = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, nid, C[nid], newcid);
			
			Qmin = MIN(Q + dQ, Qmin);
		
			if(  udist(mtrnd) <= exp(beta * dQ/ (Qbest - Qmin) ) ){// accept	
				Q = Q + dQ;
				
				update_com_param(A, W, xlist, C, Nk, Wrs, Dk, toU, M, nid, newcid);
				
				if(Q > Qbest){
					xbest = xlist;	
					lastupdate = itnum;			
					Qbest = Q;
				}	
			}
		}
		itnum +=1;	
		beta = beta * 1.5;
	}
	xlist = xbest;
}
