/*
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>     // cout, end
#include <fstream>
*/
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

/* ----------------------- */
/* Label Switching algorithm */
/* ----------------------- */
void label_switching(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
        mt19937_64& mtrnd
	){

	/* Initialise variables */	
	int N = A.size();
	int K = xlist.size();
	vector<int> ndord(N);
	vector<double>Nk(K,0.0);
	vector<vector<double>>Wrs(K,vector<double>(K, 0.0));
	vector<double>Dk(K,0.0);
	vector<double>SelfLoop(N,0.0);
	vector<vector<double>> toU(N, vector<double>(K, 0.0));
	vector<int>C(N,0); 
	double M = 0;
	
	// assign empty nodes
	for(int i = 0; i < N; i++){
		ndord[i] = i;
		C[i] = i % K;
	}
		
	random_shuffle( C.begin(), C.end());
	
	for(int i = 0; i < N; i++) xlist[C[i]][i] = true;

	init_com_param(A, W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);

	/* Label switching */
	bool updated = false;	
	int itNum = 0;
	do{
		updated = false;
		shuffle(ndord.begin(), ndord.end(), mtrnd);
		for(int idx = 0;idx <N;idx++){
			int nid = ndord[idx];
			int ncid = C[nid];
			int newcid = ncid;
			double dQ = 0;
			
			for(int k =0; k<K; k++){
				if(k==C[nid]) continue;
				double qgain = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, nid, C[nid], k);
				if(dQ < qgain){
					newcid = k;
					dQ = qgain;
				}
			}
			
			if(dQ < 0) continue;
			if( ncid== newcid  ) continue;
			if( Nk[ncid]<=2  ) continue;
			updated = true;
			
			update_com_param( A,  W, xlist, C, Nk, Wrs, Dk, toU, M, nid, newcid);
		}
		itNum++;
	}while( (updated == true) & (itNum<=100) );
	
	for(int k = 0; k < K; k++) fill(xlist[k].begin(), xlist[k].end(), false);
	for(int i = 0; i < N; i++) xlist[C[i]][i] = true;
}
