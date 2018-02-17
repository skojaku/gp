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
/* Kernighan-Lin algorithm */
/* ----------------------- */
void run_kl(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid)
	){

	/* Initialise variables */	
        uniform_real_distribution<double> udist(0.0, 1.0);	
	int N = A.size();
	int K = xlist.size();
	vector<double>Nk(K,0.0);
	vector<vector<double>>Wrs(K,vector<double>(K, 0.0));
	vector<double>Dk(K,0.0);
	vector<double>SelfLoop(N,0.0);
	vector<vector<double>> toU(N, vector<double>(K, 0.0));
	vector<int>C(N,0); 
	double M = 0;
	
	// assign empty nodes
	for(int i = 0; i < N; i++){
		int cid = -1;
	
		for(int k = 0; k < K; k++){
			if(!xlist[k][i]) continue;
			cid = k;
			C[i] = k;
			break;
		}
		if(cid<0){
			cid = std::uniform_int_distribution<>(0, K-1)(mtrnd);
			C[i] = cid;
			xlist[cid][i] = true;
		}
	}

	// best quality
	double Q = calc_q(A, W, xlist);
	double Qmax = Q;
	double Qmin = Q;
	
	/* Main routine */
	vector<vector<bool>>xbest = xlist; 
	vector<bool>fixed(N, false); 
	vector<vector<double>>qgain(N, vector<double>(K,0)); // gain yielded by a relabelling 
	vector<int>candidate(N,0); // candidate label 
	 
	for( int it = 0;it < N; it++){
		
		/* Initialize for relabelling procedure */
		fill(fixed.begin(),fixed.end(),false);
		bool updated = false;	
		Q = calc_q(A, W, xlist);
		init_com_param(A, W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);
		for(int i =0; i<N; i++){
			for(int k =0; k<K; k++){
				if(k==C[i]){
					qgain[i][k] = -1 * numeric_limits<double>::max();
				}else{
					qgain[i][k] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], k);
				}
			}
			candidate[i] = distance( qgain[i].begin(), max_element(qgain[i].begin(), qgain[i].end() ) );	
		}
	
		// move a node in turn		
		for( int it2 = 0; it2 < N; it2++){
			double qmax = -1 * std::numeric_limits<double>::max();
				
			// Find a move of a node that yields the maximum increment	
			double dQmax = -1 * numeric_limits<double>::max();
			int nid = -1;
			int newcid = -1;
			for(int i =0; i<N; i++){
				if( fixed[i] ) continue;
				if(Nk[C[i]]<=2) continue;
		
				/* propose a new label */
				int k = candidate[i];	
				double qtest = qgain[i][k];
	    			if(qtest != qtest ) continue; // if qtest is not a number, e.g., inf or NaN

				if(qtest > dQmax){
					dQmax = qtest;	
					newcid = k;
					nid = i;
				}
			}
		
		
			if(newcid < 0 | nid <0) break;
			int cid = C[nid];
			update_com_param( A,  W, xlist, C, Nk, Wrs, Dk, toU, M, nid, newcid);
			
			// update maximum
			Q = Q + dQmax;
			
			// Save X if W is the maximum so far	
			if(Qmax < Q){
				xbest = xlist;
				Qmax = Q;
				updated = true;
			}
			
			fixed[ nid ] = true; // Fix the label of node nid
		
			/* update gain qgain */
			for(int i =0; i<N; i++){
				if(fixed[i]) continue;
	
				if(C[i]!=cid){
					qgain[i][cid] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], cid);
				}else{
					qgain[i][cid] = -1 * numeric_limits<double>::max();
				}
				
				
				if(C[i]!=newcid){
					qgain[i][newcid] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], newcid);
				}else{
					qgain[i][newcid] = -1 * numeric_limits<double>::max();
				}
				
				if(qgain[i][candidate[i]] < qgain[i][cid] & qgain[i][cid] <= qgain[i][newcid]){
					candidate[i] = newcid;
				}else if(qgain[i][candidate[i]] < qgain[i][newcid] & qgain[i][newcid] <= qgain[i][cid]){
					candidate[i] = cid;
				}
			}

		}
		//cout<<" -------------"<<endl;
		if ( updated == false ){
			break;
		}
		
		xlist = xbest;
	}
	//cout<<Qmax<<endl;
	xlist = xbest;
	/*	
	int nnc = count(X.begin(), X.end(), true);
	assert(nnc == nc);
	*/	
}
