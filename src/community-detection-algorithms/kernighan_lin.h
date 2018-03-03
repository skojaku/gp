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

void init_by_Lazy_Kernighan_Lin_Algorithm(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W,
	int K, 
	vector<int>& C, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
	int num_of_scans,
        mt19937_64& mtrnd
	);

void Kernighan_Lin_Algorithm(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W,
	vector<vector<bool>>& xlist, // need to be initialised in advance 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
        mt19937_64& mtrnd
	);

/* ----------------------- */
/* Kernighan-Lin algorithm */
/* ----------------------- */
void run_kl(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
        mt19937_64& mtrnd
	){

	// Refine by the lazy Lernighan-lin algorithm 
	int K = xlist.size();
	int N = xlist[0].size();
	vector<int> C(N,0);
	init_by_Lazy_Kernighan_Lin_Algorithm(A, W, K, C, xlist, calc_q, calc_q_del, 10, mtrnd);

	/*	
	for(int i = 0; i < N; i++){
		C[i] = i % K;
	}
	random_shuffle( C.begin(), C.end());
	for(int i = 0; i < N; i++) xlist[C[i]][i] = true;
	*/

	// Run the Kernighan-Lin algorithm
	Kernighan_Lin_Algorithm(A, W, xlist, calc_q, calc_q_del, mtrnd);
	
}

void Kernighan_Lin_Algorithm(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W,
	vector<vector<bool>>& xlist, // need to be initialised in advance 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
        mt19937_64& mtrnd
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

	// best quality
	double Q = calc_q(A, W, xlist);
	double Qmax = Q;
	
	/* Main routine */
	vector<vector<bool>>xbest = xlist; 
	vector<bool>fixed(N, false); 
	vector<vector<double>>qgain(N, vector<double>(K,0)); // gain yielded by a relabelling 
	vector<int>candidate(N,0); // candidate label 
	 
	for( int it = 0;it < N; it++){
		
		/* Initialize for relabelling procedure */
		fill(fixed.begin(),fixed.end(),false);
		bool updated = false;
		int nextNode = -1;
		int nextCom = -1;
		double nextQ = -1 * numeric_limits<double>::max();
		Q = calc_q(A, W, xlist);
		init_com_param(A, W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);
		for(int i =0; i<N; i++){
			for(int k =0; k<K; k++){
				if(k==C[i]){
					//qgain[i][k] = 0;
					qgain[i][k] = -1 * numeric_limits<double>::max();
				}else{
					qgain[i][k] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], k);
					if(nextQ < qgain[i][k]) {
						candidate[i] = k;
						nextQ = qgain[i][k];
						nextNode = i;
						nextCom = k;
					}
				}
			}
		}
		// move a node in turn		
		for( int it2 = 0; it2 < N; it2++){
		//for( int it2 = 0; it2 < N; it2++){
			// Find a move of a node that yields the maximum increment	
			double dQmax = nextQ;
			int nid = nextNode;
			int newcid = nextCom;

			if(nextQ <-1 * numeric_limits<double>::max()*0.99 ) break;	
			/*	
			for(int i =0; i<N; i++){
				if( fixed[i] ) continue;
				if( Nk[C[i]]<=2 ) continue;
		
				// propose a new label
				int k = candidate[i];	
				double qtest = qgain[i][k];
				if(qtest > dQmax){
					dQmax = qtest;	
					newcid = k;
					nid = i;
				}
			}
			*/
			fixed[ nid ] = true; // Fix the label of node nid
			if(newcid == C[nid]){
				nextQ = -1 * numeric_limits<double>::max();
				for(int i =0; i<N; i++){
					if(fixed[i]) continue;
					if(nextQ < qgain[i][candidate[i]]) {
						nextQ = qgain[i][candidate[i]];
						nextNode = i;
						nextCom = candidate[i];
					}
				}
				continue;
			}

			if((newcid < 0) | (nid <0) ) break;
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
			
		
		
			/* update gain qgain */
			nextQ = -1 * numeric_limits<double>::max();
			for(int i =0; i<N; i++){
				if(fixed[i]) continue;
	
				if(C[i]!=cid){
					qgain[i][cid] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], cid);
				}else{
					//qgain[i][cid] = 0;
					qgain[i][cid] = -1 * numeric_limits<double>::max();
				}
				
				
				if(C[i]!=newcid){
					qgain[i][newcid] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], newcid);
				}else{
			//		qgain[i][newcid] = 0;
					qgain[i][newcid] = -1 * numeric_limits<double>::max();
				}
				
				if( (qgain[i][candidate[i]] < qgain[i][cid]) & (qgain[i][cid] <= qgain[i][newcid]) ){
					candidate[i] = newcid;
				}else if( (qgain[i][candidate[i]] < qgain[i][newcid]) & (qgain[i][newcid] <= qgain[i][cid]) ){
					candidate[i] = cid;
				}

				if(nextQ < qgain[i][candidate[i]]) {
					nextQ = qgain[i][candidate[i]];
					nextNode = i;
					nextCom = candidate[i];
				}
			}

		}
		//cout<<" -------------"<<endl;
		if ( updated == false ){
			break;
		}
		
		xlist = xbest;
	}
	xlist = xbest;
}

/* ----------------------- */
/* Lazy_Kernighan_Lin_Algorithm */
/* ----------------------- */
void init_by_Lazy_Kernighan_Lin_Algorithm(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W,
	int K, 
	vector<int>& C, 
	vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&x),
        double (*calc_q_del)(const vector<vector<int> >& A, const vector<vector<double> >&W, const vector<vector<bool>>& xlist, const vector<double>& Nk, const vector<vector<double>>& Wrs, const vector<double>& Dk, const vector<vector<double>>& toU, const vector<double>& SelfLoop, const double M, const int nid, const int cid, const int newcid),
	int num_of_scans,
        mt19937_64& mtrnd
	){

	/*
	* Refine by the label switching
	*/
	
	/* Initialise variables */	
	int N = A.size();
	vector<vector<bool>> xlisttmp(K, vector<bool>(N, false)); 
	xlist.clear(); xlist = xlisttmp;
	vector<int> ndord(N);
	vector<double>Nk(K,0.0);
	vector<vector<double>>Wrs(K,vector<double>(K, 0.0));
	vector<double>Dk(K,0.0);
	vector<double>SelfLoop(N,0.0);
	vector<vector<double>> toU(N, vector<double>(K, 0.0));
	double M = 0;
	
	// assign empty nodes
	for(int i = 0; i < N; i++){
		ndord[i] = i;
		C[i] = i % K;
	}
		
	random_shuffle( C.begin(), C.end());
	
	for(int i = 0; i < N; i++) xlist[C[i]][i] = true;

	init_com_param(A, W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);

	// best quality
	double Q = calc_q(A, W, xlist);
	double dQ = 0;
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
	
	/*
	* Refine by the Kernighan-Lin (lazy) 
	*/
	
	/* Initialise variables */	
        uniform_real_distribution<double> udist(0.0, 1.0);	
	
	// best quality
	Q = calc_q(A, W, xlist);
	double Qmax = Q;
	
	/* Main routine */
	vector<vector<bool>>xbest = xlist; 
	vector<bool>fixed(N, false); 
	vector<vector<double>>qgain(N, vector<double>(K,0)); // gain yielded by a relabelling 
	vector<int>candidate(N,0); // candidate label 
	for( int it = 0;it < N; it++){
		
		/* Initialize for relabelling procedure */
		fill(fixed.begin(),fixed.end(),false);
		bool updated = false;
		int nextNode = -1;
		int nextCom = -1;
		double nextQ = -1 * numeric_limits<double>::max();
		
		Q = calc_q(A, W, xlist);
		init_com_param(A, W, xlist,  C, Nk, Wrs, Dk, toU, SelfLoop, M);
		for(int i =0; i<N; i++){
			for(int k =0; k<K; k++){
				if(k==C[i]){
					qgain[i][k] = -1 * numeric_limits<double>::max();
				}else{
					qgain[i][k] = calc_q_del(A, W, xlist, Nk, Wrs, Dk, toU, SelfLoop, M, i, C[i], k);
					if(nextQ < qgain[i][k]) {
						candidate[i] = k;
						nextQ = qgain[i][k];
						nextNode = i;
						nextCom = k;
					}
				}
			}
		}
		// move a node in turn		
		for( int it2 = 0; it2 < num_of_scans; it2++){
		//for( int it2 = 0; it2 < N; it2++){
			// Find a move of a node that yields the maximum increment	
			double dQmax = nextQ;
			int nid = nextNode;
			int newcid = nextCom;
	
		
			if((newcid < 0) | (nid <0) ) break;
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
			nextQ = -1 * numeric_limits<double>::max();
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
				
				if( (qgain[i][candidate[i]] < qgain[i][cid]) & (qgain[i][cid] <= qgain[i][newcid]) ){
					candidate[i] = newcid;
				}else if( (qgain[i][candidate[i]] < qgain[i][newcid]) & (qgain[i][newcid] <= qgain[i][cid]) ){
					candidate[i] = cid;
				}

				if(nextQ < qgain[i][candidate[i]]) {
					nextQ = qgain[i][candidate[i]];
					nextNode = i;
					nextCom = candidate[i];
				}
			}

		}
		//cout<<" -------------"<<endl;
		if ( updated == false ){
			break;
		}
		
		xlist = xbest;
	}
	xlist = xbest;
}
