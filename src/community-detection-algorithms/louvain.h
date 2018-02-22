/*
#include <random>
#include <algorithm>
#include <vector>
#include <random>
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


void louvain(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const int num_of_runs,
    vector<vector<bool>>& xlist);

double calc_dQmod( double d, double degi, double D, double selfw, const double M ){
	return 2*( d - degi*D/(2.0*M)) + (selfw-degi*degi/(2.0*M)); 
}

double calc_Qmod( 
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<int>& C, 
	const double M
	 ){
	
	double retval = 0;
	int N = C.size();
	vector<double> degC(N);fill(degC.begin(),degC.end(),0);
	for(int i =0;i<N;i++){
	for(int j =0;j<A[i].size();j++){
		degC[ C[i] ]+=W[i][j];
		if(C[i] == C[A[i][j]]) {
			retval+=W[i][j];	
		}	
	}
	}
	
	for(int i =0;i<N;i++){
		retval-=degC[i]*degC[i]/(2*M);	
	}
	retval/=(2*M);
	return retval;
}




void coarsing(
	vector<vector<int>>& A, 
	vector<vector<double>>& W, 
	vector<int>& C
	){
	
	int K = 0; 
	for(int i = 0;i<C.size();i++){
		if( K<C[i] )K = C[i];
	}
	K = K +1;
	vector<vector<int>> newA(K);
	vector<vector<double>> newW(K);
	for(int i = 0;i<K;i++){
		vector<int> tmp;
		vector<double> tmp2;
		newA[i] = tmp;
		newW[i] = tmp2;
	} 
	
	for(int i = 0;i<A.size();i++){
	for(int j = 0;j<A[i].size();j++){
		int rid = C[i];
		int cid = C[A[i][j]];
		double w = W[i][j];
		
		int l = -1;
		for( int k = 0;k<newA[rid].size();k++){	
			if(newA[rid][k]==cid){
				l = k;	
				break;
			}
		}	
		if(l<0){
			newA[rid].push_back(cid);	
			newW[rid].push_back(w);
		}else{
			newW[rid][l]+=w;
		}
	}
	}
	int N = A.size();
	for(int i = 0;i<N;i++){
		A[i].clear();
		W[i].clear();
	}
	A.clear(); 
	W.clear(); 
	
	A = newA;	
	W = newW;	
}

void modularity_label_switching(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<int>& C, 
	const double M
	){
	
        int N=C.size();
	vector<int> ndord(N);
	vector<double> D(N);fill(D.begin(),D.end(),0);
	vector<double>deg(N);fill(deg.begin(),deg.end(),0);
	for(int i = 0;i < N;i++) {
		ndord[i] = i;
		deg[i] = accumulate(W[i].begin(), W[i].end(), 0.0);
		D[C[i]]+=deg[i];
	};
		
	// --------------------------------
	// Label switching  
	// --------------------------------
	vector<double> toC;
	toC.assign(N,0.0);
	bool isupdated = false;
	int itNum = 0;
	do{
		isupdated = false;
		shuffle(ndord.begin(), ndord.end(),mtrnd);
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
			int ncid = C[nid];
			
			int neighbourNum = A[nid].size();
			fill(toC.begin(),toC.end(),0.0);	
				
			int newcid = ncid;
			double selfw = 0;
			for(int j = 0;j<neighbourNum;j++){		
				int nei = A[nid][j];
				int cid = C[nei];
				if(nid==nei){
					selfw+=W[nid][j];
				}else{
					toC[cid]+= W[nid][j];
				}	
			}
			
			double dQold = calc_dQmod( toC[ncid], deg[nid], D[ncid] - deg[nid], selfw, M );
			double dQ = 0;
			for(int j = 0;j<neighbourNum;j++){		
				int nei = A[nid][j];
				int cid = C[nei];
				if(nei==nid) continue;
					
				double dQc=calc_dQmod( toC[cid], deg[nid], D[cid] - deg[nid]*(double)!!(ncid==cid), selfw, M )-dQold;
			
				if( dQc<dQ ) continue;
					
				newcid = cid;
				dQ = dQc;	
			
			}
			
			if(dQ< 0) continue;
				
			if( ncid==newcid  ) continue;
			
			
			D[ncid]-=deg[nid];
			D[newcid]+=deg[nid];
			C[nid] = newcid;	
			
			isupdated = true;
		}
		itNum++;
	}while( isupdated == true );
	// remove redundant cp
	std::vector<int> labs;
	for(int i=0;i<N;i++){
		int cid = -1;
		for(int j=0;j<labs.size();j++){
			if(labs[j]==C[i]){
				cid = j;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(C[i]);
			cid = labs.size()-1;
		}
		C[i] = cid;		
	}
}

void louvain_core(
	const vector<vector<int>>& A, 
	const vector<vector<double>>& W, 
	vector<int>& C, 
	const double M
	){
	
	
	
	vector<vector<int>> newA = A; 
	vector<vector<double>> newW = W;
	vector<int>Zt = C; 
	vector<int>Ct = C;
	int prevGraphSize = C.size();
	double Qbest = calc_Qmod(newA, newW, Zt, M); 
	do{
		prevGraphSize = newA.size();
		
	 	modularity_label_switching(newA, newW, Zt, M);
		double Qt = calc_Qmod(newA, newW, Zt, M);	
		coarsing(newA,newW,Zt);
		
		// update C
		// Ct = Ct*Zt;
		for(int i = 0;i<Ct.size();i++){
			Ct[i] = Zt[ Ct[i] ];
		}
		
		if(Qt>Qbest){
			C = Ct;
			Qbest = Qt;
		}
	}while( newA.size()!= prevGraphSize);
		
}

void louvain(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const int num_of_runs,
    vector<vector<bool>>& xlist){

    int N = A.size();
    vector<int> C(N);

    for(int k = 0;k < xlist.size();k++){
        xlist[k].clear();
    }
    xlist.clear();

    double M = 0;
    for(int i =0;i<N;i++){
    	C[i] = i;
	M+=accumulate(W[i].begin(), W[i].end(), 0.0);
    }
    M = M / 2; 

    double Q = -1;
    vector<int> cbest;
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci = C;
        double Qi = 0.0;

        louvain_core(A, W, ci, M);
	Qi = calc_Qmod(A, W, ci, M);
	
        if (Qi > Q) {
            Q = Qi;
            cbest = ci;
        }
    }
    int K = *max_element(cbest.begin(),cbest.end()) + 1; 
    for(int k = 0; k < K; k++){
	vector<bool> tmp(N);
    	for(int i = 0; i < N; i++){
		tmp[i] = cbest[i]==k;
    	}	
	xlist.push_back(tmp);	
    }
}
