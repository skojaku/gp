#include "mex.h"
#include <map>
#include <cmath>
#include <cfloat>

#include "gp.h"
#include "quality_functions.h"
#include "community-detection-algorithms/comalgorithms.h"

void init_random_number_generator(){
    	int seeds[624];
       	size_t size = 624*4; //Declare size of data
       	ifstream urandom("/dev/urandom", ios::in | ios::binary); //Open stream
       	if (urandom) //Check if stream is open
       	{
       	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
       	    urandom.close(); //close stream
       	}
       	else //Open failed
       	{
            		cerr << "Failed to open /dev/urandom" << endl;
       	}
    	seed_seq seed(&seeds[0], &seeds[624]);
    	mtrnd.seed(seed);
}



/* ---- Mex functions ----*/
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    double* edgeList = mxGetPr(prhs[0]);
    int N = (int)mxGetPr(prhs[1])[0];
    int Enum = (int)mxGetPr(prhs[2])[0];
    int K =  (int)mxGetPr(prhs[3])[0];
    int num_of_runs =  (int)mxGetPr(prhs[4])[0];
    string qfunc_name = mxArrayToString(prhs[5]); // quality function 
    string alg_name = mxArrayToString(prhs[6]); // community detection algorithm 

    /* Parse Input */
    vector<vector<int> > A(N, vector<int>(0));
    vector<vector<double> > W(N);
    for (int i = 0; i < Enum; i++) {
        int rid = (edgeList[i] - 1);
        int cid = (edgeList[i + (int)Enum] - 1);
        double w = (edgeList[i + 2*(int)Enum]);
        if (rid != cid) {
            A[rid].push_back(cid);
            A[cid].push_back(rid);
            W[rid].push_back(w);
            W[cid].push_back(w);
        }
    }
    /* Detect K communities in networks */
    init_random_number_generator();
    mcmc_qfunc = quality_functions[qfunc_name];
    mcmc_qfunc_diff = quality_functions_diff[qfunc_name];
    vector<vector<bool>> xlist;

    // repeat com detect num_of_runs times 
    double Qr = -numeric_limits<double>::max();
    for (int r = 0; r < num_of_runs; r++) {
        vector<vector<bool>> xlist_tmp(K, vector<bool>(N, false) );
    	community_detection[alg_name](A, W, xlist_tmp);
    
        double Qi = mcmc_qfunc(A, W, xlist_tmp);
        if(Qi != Qi) {continue;}
        
        if (Qi > Qr) {
            Qr = Qi;
            xlist = xlist_tmp;
        }
    }

    K = xlist.size();	
  
    /* Save results*/
    plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);
    double* cids = mxGetPr(plhs[0]);
    double* qs = mxGetPr(plhs[1]);

    for (int i = 0; i < N; i++) {
	int cid = -1;
	for(int k = 0; k < K; k++){
		if(xlist[k][i]){
			cid =k;
			break;
		}
	}
        cids[i] = cid + 1;
    }
    for(int k = 0; k < K; k++){
	double qk = quality_functions_ind[qfunc_name](A, W, xlist, k);
	qs[k] = qk;
    }
}
