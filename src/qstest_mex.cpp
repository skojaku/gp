#include "mex.h"
#include <map>
#include <cmath>
#include <cfloat>
#include "gp.h"


/* ---- Mex functions ----*/
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    double* edgeList = mxGetPr(prhs[0]);
    int N = (int)mxGetPr(prhs[1])[0];
    int Enum = (int)mxGetPr(prhs[2])[0];
    double* initC =  mxGetPr(prhs[3]);
    string qfunc_name = mxArrayToString(prhs[4]); // quality function 
    string sfunc_name = mxArrayToString(prhs[5]); // size function
    string alg_name = mxArrayToString(prhs[6]); // community detection algorithm 
    int num_of_runs = (int)mxGetPr(prhs[7])[0];
    int num_of_rand_nets = (int)mxGetPr(prhs[8])[0];

    /* Parse Input */
    vector<vector<int>> A(N, vector<int>(0));
    vector<vector<double>> W(N, vector<double>(0));
    int K = 0;
    for (int i = 0; i < N; i++) {
    	K = max(K, (int)(initC[i]));
    }

    for (int i = 0; i < Enum; i++) {
        int rid = (edgeList[i] - 1);
        int cid = (edgeList[i + (int)Enum] - 1);
        double w = edgeList[i + 2*(int)Enum];
        A[rid].push_back(cid);
        W[rid].push_back(w);
        if (rid != cid) {
            A[cid].push_back(rid);
            W[cid].push_back(w);
        }
    }
    vector<double> nlist(K,0);
    vector<double> qlist(K,0); 
    vector<vector<bool>>xlist;
    for(int k = 0; k < K; k++){
	vector<bool> tmp(N, false);
    	for(int i = 0; i < N; i++){
		tmp[i] = ((int)initC[i]-1)==k;
    	}
	xlist.push_back(tmp);
	qlist[k]=quality_functions_ind[qfunc_name](A, W, xlist, k);
	nlist[k]=size_functions[sfunc_name](A, W, xlist[k]);
    }

    vector<double> p_values;
    vector<int> nhat;
    vector<int> rgindex;
    vector<double> qhat;
	
    mcmc_qfunc = quality_functions[qfunc_name];
    mcmc_qfunc_ind = quality_functions_ind[qfunc_name];
    mcmc_qfunc_diff = quality_functions_diff[qfunc_name];
    estimate_statistical_significance(A, W, xlist, mcmc_qfunc, mcmc_qfunc_ind, size_functions[sfunc_name], community_detection[alg_name], num_of_runs, num_of_rand_nets, p_values, nhat, qhat, rgindex);

    int R = nhat.size();

    plhs[0] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix((mwSize)R, (mwSize)1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix((mwSize)R, (mwSize)1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix((mwSize)R, (mwSize)1, mxREAL);

    double* retPvals = mxGetPr(plhs[0]);
    double* retn = mxGetPr(plhs[1]);
    double* retq = mxGetPr(plhs[2]);
    double* retnhat = mxGetPr(plhs[3]);
    double* retqhat = mxGetPr(plhs[4]);
    double* retrgindex = mxGetPr(plhs[5]);

    for (int i = 0; i < R; i++) {
        retnhat[i] = nhat[i];
        retqhat[i] = qhat[i];
        retrgindex[i] = rgindex[i];
    }

    for (int k = 0; k < K; k++) {
        retPvals[k] = p_values[k];
        retn[k] = nlist[k];
        retq[k] = qlist[k];
    }
}
