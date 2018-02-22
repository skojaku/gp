/* 
*
* Optimisation algorithms 
*
* */
#include "mcmc.h"
#include "kernighan_lin.h"
#include "louvain.h"

/* ---- Algorithm for finding communities in networks ----*/
double (*mcmc_qfunc)(const vector<vector<int> >&, const vector<vector<double> >&, const vector<vector<bool>>&);
double (*mcmc_qfunc_ind)(const vector<vector<int> >&, const vector<vector<double> >&, const vector<vector<bool>>&, int);
double (*mcmc_qfunc_diff)(
	const vector<vector<int> >& A, 
	const vector<vector<double> >& W, 
	const vector<vector<bool>>& xlist, 
	const vector<double>& Nk, 
	const vector<vector<double>>& Wrs, 
	const vector<double>& Dk, 
	const vector<vector<double>>& toU, 
	const vector<double>& SelfLoop, 
	const double M, 
	const int nid, 
	const int cid, 
	const int newcid);


/* ---- Wrapper functions ----*/
void mymcmc(const vector<vector<int> >&A, const vector<vector<double> >&W, vector<vector<bool>>& xlist){
	run_mcmc(A, W, xlist, mcmc_qfunc, mcmc_qfunc_diff);
}

void mykl(const vector<vector<int> >&A, const vector<vector<double> >&W, vector<vector<bool>>& xlist){
	run_kl(A, W, xlist, mcmc_qfunc, mcmc_qfunc_diff);	
}

void mylouvain(const vector<vector<int> >&A, const vector<vector<double> >&W, vector<vector<bool>>& xlist){
        louvain(A, W, 1, xlist);
}


typedef void (*com_detect_func)(const vector<vector<int> >&, const vector<vector<double>>&, vector<vector<bool>>&);
map<string, com_detect_func> community_detection = {
	{"kl", mykl},
	{"mcmc", mymcmc},
	{"louvain", mylouvain}
};

