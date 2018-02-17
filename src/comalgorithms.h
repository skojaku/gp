/* 
*
* Optimisation algorithms 
*
* */
#include "mcmc.h"
#include "kernighan_lin.h"

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

void mymcmc(const vector<vector<int> >&A, const vector<vector<double> >&W, vector<vector<bool>>& xlist){
	run_mcmc(A, W, xlist, mcmc_qfunc, mcmc_qfunc_diff);
}

// wrapper for the Kernighan--Lin algorithm 
void mykl(const vector<vector<int> >&A, const vector<vector<double> >&W, vector<vector<bool>>& xlist){
	run_kl(A, W, xlist, mcmc_qfunc, mcmc_qfunc_diff);	
}

typedef void (*com_detect_func)(const vector<vector<int> >&, const vector<vector<double>>&, vector<vector<bool>>&);
map<string, com_detect_func> community_detection = {
	{"kl", mykl},
	{"mcmc", mymcmc}
};

/* ---- Quality functions ----*/
typedef double (*quality_func)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&);
map<string, quality_func> quality_functions = {
	{"int", calc_qint_all},
	{"rcut", calc_qext_all},
	{"ncut", calc_qcnd_all},
	{"mod", calc_qmod_all},
	{"dcsbm", calc_qdcsbm_all}
};

typedef double (*quality_func_ind)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&, int );
map<string, quality_func_ind> quality_functions_ind = {
	{"int", calc_qint},
	{"rcut", calc_qext},
	{"ncut", calc_qcnd},
	{"mod", calc_qmod},
	{"dcsbm", calc_qdcsbm}
};

typedef double (*quality_func_diff)(
		const vector<vector<int>>& A, 
		const vector<vector<double>>& W, 
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

map<string, quality_func_diff> quality_functions_diff = {
	{"int", calc_qint_diff},
	{"rcut", calc_qext_diff},
	{"ncut", calc_qcnd_diff},
	{"mod", calc_qmod_diff},
	{"dcsbm", calc_qdcsbm_diff}
};

