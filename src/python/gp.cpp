#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<algorithm>
#include <iostream>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "../lib/gp.h"

using namespace std;
namespace py = pybind11;

void readEdgeTable(py::array_t<double> edges_array_t, vector<vector<int>>& A, vector<vector<double>>& W, int& N, int& M)
{

    vector<int> edgeList;
    vector<double> wList;
    N = 0;	
    auto edges = edges_array_t.data();
    auto r = edges_array_t.request();
    M = r.shape[0];

    for(int i =0; i< M; i++){
        int sid = int(edges[3*i]) - 1;
        int did = int(edges[3*i+1]) - 1;
	double w = edges[3*i+2];
	
        if (sid == did)
            continue;

        if (N < sid)
            N = sid;
        if (N < did)
            N = did;
        edgeList.push_back(sid);
        edgeList.push_back(did);
        wList.push_back(w);
    }
    N = N + 1;
   
    vector<vector<int>>tmp(N);
    A = tmp; 
    vector<vector<double>>tmp2(N);
    W = tmp2; 

    int wid = 0; 
    int edgeListsize = edgeList.size();
    for (int i = 0; i < edgeListsize; i += 2) {
        int sid = edgeList[i];
        int did = edgeList[i + 1];
	double w = wList[wid];
	wid++;
        A[sid].push_back(did);
        A[did].push_back(sid);
        W[sid].push_back(w);
        W[did].push_back(w);
    }
}

auto detect(py::array_t<double> edges, int K, string qfunc, string algorithm, int num_of_runs, double significance_level, string sfunc, int num_of_rand_nets, py::array_t<int> Cinit_array_t ){
        int N = 0;	
        int M = 0;	
        vector<vector<int> > A;
        vector<vector<double>> W;
	readEdgeTable(edges, A, W, N, M);
    
    	mt19937_64 mtrnd;
    	random_device r;
    	seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    	mtrnd.seed(seed);
        auto Cinit = Cinit_array_t.data();
	
    	vector<vector<bool>> xlist;
    	double Qr = -numeric_limits<double>::max();
    	mcmc_qfunc = quality_functions[qfunc];
    	mcmc_qfunc_diff = quality_functions_diff[qfunc];
    	for (int r = 0; r < num_of_runs; r++) {
     		vector<vector<bool>> xlist_tmp(K, vector<bool>(N, false) );
	
		for(int i = 0; i < N; i++){
			xlist_tmp[Cinit[i]][i] = true;
		}
		
		
    		community_detection[algorithm](A, W, xlist_tmp, mtrnd);
        	double Qi = quality_functions[qfunc](A, W, xlist_tmp);
        	if(Qi != Qi) {continue;}
        
        	if (Qi > Qr) {
        	    Qr = Qi;
        	    xlist = xlist_tmp;
        	}
    	}
	K = xlist.size();
	double corrected_significance_level =1.0 - pow(1.0 - significance_level, 1.0 / (double)K); // Sidak correction.
	vector<double> p_values(K, 0.0);
	vector<int> nhat;
	vector<double> qhat;
	vector<int> rgindex;
	if(significance_level < 1){
	    	fill(p_values.begin(), p_values.end(), 0.0);
	    	mcmc_qfunc = quality_functions[qfunc];
	    	mcmc_qfunc_ind = quality_functions_ind[qfunc];
	    	mcmc_qfunc_diff = quality_functions_diff[qfunc];
		estimate_statistical_significance(A, W, xlist, mcmc_qfunc, mcmc_qfunc_ind, size_functions[sfunc], community_detection[algorithm], num_of_runs, num_of_rand_nets, p_values, nhat, qhat, rgindex);
	}

	py::array_t<double> cids_array_t(N);
	auto cids = cids_array_t.mutable_data();
	
	py::array_t<double> pvals_array_t(N);
	auto pvals = pvals_array_t.mutable_data();
	
	py::array_t<bool> sigs_array_t(N);
	auto sigs = sigs_array_t.mutable_data();	

	K = xlist.size();	
	for(int i = 0; i < N; i++){
		
		for(int k = 0; k < K; k++){
			if(xlist[k][i]){
				cids[i] = k;
				break;
			}	
		}
		pvals[i] = p_values[cids[i]];
		sigs[i] = p_values[cids[i]]<=corrected_significance_level;
	}

	py::list results(3);
	results[0] = cids_array_t;
	results[1] = pvals_array_t;
	results[2] = sigs_array_t;
	return results;
}

PYBIND11_MODULE(gp, m){
	m.doc() = "Community detection in networks";
	m.def("detect", &detect, "Detect communities in networks",
	py::arg("edges"),
	py::arg("K") = 2, 
	py::arg("qfunc") = "mod", 
	py::arg("algorithm") = "louvain", 
	py::arg("num_of_runs") = 10,
	py::arg("significance_level") = 1.0, 
	py::arg("sfunc") = "edges",
	py::arg("num_of_rand_nets") = 500,
	py::arg("Cinit")
	);
}
