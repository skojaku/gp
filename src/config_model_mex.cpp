#include "config_model.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
   	double* degList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
	int num_of_rand_nets = (int) mxGetPr(prhs[2])[0]; 
	int noSelfLoop = (int) mxGetPr(prhs[3])[0]; 
	int isunweighted = (int) mxGetPr(prhs[4])[0]; 
	
	vector<double> deg(N);
	for(int i = 0;i<N;i++) deg[i] = degList[i];

	init_random_number_generator();
        
	vector<vector<vector<int>>> Alist;
        vector<vector<vector<double>>> Wlist;
	generate_randomised_nets(deg, num_of_rand_nets, Alist, Wlist, noSelfLoop==1, isunweighted==1);

	// pack to a list
	vector<vector<double>> ret;
	for(int i = 0; i < Alist.size(); i++){
		for(int j = 0; j < Alist[i].size(); j++){
			for(int k = 0; k < Alist[i][j].size(); k++){
				vector<double> tmp(4);
				tmp[0] = j + 1;
				tmp[1] = Alist[i][j][k] + 1;
				tmp[2] = Wlist[i][j][k];
				tmp[3] = i + 1;
				ret.push_back(tmp);
			}
		}
	}

	int S = ret.size(); 
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)S, (mwSize)4, mxREAL); 
	
    	double* R = mxGetPr(plhs[0]);
	
	for(int i=0;i<S;i++){
		R[i] = ret[i][0];	
		R[i + S] = ret[i][1];	
		R[i + 2 * S] = ret[i][2];	
		R[i + 3 * S] = ret[i][3];	
	}
}
