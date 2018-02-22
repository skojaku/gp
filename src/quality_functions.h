/* 
*
 * Quaity function for a community 
 *
 */


/* Internal degree */
double calc_qint(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const int cid){
	
    int N = A.size();
    int Nu = 0;
    double retval = 0;
    for(int i = 0;i < N;i++) {
	if(xlist[cid][i]==false) continue;

	for(int j = 0;j < A[i].size();j++) {
		if(xlist[cid][i] & xlist[cid][A[i][j]]) retval+=W[i][j];
        }
	Nu+=!!(xlist[cid][i]);
    }
   if(Nu==1 | Nu ==N) return 0; 
   return retval / (double) Nu;
}

double calc_qint_diff(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const vector<double>& Nk,
    const vector<vector<double> >& Wrs,
    const vector<double>& Dk,
    const vector<vector<double> >& toU,
    const vector<double>& SelfLoop,
    const double M,
    const int nid,
    const int cid,
    const int newcid
    ){
    
	
    double toC_old = toU[nid][cid];
    double toC_new = toU[nid][newcid];
    double selfw = SelfLoop[nid];
	
   double dQ = (Wrs[cid][cid] - 2.0 * toC_old - selfw) / (Nk[cid] - 1.0) - Wrs[cid][cid] / Nk[cid] 
		+ (Wrs[newcid][newcid] + 2.0 * toC_new + selfw) / (Nk[newcid] + 1.0) - Wrs[newcid][newcid] / Nk[newcid];
    
   if(Nk[cid]-1<1e-30 | Nk[newcid]+1 >= A.size()-1e-30) return -numeric_limits<double>::max(); 
   return dQ;
}

double calc_qint_all(const vector<vector<int>>&A, const vector<vector<double>>&W, const vector<vector<bool>>&xlist){
	double retval = 0;
	for(int cid = 0; cid < xlist.size(); cid++){
		retval+=calc_qint(A, W, xlist, cid);	
	}
	return retval;
} 



/* Expansion */
double calc_qext(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const int cid 
){
	
    int N = A.size();
    int Nu = 0;
    double retval = 0;
    for(int i = 0;i < N;i++) {
	if(xlist[cid][i]==false) continue;
        
	for(int j = 0;j < A[i].size();j++) {
		if(xlist[cid][i] & xlist[cid][A[i][j]]==false) retval+=W[i][j];
        }
	Nu+=1;
    }
   return  - retval / (double) Nu;
}

double calc_qext_diff(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const vector<double>& Nk,
    const vector<vector<double> >& Wrs,
    const vector<double>& Dk,
    const vector<vector<double> >& toU,
    const vector<double>& SelfLoop,
    const double M,
    const int nid,
    const int cid,
    const int newcid
    ){
    
	
    double toC_old = toU[nid][cid];
    double toC_new = toU[nid][newcid];
    double selfw = SelfLoop[nid];
    double deg = accumulate(W[nid].begin(), W[nid].end(), 0.0);
    double Okcid = accumulate(Wrs[cid].begin(), Wrs[cid].end(), 0.0);
    double Oknewcid = accumulate(Wrs[newcid].begin(), Wrs[newcid].end(), 0.0);
   
    double dQ = (Okcid + toC_old - ( deg - toC_old - selfw)) / (Nk[cid] - 1.0) - Okcid / Nk[cid] 
		+ (Oknewcid - toC_new + (deg - toC_new - selfw )) / (Nk[newcid] + 1.0) - Oknewcid / Nk[newcid];
   
    dQ = dQ * (-1); 
   if(Nk[cid]-1<1e-30 | Nk[newcid]+1 >= A.size()-1e-30) return -numeric_limits<double>::max(); 
    return  dQ;
}

double calc_qext_all(const vector<vector<int>>&A, const vector<vector<double>>&W, const vector<vector<bool>>&xlist){
	double retval = 0;
	for(int cid = 0; cid < xlist.size(); cid++){
		retval+=calc_qext(A, W, xlist, cid);	
	}
	return retval;
} 

/* Conductance */
double calc_qcnd(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const int cid 
){

    int N = A.size();	
    double Du = 0;
    double retval = 0;
    for(int i = 0;i < N;i++) {
	if(xlist[cid][i]==false) continue;
       
	double deg = 0; 
	for(int j = 0;j < A[i].size();j++) {
		if(xlist[cid][i] & xlist[cid][A[i][j]]==false) retval+=W[i][j];
		deg+=W[i][j];
        }
	Du+= deg;
    }
    return  - retval / (double) Du;
}


double calc_qcnd_diff(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const vector<double>& Nk,
    const vector<vector<double> >& Wrs,
    const vector<double>& Dk,
    const vector<vector<double> >& toU,
    const vector<double>& SelfLoop,
    const double M,
    const int nid,
    const int cid,
    const int newcid
    ){
    
	
    double toC_old = toU[nid][cid];
    double toC_new = toU[nid][newcid];
    double selfw = SelfLoop[nid];
    double deg = accumulate(W[nid].begin(), W[nid].end(), 0.0);

    double Okcid = accumulate(Wrs[cid].begin(), Wrs[cid].end(), 0.0);
    double Oknewcid = accumulate(Wrs[newcid].begin(), Wrs[newcid].end(), 0.0);

   
    double dQ = (Okcid + toC_old - ( deg - toC_old -selfw)) / (Dk[cid] - deg) - Okcid / Dk[cid] 
		+ (Oknewcid - toC_new + (deg - toC_new - selfw)) / (Dk[newcid] + deg) - Oknewcid / Dk[newcid];


    dQ = dQ * (-1); 
    if(Dk[cid]-deg<1e-30 | Dk[newcid]+deg >= 2*M-1e-30) return -numeric_limits<double>::max(); 
    return dQ;
}

double calc_qcnd_all(const vector<vector<int>>&A, const vector<vector<double>>&W, const vector<vector<bool>>&xlist){
	double retval = 0;
	for(int cid = 0; cid < xlist.size(); cid++){
		retval+=calc_qcnd(A, W, xlist, cid);	
	}
	return retval;
} 

/* Modularity */
double calc_qmod(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const int cid){
    
    	
    double retval = 0;
    int N = A.size();	
    double D = 0;
    double M = 0;
    for(int i =0;i<N;i++){
    	double deg = accumulate(W[i].begin(), W[i].end(), 0.0);
	M+= deg;
	if(xlist[cid][i]==false) continue;

	D+=!!(xlist[cid][i]) * deg;
        for(int j =0;j<A[i].size();j++){
		if(xlist[cid][i] & xlist[cid][A[i][j]]) retval+=W[i][j];
	}
    }
    M = M / 2;
    
    retval = ( retval - D * D / (double)(2 * M) )/ (double)(2 * M);
    return retval;
}

double calc_qmod_diff(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<vector<bool>>& xlist,
    const vector<double>& Nk,
    const vector<vector<double> >& Wrs,
    const vector<double>& Dk,
    const vector<vector<double> >& toU,
    const vector<double>& SelfLoop,
    const double M,
    const int nid,
    const int cid,
    const int newcid
    ){
    
    double toC_old = toU[nid][cid];
    double toC_new = toU[nid][newcid];
    double selfw = SelfLoop[nid];
    double deg = accumulate(W[nid].begin(), W[nid].end(), 0.0);
   
   double dQ = (Wrs[cid][cid] - 2* toC_old - pow(Dk[cid] - deg, 2)/(double)(2*M) ) - (Wrs[cid][cid] - pow(Dk[cid], 2)/(double)(2*M) ); 

   dQ+= (Wrs[newcid][newcid] + 2* toC_new - pow(Dk[newcid] + deg, 2)/(double)(2*M) ) - (Wrs[newcid][newcid] - pow(Dk[newcid], 2)/(double)(2*M) ); 
   dQ/=(double)(2.0*M);
   return dQ;
}

double calc_qmod_all(const vector<vector<int>>&A, const vector<vector<double>>&W, const vector<vector<bool>>&xlist){
	double retval = 0;
	for(int cid = 0; cid < xlist.size(); cid++){
		retval+=calc_qmod(A, W, xlist, cid);	
	}
	return retval;
} 

/* degree-corrected SBM */
double klent(double x, double y){
	double retval = 0;
	
	if(x > 1e-20){
	    retval+=x*log(x);
	}
	if(y > 1e-20){
	    retval-=x*log(y);
	}
	return retval;
}


double _a(double x){
	if(x<1e-30) return 0.0;
	return 2 * x * log(x);
}
double _b(double x){
	if(x<1e-30) return 0.0;
	return x * log(x);
}

double calc_qdcsbm_diff(
    const vector<vector<int> >& A,
    const vector<vector<double>>& W,
    const vector<vector<bool>>& xlist,
    const vector<double>& Nk,
    const vector<vector<double>>& Wrs,
    const vector<double>& Dk,
    const vector<vector<double> >& toU,
    const vector<double>& SelfLoop,
    const double M,
    const int nid,
    const int r, // from
    const int s // to
    ){
   
   int K = Wrs.size(); 
   double retval = 0;
   double ki = accumulate(W[nid].begin(), W[nid].end(), 0.0);
   for(int t = 0; t < K; t++){
	if(t==s | t==r) continue;
	retval += _a(Wrs[r][t] - toU[nid][t]) - _a(Wrs[r][t]) + _a(Wrs[s][t] + toU[nid][t]) - _a(Wrs[s][t]);	
   }
   retval+= _a(Wrs[r][s] + toU[nid][r] - toU[nid][s]) - _a(Wrs[r][s]) + _b(Wrs[r][r] - 2 *(toU[nid][r] + SelfLoop[nid])) - _b(Wrs[r][r])
	 + _b(Wrs[s][s] + 2*(toU[nid][s] + SelfLoop[nid])) - _b(Wrs[s][s]) -_a(Dk[r] - ki) + _a(Dk[r]) - _a(Dk[s] + ki) + _a(Dk[s]);
   return retval;
}

double calc_qdcsbm_all(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<vector<bool>>& xlist){

    int K = xlist.size();
    int N = A.size();
	
    vector<int> C(N, 0.0);
    for(int i = 0; i < N; i++){
    	for(int k = 0; k < K; k++){
		if(xlist[k][i]){
			C[i] = k;
			break;
		}
	}
    }
    vector<vector<double>> Wrs(K, vector<double>(K, 0.0));
    vector<double> Dk(K, 0.0);
    double M = 0.0;
    for(int i = 0; i < N; i++){
    	for(int j = 0; j < A[i].size(); j++){
		Wrs[C[i]][C[j]] +=W[i][j];
		Dk[ C[i] ]+=W[i][j];
		M+=W[i][j];
    	}
    }
    M = M / 2;
    
    double retval = 0;
    for(int r = 0; r < K; r++){
    for(int s = r; s < K; s++){
	if(r==s){
		retval+= klent(Wrs[r][s]/(2.0*M), Dk[r] * Dk[s] /(4.0*M*M) );
	}else{
		retval+= 2*klent(Wrs[r][s]/(2.0*M), Dk[r] * Dk[s] /(4.0*M*M) );
	}
    }
    }
    return retval;	
}

double calc_qdcsbm(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<vector<bool>>& xlist, int cid){

    int K = xlist.size();
    int N = A.size();
	
    vector<int> C(N, 0.0);
    for(int i = 0; i < N; i++){
    	for(int k = 0; k < K; k++){
		if(xlist[k][i]){
			C[i] = k;
			break;
		}
	}
    }
    vector<vector<double>> Wrs(K, vector<double>(K, 0.0));
    vector<double> Dk(K, 0.0);
    double M = 0.0;
    for(int i = 0; i < N; i++){
    	for(int j = 0; j < A[i].size(); j++){
		Wrs[C[i]][C[j]] +=W[i][j];
		Dk[ C[i] ]+=W[i][j];
		M+=W[i][j];
    	}
    }
    M = M / 2;
    
    double retval = 0;
    for(int r = 0; r < K; r++){
    for(int s = r; s < K; s++){
	if(r==cid | s==cid){
	if(r==s){
		retval+= klent(Wrs[r][s]/(2.0*M), Dk[r] * Dk[s] /(4.0*M*M) );
	}else{
		retval+= 2*klent(Wrs[r][s]/(2.0*M), Dk[r] * Dk[s] /(4.0*M*M) );
	}
	}
    }
    }
    return retval;
}


void init_com_param(const vector<vector<int> >& A, const vector<vector<double> >& W, const vector<vector<bool>>& xlist, vector<int>& C, vector<double>& Nk, 
    vector<vector<double> >& Wrs, vector<double>& Dk, vector<vector<double>>& toU, vector<double>& SelfLoop, double& M){
	
	int N = A.size();
	int K = xlist.size();

	fill(Wrs.begin(), Wrs.end(), vector<double>(K, 0.0));
	fill(Nk.begin(), Nk.end(),0.0);
	fill(Dk.begin(), Dk.end(),0.0);
	fill(SelfLoop.begin(), SelfLoop.end(),0.0);
	fill(C.begin(), C.end(),0);
	M = 0;
	for(int i = 0; i < N; i++){
		for(int k = 0; k < K; k++){
			if(xlist[k][i]) C[i] = k;
		}
		Nk[C[i]]++;	
	}
	for(int i = 0; i < N; i++){
		double selfw = 0; 
		for(int j = 0; j < A[i].size() ; j++){
			int nei = A[i][j];
			int r = C[i];
			int s = C[nei];
			
			Wrs[r][s]+= W[i][j];	
			
			if(i==nei){
				selfw+=W[i][j];
				continue;
			}
		}
		SelfLoop[i] = selfw;
	}
	
	for(int k = 0; k < K; k++){
		for(int l = 0; l < K; l++){
			Dk[k]+= Wrs[k][l];
		}
		M+=Dk[k];
	}
	M = M  / 2.0;	
	
	for(int i = 0; i < N; i++){
		fill(toU[i].begin(), toU[i].end(),0.0);
		int deg = A[i].size();	
		for(int j = 0; j < deg; j++){
			toU[i][ C[ A[i][j] ] ]+=W[i][j];
		}
	}
}

void update_com_param(const vector<vector<int> >& A, const vector<vector<double>>& W, vector<vector<bool>>& xlist, vector<int>& C, vector<double>& Nk, vector<vector<double>>& Wrs, vector<double>& Dk, vector<vector<double>>& toU, const double& M, const int nid, const int newcid){

	int oldcid = C[nid];
	
	// update params
	double deg = accumulate(W[nid].begin(), W[nid].end(), 0.0);
	Dk[oldcid]-=deg;
	Dk[newcid]+=deg;
	Nk[oldcid]--;
	Nk[newcid]++;
	
	double win_old = 0;
	double win_new = 0;
	double selfw = 0;
	for(int j = 0; j < A[nid].size() ;j++) {
	    int nei = A[nid][j];
	    if(nid==nei){
		selfw+=W[nid][j];
		continue;
	    }
	    int s = C[nei];
	    Wrs[oldcid][s]-=W[nid][j];
	    Wrs[s][oldcid]-=W[nid][j];
	    Wrs[newcid][s]+=W[nid][j];
	    Wrs[s][newcid]+=W[nid][j];
		
	    toU[nei][ oldcid ]-=W[nid][j]; 
	    toU[nei][ newcid ]+=W[nid][j]; 
	}
	// update to xlist
	xlist[oldcid][nid] = false;
	xlist[newcid][nid] = true;
	C[nid] = newcid;	
}

/* Compute the sum of qualities */
double calc_sum_q(
	const vector<vector<int>>& A,
	const vector<vector<double>>& W,
	const vector<vector<bool>>& xlist, 
        double (*calc_q)(const vector<vector<int>>&, const vector<vector<double>>&, const vector<vector<bool>>&, int cid) 
){
	int K = xlist.size();
	double Q = 0;
	for(int k = 0; k < K; k++){
		Q+=calc_q(A, W, xlist, k);	
	}
	return Q;	
}

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

