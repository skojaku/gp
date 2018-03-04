/*
*
* Command-line client for the Graph partitioning algorithms
*
*
* Please do not distribute without contacting the authors above.
* If you find a bug in this code, please contact the authors.
*
*
* AUTHOR - Sadamori Kojaku
*
*
* DATE - 17 Feb, 2018
*/


#include <map>
#include <cmath>
#include <cfloat>
#include "../lib/gp.h"

#include <string>
#include <stdio.h>
#include <unistd.h>

char delimiter = '\t';

void split(const string& s, char c,
    vector<string>& v);

void readEdgeTable(string filename, vector<vector<int> >& A, vector<vector<double>>& W, int& N);

void writeLabels(const string filename, const vector<vector<bool>>& xlist, const vector<double>p_values, double alpha);

void usage();


int main(int argc, char* argv[])
{
    if (argc <= 2) {
        usage();
        return 1;
    }

    string linkfile = argv[1];
    string outputfile = argv[2];
	
    int num_of_rand_nets = 500; 
    int num_of_runs = 10;
    int K = 2;
    double alpha = 0.05;

    int opt;
    opterr = 0;
    string tmp;
    string alg_name = "kl";
    string qfunc_name = "dcsbm";
    string sfunc_name = "edges";
    while ((opt = getopt(argc, argv, "h:k:r:a:q:s:o:l:")) != -1) {
        switch (opt) {
        case 'h':
            usage();
            break;
        case 's':
            tmp.assign(optarg);
            sfunc_name = tmp.c_str();
            break;
        case 'a':
            tmp.assign(optarg);
            alpha = atof(tmp.c_str());
            break;
        case 'o':
            tmp.assign(optarg);
            alg_name = tmp.c_str();
            break;
        case 'q':
            tmp.assign(optarg);
            qfunc_name = tmp.c_str();
            break;
        case 'r':
            tmp.assign(optarg);
            num_of_runs = atoi(tmp.c_str());
            break;
        case 'l':
            tmp.assign(optarg);
            num_of_rand_nets = atoi(tmp.c_str());
            break;
        case 'k':
            tmp.assign(optarg);
            K = atoi(tmp.c_str());
            break;
        default: /* '?' */
            usage();
            break;
        }
    }

    cout << "===================" << endl;
    cout << "# input file:  "<< linkfile<< endl;
    cout << "# output file: "<< outputfile<< endl;
    cout << "" << endl;
  
    cout << "# Status: "<< endl;
    cout << "   Reading input file...";
    /* Read file */
    int N;
    vector<vector<int> > A;
    vector<vector<double>> W;
    readEdgeTable(linkfile, A, W, N);
    cout <<"end"<<endl<<endl;
      

    /* Detect communities in networks */
    cout << "   Seeking communities..."<<endl;
    cout << "      - algorithm: "<< alg_name<< endl;
    cout << "      - quality function: "<< qfunc_name<< endl;
    cout << "      - Number of runs: "<< num_of_runs<< endl;
    cout << "      - Number of communities: "<< K<< endl;
    mt19937_64 mtrnd;
    random_device r;
    seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    mtrnd.seed(seed);
    vector<vector<bool>> xlist;
    double Qr = -numeric_limits<double>::max();
    mcmc_qfunc = quality_functions[qfunc_name];
    mcmc_qfunc_diff = quality_functions_diff[qfunc_name];
    for (int r = 0; r < num_of_runs; r++) {
        vector<vector<bool>> xlist_tmp(K, vector<bool>(N, false) );
    	community_detection[alg_name](A, W, xlist_tmp, mtrnd);
        double Qi = quality_functions[qfunc_name](A, W, xlist_tmp);
        if(Qi != Qi) {continue;}
        
        if (Qi > Qr) {
            Qr = Qi;
            xlist = xlist_tmp;
        }
    }
    cout <<"   end"<<endl<<endl;
	
    /* Significance test */
    K = xlist.size();
    double corrected_alpha =1.0 - pow(1.0 - alpha, 1.0 / (double)K); // Sidak correction.
    vector<double> p_values(K, 0.0);
    vector<int> nhat;
    vector<double> qhat;
    vector<int> rgindex;
    if(alpha < 1){
    	cout << "   Significance test..."<<endl;
    	cout << "      - size function: "<< sfunc_name<< endl;
    	cout << "      - number of random networks: "<< num_of_rand_nets<< endl;
    	cout << "      - significance level: "<< alpha<< endl;
    	cout << "      - corrected-significance level: "<< corrected_alpha<< endl;
    	cout << "      - number of communities under testing: "<< K<< endl;
        fill(p_values.begin(), p_values.end(), 0.0);
        mcmc_qfunc = quality_functions[qfunc_name];
        mcmc_qfunc_ind = quality_functions_ind[qfunc_name];
        mcmc_qfunc_diff = quality_functions_diff[qfunc_name];
    	estimate_statistical_significance(A, W, xlist, mcmc_qfunc, mcmc_qfunc_ind, size_functions[sfunc_name], community_detection[alg_name], num_of_runs, num_of_rand_nets, p_values, nhat, qhat, rgindex);
       	cout <<"   end"<<endl<<endl;
    }
   
	

    /* Save results */
    cout << "   Saving..."<<endl;
    writeLabels(outputfile, xlist, p_values, corrected_alpha);
    cout <<"   end"<<endl;
    cout << "===================" << endl;

    return 0;
}

void usage()
{

    cout << endl
         << "\e[1mNAME:\e[0m" << endl
         << endl;
    cout << "	gp - Graph partitioning algorithms" << endl
         << endl
         << endl;

    cout << "\e[1mSYNOPSIS:\e[0m" << endl
         << endl;
    cout << "	\e[1m./gp\e[0m [input-file] [output-file] [options]" << endl
         << endl
         << endl;
    cout << "\e[1mDESCRIPTION:\e[0m" << endl
         << endl;
    cout << "	\e[1m./gp\e[0m seeks non-overlapping communities in the network given by [input-file] and" << endl;
    cout << "	saves the detected communities in [output-file]." << endl
         << endl;

    cout << "	\e[1m[input-file]\e[0m" << endl;
    cout << "	    The file is the list of all pairs of adjacent nodes (tab-separated value file)." << endl;
    cout << "	    The first and second columns are the IDs of a node and the adjacent node, respectively." << endl;
    cout << "	    The third column is the weight of the edge." << endl;
    cout << "	    The node's ID is assumed to start from 1." << endl
         << endl;

    cout << "	\e[1m[output_file]\e[0m" << endl;
    cout << "	    This file describes the detected communities (tab-separated value file)." << endl;
    cout << "	    The first column is the ID of each node." << endl;
    cout << "	    The second column is the index of the community to which each node belongs." << endl
         << endl
         << endl;
    cout << "\e[1mOPTIONS:\e[0m" << endl
         << endl;
    cout << "	\e[1m-r=[R]\e[0m  Run the algorithm R times. (Default: 10)" << endl
         << endl
         << "	\e[1m-o=[ALG]\e[0m Specify one of the following optimisation algorithms" << endl
         << "	    - kl: Kernighan-Lin algorithm" << endl
         << "	    - ls: Label switching algorithm" << endl
         << "	    - mcmc: Markov chain monte carlo algorithm" << endl
         << "	    - louvain: Louvain algorithm" << endl
         << endl\
         << "	\e[1m-q=[Q]\e[0m Specify one of the following quality functions:" << endl
         << "	    - qint: average degree of each community" << endl
         << "	    - qext: ratio cut" << endl
         << "	    - qcnd: normalised cut" << endl
         << "	    - qmod: modularity" << endl
         << "	    - dcsbm: dcSBM" << endl
         << endl
         << "	\e[1m-k=[K]\e[0m  Set the number of communities to K. (Default: 2)" << endl
         << endl
         << "	\e[1m-a=[ALPHA]\e[0m  Set significance level ALPHA. (Default: 0.05. Set 1 to disable the statistical test)" << endl
         << endl
         << "	\e[1m-l=[NUM]\e[0m  Set the number of randomised networks to NUM. (Default: 500)" << endl
         << endl
         << "	\e[1m-s=[S]\e[0m Specify one of the following size functions:" << endl
         << "	    - nodes: number of nodes in a community" << endl
         << "	    - edges: number of edges incident to a community" << endl
         << endl;
}

void split(const string& s, char c,
    vector<string>& v)
{
    string::size_type i = 0;
    string::size_type j = s.find(c);

    while (j != string::npos) {
        v.push_back(s.substr(i, j - i));
        i = ++j;
        j = s.find(c, j);

        if (j == string::npos)
            v.push_back(s.substr(i, s.length()));
    }
}

void readEdgeTable(string filename, vector<vector<int>>& A, vector<vector<double>>& W, int& N)
{

    std::ifstream ifs(filename);
    vector<int> edgeList;
    vector<double> wList;
    string str;
    N = 0;
    edgeList.clear();
    while (getline(ifs, str)) {
        vector<string> v;
        split(str, delimiter, v);

        int sid = stoi(v[0]) - 1;
        int did = stoi(v[1]) - 1;
	double w = 1;
	
	if(v.size()>2){
        	w = stof(v[2]);
	}

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

void writeLabels(const string filename, const vector<vector<bool>>& xlist, const vector<double>p_values, double alpha)
{
    int K = xlist.size();
    int N = xlist[0].size();
    vector<int> C(N);
    for (int i = 0; i < N; i ++) {
    for (int k = 0; k < K; k ++) {
	if(xlist[k][i]){
        	C[i] = k;
		break;
	}
    }
    }
    FILE* fid = fopen(filename.c_str(), "w");
    fprintf(fid, "id%cc%cs%cpval", delimiter, delimiter, delimiter);
    for (int i = 0; i < N; i++) {
        if(p_values[C[i]] < alpha){
        	fprintf(fid, "\n%d%c%d%c%d%c%f", i + 1, delimiter, C[i] + 1, delimiter, 1, delimiter, p_values[C[i]]);
	}else{
        	fprintf(fid, "\n%d%c%d%c%d%c%f", i + 1, delimiter, C[i] + 1, delimiter, 0, delimiter, p_values[C[i]]);
	}
    }
    fclose(fid);
}
