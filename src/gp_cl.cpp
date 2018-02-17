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
#include "gp.h"
#include "quality_functions.h"
#include "comalgorithms.h"

#include <string>
#include <stdio.h>
#include <unistd.h>

char delimiter = '\t';

void split(const string& s, char c,
    vector<string>& v);

void readEdgeTable(string filename, vector<vector<int> >& A, vector<vector<double>>& W, int& N);

void writeLabels(const string filename, const vector<vector<bool>>& xlist);

void usage();

// initialise mtrnd
void init_random_number_generator()
{
    random_device r;
    seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    mtrnd.seed(seed);

    /* Use the following code if you cannot initialise mtrnd with random_device */
    /*
	int seeds[624];
	size_t size = 624*4; //Declare size of data
	std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary); //Open stream
	if (urandom) //Check if stream is open
	{
	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
	    urandom.close(); //close stream
	}
	else //Open failed
	{
	    		std::cerr << "Failed to open /dev/urandom" << std::endl;
	}
	std::seed_seq seed(&seeds[0], &seeds[624]);
	mtrnd.seed(seed);
*/
}

int main(int argc, char* argv[])
{
    if (argc <= 2) {
        usage();
        return 1;
    }

    string linkfile = argv[1];
    string outputfile = argv[2];

    int num_of_runs = 10;
    int K = 2;

    int opt;
    opterr = 0;
    string tmp;
    string alg_name = "kl";
    string qfunc_name = "dcsbm";
    while ((opt = getopt(argc, argv, "h:k:r:a:q:")) != -1) {
        switch (opt) {
        case 'h':
            usage();
            break;
        case 'a':
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
    cout << "# algorithm: "<< alg_name<< endl;
    cout << "# quality function: "<< qfunc_name<< endl;
    cout << "# Number of runs: "<< num_of_runs<< endl;
    cout << "# Number of communities: "<< K<< endl;
    cout << "# input file:  "<< linkfile<< endl;
    cout << "# output file: "<< outputfile<< endl;
    cout << "" << endl;
  
    cout << "# Status: "<< endl;
    cout << "    Reading input file...";
    /* Read file */
    int N;
    vector<vector<int> > A;
    vector<vector<double>> W;
    readEdgeTable(linkfile, A, W, N);
    cout <<"end"<<endl;
      

    /* Detect communities in networks */
    cout << "    Seeking communities...";
    init_random_number_generator();
    vector<vector<bool>> xlist;
    double Qr = -numeric_limits<double>::max();
    mcmc_qfunc = quality_functions[qfunc_name];
    mcmc_qfunc_diff = quality_functions_diff[qfunc_name];
    for (int r = 0; r < num_of_runs; r++) {
        vector<vector<bool>> xlist_tmp(K, vector<bool>(N, false) );
    	community_detection[alg_name](A, W, xlist_tmp);
        double Qi = quality_functions[qfunc_name](A, W, xlist_tmp);
        if(Qi != Qi) {continue;}
        
        if (Qi > Qr) {
            Qr = Qi;
            xlist = xlist_tmp;
        }
    }
    cout <<"end"<<endl;

    /* Save results */
    cout << "    Saving...";
    writeLabels(outputfile, xlist);
    cout <<"end"<<endl;
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
    cout << "	    This file describes the detected core-periphery pairs (tab-separated value file)." << endl;
    cout << "	    The first column is the ID of each node." << endl;
    cout << "	    The second column is the index of the community to which each node belongs." << endl
         << endl
         << endl;
    cout << "\e[1mOPTIONS:\e[0m" << endl
         << endl;
    cout << "	\e[1m-r=[R]\e[0m  Run the algorithm R times. (Default: 10)" << endl
         << endl
         << "	\e[1m-a=[ALG]\e[0m Specify one of the following algorithms" << endl
         << "	    - kl: Kernighan-Lin algorithm" << endl
         << "	    - mcmc: Markov chain monte carlo algorithm" << endl
         << endl\
         << "	\e[1m-a=[Q]\e[0m Specify one of the following quality functions:" << endl
         << "	    - qint: average degree of each community" << endl
         << "	    - qext: ratio cut" << endl
         << "	    - qcnd: normalised cut" << endl
         << "	    - qmod: modularity" << endl
         << "	    - dcsbm: dcSBM" << endl
         << endl
         << "	\e[1m-k=[K]\e[0m  Set the number of communities to K. (Default: 2)" << endl
         << endl
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
    for (int i = 0; i < edgeList.size(); i += 2) {
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

void writeLabels(const string filename, const vector<vector<bool>>& xlist)
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
    fprintf(fid, "id%cc\n", delimiter);
    for (int i = 0; i < N; i++) {
        fprintf(fid, "%d%c%d\n", i + 1, delimiter, C[i] + 1);
    }
    fclose(fid);
}
