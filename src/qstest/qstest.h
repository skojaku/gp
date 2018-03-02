/*
*
* Header file of the Kojaku-Masuda algorithm (C++ version)
*
*
* An algorithm for finding multiple core-periphery pairs in networks
* "???"
* Sadamori Kojaku and Naoki Masuda
* Preprint arXiv:???
* 
*
* Please do not distribute without contacting the authors above.
* If you find a bug in this code, please contact the authors.
*
*
* AUTHOR - Sadamori Kojaku
*
*
* DATE - 11 Oct 2017
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>
#ifdef _OEPNMP
#include <omp.h>
#endif

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
using namespace std;

/* Global variables */
//std::mt19937_64 mtrnd;

/* ---- Estimating statistical significance of a community ----*/
void estimate_statistical_significance(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<vector<bool>>& xlist,
    double (*calc_s)(const vector<vector<int>>&, const vector<bool>&),
    double (*calc_q)(const vector<vector<int>>&, const vector<bool>&),
    void (*find_communities)(const vector<vector<int> >&, vector<vector<bool>>&),
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values,
    vector<int>& nhat,
    vector<double>& qhat,
    vector<int>& rgindex);

/* 
*
* Impliment the following functions
*
* */
mt19937_64 init_random_number_generator();
