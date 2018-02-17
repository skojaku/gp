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
std::mt19937_64 mtrnd;

void init_random_number_generator();
