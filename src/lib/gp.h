#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>
#include <omp.h>

using namespace std;

#include "quality_functions.h"
#include "community-detection-algorithms/comalgorithms.h"
#include "qstest/qstest.cpp"

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

