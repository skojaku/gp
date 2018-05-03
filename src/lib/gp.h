#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif


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


// initialise mtrnd
mt19937_64 init_random_number_generator()
{
    mt19937_64 mtrnd;
    random_device r;
    seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    mtrnd.seed(seed);
    return mtrnd;

/* Use the following code if you cannot initialise mtrnd with random_device */
/*
        mt19937_64 mtrnd;
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
	return mtrnd;
*/
}
