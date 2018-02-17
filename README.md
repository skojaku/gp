# gp
C++ and MATLAB codes for community detetion algorithms in networks.

Available algorithms
 * Kernighan-Lin algorithm  
 * Markov chain monte carlo 

Available quality functions
 * Internal average degree 
 * Ratio cut criterion
 * Normalised cut criterion
 * Modularity 
 * Degree-corrected stochastic block model 

# Installation

To install, run 

```bash 
make 
```


This creates the following two files:
 * gp - command-line crient
 * gp_mex.mexa64 - mex file (the file extension may be different depending on OS)

For MATLAB user, copy gp.m and gp_mex.mexa64 to your working directly. 

# Usage

## Matlab
 
```Matlab
cids = gp(adjmat, qfunc='dcsbm', K=2, num_of_runs=1, algorithm='kl')
```
 
#### Input 

 * `adjmat` - Adjacency matrix 
 * `K` (optional) - Number of communities. Default K=2. 
 * `qfunc` (optional) - Quality function for communities. Default qfunc= 'dcsbm'. The following quality functions are available:
   * qmod - Contribution of a community to the modularity 
   * qint - Internal average degree 
   * qexp - Expansion　
   * qcnd - Conductance
   * dcsbm - Degree-corrected stochastic block model
 * `num_of_runs` (optional) - Number of runs. Default num_of_runs = 1. 
 * `algorithm` (optional) - Optimisation algorithm. Default algorithm = 'kl'. The following algorithms are available: 
   * kl - Kernighan-Lin algorithm 
   * mcmc - Markov chain monte carlo 
　
  
#### Output 

 * `cids` - Column vector of length N, where N is the number of nodes. cids[i] is the index of the community to which node i belongs. 
  
#### Example example.m
  
```Matlab
T = readtable('links_karate.dat', 'Delimiter', '\t', 'HeaderLines',0);
T = table2array(T);
N = max([max(T(:,2)),max(T(:,1))]);
A = sparse(T(:,1), T(:,2), T(:,3), N,N);
[r,c,v] = find(triu(A,1));

cids = gp(A, 'dcsbm',2)
```

## C++
 
``` c++
./gp [input-file] [output-file] [options]
```
 
./gp seeks non-overlapping communities in the network given by [input-file] and saves the detected communities in [output-file].

#### Input 
 
 * `[inputfile]` - List of edges 
   * The file is the list of all pairs of adjacent nodes (tab-separated value file).
   * The first and second columns are the IDs of a node and the adjacent node, respectively.
   * The third column is the weight of the edge.
   * The node's ID is assumed to start from 1.
 * `[outputfile]` - Output file 
   * This file describes the detected communities (tab-separated value file).
   * The first column is the ID of each node.
   * The second column is the index of the community to which each node belongs.
 * `[options]` - Options 
   * -k=[K] Set the number of communities to K. (Default: 2)
   * -r=[R] - Run the algorithm R times. (Default: 10)
   * -a=[A] - Specify one of the following algorithms. (Default: 'kl')
     * 'kl': Kernighan-Lin algorithm
     * 'mcmc': Markov chain monte carlo algorithm
   * q: Quality function. (Default: 'dcsbm') 
	    - 'qint': average degree of each community
	    - 'qext': ratio cut
	    - 'qcnd': normalised cut
	    - 'qmod': modularity
	    - 'dcsbm': dcSBM
  
#### Example
  
```bash
./gp links_karate.dat result.dat -k 3 -a 'kl' -q 'dcsbm' -r 1
```


# Requirements

 * g++4.9 or later 
 * Matlab 2012 or later 

