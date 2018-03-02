# gp
C++ and MATLAB codes for community-detetion algorithms in networks.

I implemented the following optimisation algorithms for community detection:
 * Kernighan-Lin algorithm \[[1](https://en.wikipedia.org/wiki/Kernighan%E2%80%93Lin_algorithm), [2](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.016107)\] 
 * Markov chain monte carlo algorithm \[[3](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)\] 
 * Louvain algorithm \[[4](http://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta)\]

These algorithms seek communities in networks by optimising the following quality functions: 
 * Internal average degree \[[5](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Ratio cut criterion \[[5](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Normalised cut criterion \[[5](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Modularity \[[6](http://www.pnas.org/content/103/23/8577)\]
 * Log likelihood for the degree-corrected stochastic block model \[[2](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.016107)\]

Some detected communities may be insignificant.
Indeed, the community-detection algorithms detect communities even if there is no community in networks.
Therefore, it is important to check the significance of each detected community.
I implemented the (q,s)-test \[[7](https://arxiv.org/abs/1712.00298)\] to compute the statistical significance of the detected communities.  

# Installation

To install, run 

```bash 
make 
```


This creates the following two files:
 * gp - command-line crient
 * gp_mex.mexa64 - mex file (the file extension may be different depending on OS)

For MATLAB user, copy gp.m, gp_mex.mexa64 and qstest_mex.mexa64 to your working directly. 

You may have the following message:

```bash
Warning: You are using gcc version '5.4.0'. The version of gcc is not supported. 
The version currently supported with MEX is '4.9.x'. 
```

This means you are required to change the version of g++ compiler. 
To remedy this, modify a line in ''./makefile'' as follows: 

```bash
MEXCOMPILER := g++-(the version compatible with your mex compiler, e.g., g++-4.9) 
```

See https://uk.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html for detail.

# Usage

## Matlab
 
```Matlab
g = gp(); 
param = g.init(); % initialise
[cids, qs, h, pvals] = g.detect(adjmat, param) % Community detection
```
 
#### Input 

 * `adjmat` - Adjacency matrix
 * `param` - Configurations of community detection and statistical test. 
   * `param.qfunc` - Quality function for communities. Default qfunc= 'dcsbm'. The following quality functions are available:
     * mod - Contribution of a community to the modularity 
     * int - Internal average degree 
     * rcut - Expansion　
     * ncut - Conductance
     * dcsbm - Degree-corrected stochastic block model
   * `param.algorithm` - Optimisation algorithm. Default algorithm = 'kl'. The following algorithms are available: 
     * kl - Kernighan-Lin algorithm 
     * mcmc - Markov chain monte carlo 
     * louvain - Markov chain monte carlo 
   * `param.K` - Number of communities. Default K=2. 
   * `param.num_of_runs` - Number of runs. Default num_of_runs = 1. 
   * `param.significance_level` - Significance level. Default significance_level = 0.05.
   * `param.sfunc` - Size of community. The following size functions are available:
     * nodes - Number of nodes in a community 
     * edges - Number of edges incident to a community 
   * `param.num_of_rand_nets` - Number of randomised networks. Default num_of_rand_nets = 500. 
　
  
#### Output 

 * `cids` - Column vector of length N, where N is the number of nodes. cids[i] is the index of the community to which node i belongs. 
 * `qs` - Column vector of length K, where qs[k] is the quality value for the kth community. 
 * `h` - Column vector of length N, where h[i] = 1 or h[i] = 0 indicates that the node i belongs to the significant community, respectively.  
 * `pvals` - Column vector of length K, pvals[k] is the p-value of the kth community.  
  
#### Example example.m
  
```Matlab
T = readtable('links_karate.dat', 'Delimiter', '\t', 'HeaderLines',0);
T = table2array(T);
N = max([max(T(:,2)),max(T(:,1))]);
A = sparse(T(:,1), T(:,2), T(:,3), N, N);

g = gp();
param = g.init(); % initialise
cids = g.detect(A, param);
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
   * -r=[R] Run the algorithm R times. (Default: 10)
   * -l=[A] Specify one of the following algorithms. (Default: 'kl')
     * 'kl': Kernighan-Lin algorithm
     * 'mcmc': Markov chain monte carlo algorithm
     * 'louvain': Louvain algorithm 
   * -q: Quality function. (Default: 'dcsbm') 
            - 'int': Internal average degree 
            - 'rcut': Expansion　
            - 'ncut': Conductance
	    - 'mod': Modularity
	    - 'dcsbm': dcSBM
   * -a=[ALPHA] Set significance level ALPHA. (Default: 0.05. Set 1 to disable the statistical test)
   * -l=[NUM] Set the number of randomised networks to NUM. (Default: 500)
   * -s=[S] Specify one of the following size functions:
	    - 'nodes': number of nodes in a community
	    - 'edges': number of edges incident to a community
  
#### Example
  
```bash
./gp links_karate.dat result.dat -k 3 -a 'kl' -q 'dcsbm' -r 1 -a 0.05 
```


# Requirements

 * Matlab 2012 or later 

