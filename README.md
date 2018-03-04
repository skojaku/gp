# gp
C++, Matlab and Python code for community-detection algorithms in networks.

The following optimisation algorithms for community detection are implemented:
 * Kernighan-Lin algorithm \[[1](https://en.wikipedia.org/wiki/Kernighan%E2%80%93Lin_algorithm), [2](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.016107)\] 
 * Label switching algorithm \[[3](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.036106)\] 
 * Markov chain monte carlo algorithm \[[4](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)\] 
 * Louvain algorithm \[[5](http://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta)\]

These algorithms seek communities in networks by optimising the following quality functions: 
 * Internal average degree \[[6](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Ratio cut criterion \[[6](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Normalised cut criterion \[[6](http://www.tandfonline.com/doi/abs/10.1080/15427951.2009.10129177)\]
 * Modularity \[[7](http://www.pnas.org/content/103/23/8577)\]
 * Log likelihood for the degree-corrected stochastic block model \[[3](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.016107)\]

Some detected communities may be insignificant.
Indeed, the community-detection algorithms detect communities even if there is no community in networks.
Therefore, it is important to check the significance of each detected community.
I implemented the (q,s)-test \[[8](https://arxiv.org/abs/1712.00298)\] to compute the statistical significance of the detected communities.  

Table of Contents
=================

* [C\+\+](#c)
* [Matlab](#matlab)
* [Python](#python)
* [Requirements](#requirements)


## C++ 

### Compile 

This package contains a command-line client.
To compile, run 

```bash 
make cpp 
```

This creates the command-line client ''gp'' in the current directory.
To test, run ``./gp''

### Usage 
 
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
     * 'ls' - Label Switching algorithm 
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
  
### Example (src/cpp/example.sh)
  
```bash
./gp links_karate.dat result.dat -o louvain -a 0.01 -q mod 
./gp links_karate.dat result.dat -o kl -a 0.01 -q dcsbm -k 5 
```


## Matlab

### Compile 

To compile, run 

```bash 
make matlab
```

This creates a mex file, gp_mex.mexa64 (the file extension may be different depending on OS).

Copy src/matlab/gp.m, src/matlab/gp_mex.mexa64 and src/matlab/qstest_mex.mexa64 to your working directory. 

You may have the following message:

```bash
Warning: You are using gcc version '5.4.0'. The version of gcc is not supported. 
The version currently supported with MEX is '4.9.x'. 
```

This means you are required to change the version of g++ compiler. 
To remedy this, modify a line in ''./src/matlab/makefile'' as follows: 

```bash
MEXCOMPILER := g++-(the version compatible with your mex compiler, e.g., g++-4.9) 
```

See https://uk.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html for detail.

### Usage 

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
     * ls - Label Switching algorithm 
     * mcmc - Markov chain monte carlo 
     * louvain - Louvain algorithm 
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
  
### Example src/matlab/example.m
  
```Matlab
T = readtable('links_karate.dat', 'Delimiter', '\t', 'HeaderLines',0);
T = table2array(T);
N = max([max(T(:,2)),max(T(:,1))]);
A = sparse(T(:,1), T(:,2), T(:,3), N, N);

g = gp();
param = g.init(); % initialise
cids = g.detect(A, param);
```


## Python

### Compile

To compile, run 

```bash 
make python
```

This creates a shared library ''src/python/gp.cpython-36m-x86_64-linux-gnu.so'' callable from python. 
Copy the shared library to your working directory. 

### Usage
 
```python
import gp as g
communities = g.detect(edges, K, qfunc, algorithm, num_of_runs, significance_level, sfunc, num_of_rand_nets)
```
 
#### Input 

 * `edges` - Mx3 Numpy array, where M is the number of edges. The first and second columns indicate the IDs of nodes connected by an edge. The third column indicates the weight of the edge.
   * `K` - Number of communities. Default K=2. 
   * `qfunc` - Quality function for communities. Default qfunc= 'dcsbm'. The following quality functions are available:
     * mod - Contribution of a community to the modularity 
     * int - Internal average degree 
     * rcut - Expansion　
     * ncut - Conductance
     * dcsbm - Degree-corrected stochastic block model
   * `algorithm` - Optimisation algorithm. Default algorithm = 'kl'. The following algorithms are available: 
     * kl - Kernighan-Lin algorithm 
     * ls - Label Switching algorithm 
     * mcmc - Markov chain monte carlo 
     * louvain - Louvain algorithm 
   * `num_of_runs` - Number of runs. Default num_of_runs = 1. 
   * `significance_level` - Significance level. Default significance_level = 0.05.
   * `sfunc` - Size of community. The following size functions are available:
     * nodes - Number of nodes in a community 
     * edges - Number of edges incident to a community 
   * `num_of_rand_nets` - Number of randomised networks. Default num_of_rand_nets = 500. 
　
  
#### Output 

 * `communities` - List of length 3. 
   * communities[0] - Numpy array of length N, where N is the number of nodes. communities[0][i] indicates the index of the community to which node i belongs.
   * communities[1] - Numpy array of length N. communities[1][i] indicates the p-value of the community to which node i belongs.
   * communities[2] - Numpy array of length N. communities[2][i] = True or False indicates the significant or insignificant communities, respectively.
  
### Example src/python/example.py
  
```python
import csv
import numpy as np
import gp as g


linkfilename='links_karate.dat'
edges = np.genfromtxt(linkfilename, delimiter='\t', skip_header = 0)

communities = g.detect(edges, algorithm = 'kl', qfunc = 'dcsbm', significance_level = 0.05)

print(communities)
```

## Requirements

 * Matlab 2012 or later 
 * Python3.4 or later
 * Cmake2.8 or later

