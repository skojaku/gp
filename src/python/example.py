import csv
import numpy as np
import gp as g


linkfilename='links_karate.dat'
edges = np.genfromtxt(linkfilename, delimiter='\t', skip_header = 0)

communities = g.detect(edges, K = 3, algorithm = 'mcmc', qfunc = 'mod', significance_level = 1.00)

print(communities)
