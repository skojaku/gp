import csv
from scipy import sparse
import numpy as np
import gp as g


linkfilename='links_karate.dat'
edges = np.genfromtxt(linkfilename, delimiter='\t', skip_header = 0)
N = int(np.max(np.amax(edges[:,[0,1]], 0))) # number of nodes

communities = g.detect(edges, algorithm = 'kl', qfunc = 'dcsbm', significance_level = 0.05)

print(communities)
