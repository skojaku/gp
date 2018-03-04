import csv
import numpy as np
import gp as g


linkfilename='links_karate.dat'
edges = np.genfromtxt(linkfilename, delimiter='\t', skip_header = 0)

communities = g.detect(edges, algorithm = 'kl', qfunc = 'dcsbm', significance_level = 0.05)

print(communities)
