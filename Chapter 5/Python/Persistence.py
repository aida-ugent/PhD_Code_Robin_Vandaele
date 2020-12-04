### library imports ###

from ripser import ripser
import sys
import os
import json
import codecs
import numpy as np
import matplotlib.pyplot as plt
import time



# Extract arguments

path = sys.argv[1]

with open(os.path.join(path,'arguments.json'), 'r') as f:
    arguments = json.load(f)
    
X = np.empty([len(arguments['X']), len(arguments['X'][0])])
dims = list(arguments['X'][0].keys())
for i in range(len(arguments['X'])):
    X[i,:] = [arguments['X'][i][dim] for dim in dims]
maxdim = arguments['maxdim']
thresh = np.Inf if arguments['thresh'] == 'Inf' else arguments['thresh']
coeff = arguments['coeff']
distance_matrix = arguments['distance_matrix']
do_cocycles = arguments['do_cocycles']
metric = arguments['metric']
n_perm = arguments['n_perm']

del arguments



### compute persistence ###
 
tic = time.time()
result = ripser(X, 
                maxdim = maxdim, 
                thresh = thresh, 
                coeff = coeff,
                distance_matrix = distance_matrix,
                do_cocycles = do_cocycles,
                metric = metric,
                n_perm = n_perm)
toc = time.time()



### store the results ###
    
result_json = json.dumps({'dgms': [dgm.tolist() for dgm in result['dgms']],
                          'cocycles': [[cc.tolist() for cc in ccs] for ccs in result['cocycles']],
                          'num_edges': result['num_edges'],
                          'dperm2all': result['dperm2all'].tolist(),
                          'idx_perm': result['idx_perm'].tolist(),
                          'r_cover': result['r_cover'],
                          'elapsed_time':  toc-tic})

with open(os.path.join(path,'results.json'),"w") as json_file:
    json_file.write(result_json)


