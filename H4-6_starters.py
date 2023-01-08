import sys
from typing import List, Union, Tuple

import vqe_methods_add_by_one
import numpy
import pickle
import numpy as np

dist = float(0.9)
print(dist)
if __name__ == '__main__':
    

    six_starters_pool=['YIIXYZZY', 'IXXZXIZY', 'ZZXYYYII', 'YZZYXZZY', 'IXXZIXYI', 'YYZIIIXY', 'ZIXIIZYI', 'IYYYIYXY', 'YIZIXZII', 'ZZXIZZYI', 'ZXZXXIYI']
    Resultat=[]
    geometry = [('H', (0, 0, 0)), ('H', (0, 0, dist)), ('H', (0, 0, 2*dist)), ('H', (0, 0, 3*dist))]
    print(geometry)
    vqe_methods_add_by_one.adapt_vqe(geometry,
    	                  adapt_thresh    = 1e-8,                        #gradient threshold
                          adapt_maxiter   = 400,                       #maximum number of ops                   
                          Pool            = six_starters_pool,
                          Resultat        = Resultat,
                          bond_legth      = dist
                          ) 
            
    with open('Bond_length_dependence.six_starters_{}'.format(dist) , 'wb') as handle:
        pickle.dump(Resultat, handle, protocol=pickle.HIGHEST_PROTOCOL)                        

