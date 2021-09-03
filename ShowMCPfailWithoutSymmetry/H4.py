import sys
import vqe_methods_add_by_one
import numpy
import pickle
import numpy as np

dist = float(sys.argv[1])
print(dist)
if __name__ == '__main__':
    

    PoolOfHope=['YYZZXYII', 'IXXZIXYI', 'IXXIZZXY', 'XYZZXXII', 'XXZIIIXY', 'ZZYYIZXY', 'XZIYXXII', 'XZIZZZZY', 'IYYXIIZY', 'ZXZZIIZY', 'XXYXYIZY', 'XIYYXIYI', 'IXYYZZXY', 'YZIIYIYI']
    Resultat=[]
        
        
    
        
    geometry = [('H', (0, 0, 0)), ('H', (0, 0, dist)), ('H', (0, 0, 2*dist)), ('H', (0, 0, 3*dist))]
    print(geometry)
    vqe_methods_add_by_one.adapt_vqe(geometry,
    	                  adapt_thresh    = 1e-8,                        #gradient threshold
                          adapt_maxiter   = 400,                       #maximum number of ops                   
                          Pool            = PoolOfHope,
                          Resultat        = Resultat,
                          bond_legth      = dist
                          ) 
            
    with open('Bond_length_dependence.LiH_dissociation_curve_pickle_min_pool_{}'.format(dist) , 'wb') as handle:
        pickle.dump(Resultat, handle, protocol=pickle.HIGHEST_PROTOCOL)                        

