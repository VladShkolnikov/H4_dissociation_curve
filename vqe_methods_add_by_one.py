import scipy
import math

import pickle
import numpy as np

# import cirq
# import openfermioncirq
# from openfermioncirq import trotter
import scipy.sparse.linalg

import operator_pools
from tVQE import tUCCSD

import openfermion, openfermionpsi4
from openfermion import *
from pyscf_backend import *
from of_translator import *

QubitNumber=8



def HartreeFock(molecule):

    #exps=rotate_me()
    
    hamiltonian_op = molecule.get_molecular_hamiltonian()  
    

    hamiltonian_f = openfermion.transforms.get_fermion_operator(hamiltonian_op)
    hamiltonian_jw = openfermion.transforms.jordan_wigner(hamiltonian_f)
    hamiltonian = openfermion.transforms.get_sparse_operator(hamiltonian_jw)
    reference_ket = scipy.sparse.csc_matrix(
    openfermion.jw_configuration_state(
    list(range(0,molecule.n_electrons)), molecule.n_qubits)).transpose()


    
    
    return reference_ket, hamiltonian

def adapt_vqe(geometry,
        adapt_conver    = 'norm',
        adapt_thresh    = 1e-7,
        theta_thresh    = 1e-10,
        adapt_maxiter   = 400,
        Pool            =['YIIIIIII'],
        Resultat        =[], 
        bond_legth      =0

        ):

    molecule = openfermion.hamiltonians.MolecularData(geometry, "sto-3g", 1)
    molecule = openfermionpsi4.run_psi4(molecule, 
                run_scf = 1, 
                run_mp2=0, 
                run_cisd=0, 
                run_ccsd = 0, 
                run_fci=1, 
                delete_input=1)

    print(" adapt threshold:", adapt_thresh)
    print(" theta threshold:", theta_thresh)
    print(" initial state: HF state")

    print('molecule.n_atoms:',molecule.n_atoms)
    print('molecule.n_electrons:',molecule.n_electrons)
    print('molecule.protons:',molecule.protons)
    print('molecule.nuclear_repulsion:',molecule.nuclear_repulsion)
    print('molecule.n_orbitals:',molecule.n_orbitals)
    print('molecule.n_qubitss:',molecule.n_qubits)
    print('molecule.orbital_energies:',molecule.orbital_energies)
    
    #pool1=operator_pools.OperatorPool(molecule,Pool[0])
    #pool2=operator_pools.OperatorPool(molecule,Pool[1])
    pool=operator_pools.OperatorPool(Pool)

    print(' HF energy      %20.16f au' %(molecule.hf_energy))
    #print(' CISD energy    %20.16f au' %(molecule.cisd_energy))
    #print(' CCSD energy    %20.16f au' %(molecule.ccsd_energy))
    print(' FCI energy     %20.16f au' %(molecule.fci_energy))

    
    reference_ket, hamiltonian=HartreeFock(molecule)
    












    w, v = scipy.sparse.linalg.eigs(hamiltonian,k=1, which='SR')   
    GSE = min(w).real
    print('Ground state energy:', GSE)

    #Thetas
    parameters = []
    ansatz_ops = []     #SQ operator strings in the ansatz
    ansatz_mat = []     #Sparse Matrices for operators in ansatz
    

    E = reference_ket.transpose().conj().dot(hamiltonian.dot(reference_ket))[0,0].real
    print('initial energy', E)

    print(" Start ADAPT-VQE algorithm")

    curr_state = 1.0*reference_ket

    fermi_ops = pool.fermi_ops
    spmat_ops = pool.spmat_ops
    n_ops = pool.n_ops


    error=[]
    trial_model=None
    curr_energy=float(100)
    print(" Now start to grow the ansatz")
    flag=0
    for n_iter in range(0,adapt_maxiter):
        
    
        print("\n\n\n")
        print(" --------------------------------------------------------------------------")
        print("                         ADAPT-VQE iteration: ", n_iter)                 
        print(" --------------------------------------------------------------------------")

        min_options = {'gtol': theta_thresh, 'disp':False}
        
        if flag==0:
            max_derivative=-math.inf
            trial_model = tUCCSD(hamiltonian, ansatz_mat, reference_ket, parameters)
            for op in range(n_ops):
                pos=0
                der=trial_model.derivative(parameters,pos,spmat_ops[op])
                #print(der)
                if der > max_derivative:
                    max_derivative=der
                    ansatz_ops_res=ansatz_ops.copy()
                    ansatz_mat_res=ansatz_mat.copy()
                    parameters_res=parameters.copy()
                    ansatz_ops_res.insert(pos,Pool[op])
                    ansatz_mat_res.insert(pos,spmat_ops[op])
                    parameters_res.insert(pos,0)
                        
                        
       
            
            ansatz_mat=ansatz_mat_res
            ansatz_ops=ansatz_ops_res
            parameters=parameters_res
            trial_model = tUCCSD(hamiltonian, ansatz_mat, reference_ket, parameters)
            opt_result = scipy.optimize.minimize(trial_model.energy, parameters, jac=trial_model.gradient, 
            options = min_options, method = 'BFGS', callback=trial_model.callback)
            
            parameters=list(opt_result['x'])
            print(ansatz_ops)
            
            curr_state = trial_model.prepare_state(parameters)
            curr_energy = trial_model.energy(parameters)
            if  abs(curr_energy-GSE)>adapt_thresh:
                #print(" new state: ",curr_state)
                print(" Finished: %20.12f" % curr_energy)
                print(" Error: %20.12f" % abs(curr_energy-GSE))
                print(bond_legth)
                error.append(abs(curr_energy-GSE))
                continue
            else:
                print(" Finished: %20.12f" % curr_energy)
                print(" Error: %20.12f" % abs(curr_energy-GSE))
                error.append(abs(curr_energy-GSE))
                flag=1
                break

    Resultat.append({'bond_length:':bond_legth, 'Pool:':Pool,'ansatz:':ansatz_ops,'parameters:':parameters,'final_error:': abs(curr_energy-GSE),'GSE:':GSE})
    print('\n','\n','------------------------------------------','\n')
                        
                        
            
        

