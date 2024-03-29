from pyscf_backend import *
from of_translator import *

#Geometry in angstroms
geometry = "Li 0 0 0; H 0 0 1"

#Basis set has to be something recognized by PySCF.  (sto-3g, 6-31g, etc.)
basis = "sto-3g"

#I only trust RHF right now- I can debug open-shell refs if you need them though.
reference = "rhf"

#Number of frozen electrons.  Probably an even number.
frozen_core = 2

#Number of frozen virtual orbitals.
frozen_vir = 2


'''
In this case you have reference:

|111100000000>

becoming

Frozen |11> + Active |11000000> + Frozen |00>
= |core> + |active> + |vir>

where you've frozen the 2 electrons in the lowest energy orbitals,
and simply removed the highest 2 spin-orbitals from consideration.
For RHF, this is equivalent to freezing the highest and lowest spatial orbitals.
'''

#Use PySCF to compute spatial orbital integrals, then put them into inefficient spin-orbital tensors that OpenFermion can read
_0body_H, _1body_H, _2body_H, _1rdm, hf_energy = get_integrals(geometry, basis, reference, frozen_core = frozen_core, frozen_vir = frozen_vir)

#_0body_H is the nuclear repulsion + <core|H|core>
#_1body_H is the 1-body part of <active|H|active> + <core|H|active> 
#_2body_H is the 2-body part of <core|H|core>
#_1rdm is the MO basis rdm.
#hf_energy is <ref|H|ref> and doesn't care about your orbital freezing
#These are all in the weird, proprietary format of pyscf's tensor construction

#Get number of electrons:
N_e = int(np.trace(_1rdm))
print("How many: " ,N_e)
#Use OpenFermion to build sparse matrix representations of JW-transformed operators:
H, ref, N_qubits, S2, Sz, Nop = of_from_arrays(_0body_H, _1body_H, _2body_H, N_e)
print(H)
#, ref, N_qubits, S2, Sz, Nop)
#H is the Hamiltonian, ref is the HF reference in the active space.  S2, S_z, and N_op are the S^2, S_z, and number operators if you happen to want them.

ci = np.linalg.eigh(H.todense())[0][0]
print(f"Problem reduced to {N_qubits} qubits.")
print(f"SCF energy:                {hf_energy:20.16f}")
print(f"Active space FCI energy:   {ci:20.16f}")


