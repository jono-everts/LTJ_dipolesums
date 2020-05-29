import numpy as np
import DipoleLatticeField
import pickle
import gzip


def tup_to_array(tup):
    mat = np.array([[tup[0], tup[1], tup[2]],[tup[3], tup[4], tup[5]],[tup[6], tup[7], tup[8]]])
    return mat

#Define SI units
pi = np.pi
mB=9.274*10**(-24)
k=1.380*10**(-23)
NA=6.022*10**23
mu0=4*pi*10**(-7)

NdGaO3 = DipoleLatticeField.MonoclinicLattice()
NdGaO3.axes(5.4276, 5.4979, 7.7078, 0)
NdGaO3.g_tensor(2.7, 2.02,0)

#Define the ions to be in the standard two sublattice antiferromagnetic ground state, 
# aligned along the z-axis (c axis). Ion orientation doesn't matter when calculating the dipole tensors
NdGaO3.add_position(0.0, 0.0, 0.25, (0,0,1))
NdGaO3.add_position(0.5, 0.5, 0.25, (0,0,1))
NdGaO3.add_position(0.5, 0.5, 0.75, (0,0,1))
NdGaO3.add_position(0.0, 0.0, 0.75, (0,0,1))
# ^ see Fig.3 of Luis et al - Magnetic susceptibility of NGO

try:
    f = gzip.open('NdGaO3_Jmatricies.dat.gz', 'rb')
    print("Loading matricies")
    J00mat, J01mat, J02mat, J03mat, J10mat, J11mat, J12mat, J13mat, J20mat, J21mat, J22mat, J23mat, \
    J30mat, J31mat, J32mat, J33mat = pickle.load(f)
    
except FileNotFoundError:
    print("Making matricies")

    J00 = NdGaO3.J_terms(200,0,0)
    J01 = NdGaO3.J_terms(200,0,1)
    J02 = NdGaO3.J_terms(200,0,2)
    J03 = NdGaO3.J_terms(200,0,3)
    J10 = NdGaO3.J_terms(200,1,0)
    J11 = NdGaO3.J_terms(200,1,1)
    J12 = NdGaO3.J_terms(200,1,2)
    J13 = NdGaO3.J_terms(200,1,3)
    J20 = NdGaO3.J_terms(200,2,0)
    J21 = NdGaO3.J_terms(200,2,1)
    J22 = NdGaO3.J_terms(200,2,2)
    J23 = NdGaO3.J_terms(200,2,3)
    J30 = NdGaO3.J_terms(200,3,0)
    J31 = NdGaO3.J_terms(200,3,1)
    J32 = NdGaO3.J_terms(200,3,2)
    J33 = NdGaO3.J_terms(200,3,3)
    
    
    J00mat = tup_to_array(J00)
    J01mat = tup_to_array(J01)
    J02mat = tup_to_array(J02)
    J03mat = tup_to_array(J03)
    J10mat = tup_to_array(J10)
    J11mat = tup_to_array(J11)
    J12mat = tup_to_array(J12)
    J13mat = tup_to_array(J13)
    J20mat = tup_to_array(J20)
    J21mat = tup_to_array(J21)
    J22mat = tup_to_array(J22)
    J23mat = tup_to_array(J23)
    J30mat = tup_to_array(J30)
    J31mat = tup_to_array(J31)
    J32mat = tup_to_array(J32)
    J33mat = tup_to_array(J33)
    
    NGOJmats = [ J00mat, J01mat, J02mat, J03mat, J10mat, J11mat, J12mat, J13mat, J20mat, J21mat, J22mat, J23mat, J30mat, J31mat, J32mat, J33mat]
    pickle.dump(NGOJmats, gzip.open('NdGaO3_Jmatricies.dat.gz', 'wb'))
    
## Jonos test stuff

# Calculate E1 
 
Sx = np.array([1, 0, 0])
Sy = np.array([0, 1, 0])
Sz = np.array([0, 0, 1])

S1 = 0.5*Sz
S2 = -0.5*Sz
S3 = 0.5*Sz
S4 = -0.5*Sx

V = NdGaO3.a[0]*NdGaO3.b[1]*NdGaO3.c[2]*(10**(-10))**3 
E1 = (-mu0*mB**2)/(8*pi*V)*S1.dot(J00mat.dot(S1) + J01mat.dot(S2) + J02mat.dot(S3) + J03mat.dot(S4))

print('Ground state energy (K) = ')
print(E1/k)