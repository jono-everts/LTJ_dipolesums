import numpy as np
import DipoleLatticeField

def tup_to_array(tup):
    mat = np.array([[tup[0], tup[1], tup[2]],[tup[3], tup[4], tup[5]],[tup[6], tup[7], tup[8]]])
    return mat

#Define SI units
pi = np.pi
mB=9.274*10**(-24)
k=1.380*10**(-23)
NA=6.022*10**23
mu0=4*pi*10**(-7)


ErLiF4 = DipoleLatticeField.Tetragonal()
ErLiF4.axes(5.162, 10.70)
ErLiF4.g_tensor(8.105, 3.147)
#Define the ions to be in the bilayered antiferromagnetic ground state, aligned along the x-axis
#Ion orientation doesn't matter when calculating the dipole tensors
ErLiF4.add_position(0.5, 0.0, 0.00, (0,0,1))
ErLiF4.add_position(0.5, 0.5, 0.25, (0,0,1))
ErLiF4.add_position(0.0, 0.5, 0.50, (0,0,1))
ErLiF4.add_position(0.0, 0.0, 0.75, (0,0,1))

try:
    f = gzip.open('ErLiF4_Jmatricies.dat.gz', 'rb')
    print("Loading matricies")
    J00mat, J01mat, J02mat, J03mat, J10mat, J11mat, J12mat, J13mat, J20mat, J21mat, J22mat, J23mat, \
    J30mat, J31mat, J32mat, J33mat = pickle.load(f)
    
except FileNotFoundError:
    print("Making matricies")

    J00 = ErLiF4.J_terms(200,0,0)
    J01 = ErLiF4.J_terms(200,0,1)
    J02 = ErLiF4.J_terms(200,0,2)
    J03 = ErLiF4.J_terms(200,0,3)
    J10 = ErLiF4.J_terms(200,1,0)
    J11 = ErLiF4.J_terms(200,1,1)
    J12 = ErLiF4.J_terms(200,1,2)
    J13 = ErLiF4.J_terms(200,1,3)
    J20 = ErLiF4.J_terms(200,2,0)
    J21 = ErLiF4.J_terms(200,2,1)
    J22 = ErLiF4.J_terms(200,2,2)
    J23 = ErLiF4.J_terms(200,2,3)
    J30 = ErLiF4.J_terms(200,3,0)
    J31 = ErLiF4.J_terms(200,3,1)
    J32 = ErLiF4.J_terms(200,3,2)
    J33 = ErLiF4.J_terms(200,3,3)
    
    
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
    
    ELF_Jmats = [ J00mat, J01mat, J02mat, J03mat, J10mat, J11mat, J12mat, J13mat, J20mat, J21mat, J22mat, J23mat, J30mat, J31mat, J32mat, J33mat]
    pickle.dump(ELF_Jmats, gzip.open('ErLiF4_Jmatricies.dat.gz', 'wb'))

## Jonos test stuff

# Calculate E1 
 
xhat = np.array([1, 0, 0])
yhat = np.array([0, 1, 0])
zhat = np.array([0, 0, 1])

S1 = 0.5*xhat
S2 = -0.5*xhat
S3 = -0.5*xhat
S4 = 0.5*xhat

V = ErLiF4.a[0]*ErLiF4.b[1]*ErLiF4.c[2]*(10**(-10))**3 
E1 = (-mu0*mB**2)/(8*pi*V)*S1.dot(J00mat.dot(S1) + J01mat.dot(S2) + J02mat.dot(S3) + J03mat.dot(S4))
print('E1 = ' + str(E1/k) + ' K')
