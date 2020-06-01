#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', 'DipoleLatticeField.ipynb')


# In[2]:


ErLiF4 = Tetragonal()

#ErLiF4 parameters

#crystal axes
a_len = 5.162
c_len = 10.700
ErLiF4.axes(a_len,c_len)

V = a_len**2*c_len*(10**(-10))**3   #unit cell volume (SI units)

#g-tensor
gperp = 3.147   #g-tensor perpendicular principle value 
gpara = 8.105   #g-tensor parallel principle value
ErLiF4.g_tensor(gpara,gperp)

#Er ion positions
ErLiF4.add_position(0.5, 0.0, 0.00, (0,0,1))
ErLiF4.add_position(0.5, 0.5, 0.25, (0,0,1))
ErLiF4.add_position(0.0, 0.5, 0.50, (0,0,1))
ErLiF4.add_position(0.0, 0.0, 0.75, (0,0,1))


# In[3]:


R = 50   #radius of dipole to sum over (Angstroms)

def tup_to_array(tup):
    mat = np.array([[tup[0], tup[1], tup[2]],[tup[3], tup[4], tup[5]],[tup[6], tup[7], tup[8]]])
    return mat

#calculating J terms and converting from tuple into numpy arrays
J_dict = {}
J_dict[(0,0)] = tup_to_array(ErLiF4.J_terms(R,0,0))
J_dict[(0,1)] = tup_to_array(ErLiF4.J_terms(R,0,1))
J_dict[(0,2)] = tup_to_array(ErLiF4.J_terms(R,0,2))
J_dict[(0,3)] = tup_to_array(ErLiF4.J_terms(R,0,3))
J_dict[(1,0)] = tup_to_array(ErLiF4.J_terms(R,1,0))
J_dict[(1,1)] = tup_to_array(ErLiF4.J_terms(R,1,1))
J_dict[(1,2)] = tup_to_array(ErLiF4.J_terms(R,1,2))
J_dict[(1,3)] = tup_to_array(ErLiF4.J_terms(R,1,3))
J_dict[(2,0)] = tup_to_array(ErLiF4.J_terms(R,2,0))
J_dict[(2,1)] = tup_to_array(ErLiF4.J_terms(R,2,1))
J_dict[(2,2)] = tup_to_array(ErLiF4.J_terms(R,2,2))
J_dict[(2,3)] = tup_to_array(ErLiF4.J_terms(R,2,3))
J_dict[(3,0)] = tup_to_array(ErLiF4.J_terms(R,3,0))
J_dict[(3,1)] = tup_to_array(ErLiF4.J_terms(R,3,1))
J_dict[(3,2)] = tup_to_array(ErLiF4.J_terms(R,3,2))
J_dict[(3,3)] = tup_to_array(ErLiF4.J_terms(R,3,3))


# In[4]:


C = (-mu0*muB**2)/(8*pi*V)

#unit vectors
xhat = np.array([1, 0, 0])
yhat = np.array([0, 1, 0])
zhat = np.array([0, 0, 1])

def sin(angle):
    """Outputs the sine of the angle (in degrees).
    """
    return math.sin(math.radians(angle))

def cos(angle):
    """Outputs the cosine of the angle (in degrees).
    """
    return math.cos(math.radians(angle))

def Spin(phi,theta):
    """Outputs a spin vector
       Args: phi = azimuthal angle of the spin
             theta = polar angle of the spins
    """
    return 1/2*(cos(phi)*sin(theta)*xhat + sin(phi)*sin(theta)*yhat + cos(theta)*zhat)

def Single_spin_energy(sublattice,angles):
    """Outputs the energy of a single spin in the given sublattice as a function of spin configuration.
       Args: angles = list with elements phi1,theta1,phi2,theta2,phi3,theta3,phi4,theta4
             sublattice = number of the sublattice (0,1,2,3) of the spin for which the energy is calculated.
    """
    s1 = Spin(angles[0],angles[1])
    s2 = Spin(angles[2],angles[3])
    s3 = Spin(angles[4],angles[5])
    s4 = Spin(angles[6],angles[7])
    s_dict = {0:s1,1:s2,2:s3,3:s4}
    energy = 0
    
    for i in range(0,4):
        energy += C*s_dict[sublattice] @ J_dict[(sublattice,i)] @ s_dict[i]
    return energy

def Total_energy(angles):
    """Outputs the total energy of the given spin configuration
       Args: angles = list with elements phi1,theta1,phi2,theta2,phi3,theta3,phi4,theta4
    """
    energy = 0
    
    for i in range(0,4):
        energy += Single_spin_energy(i,angles)
    return energy


# In[5]:


E200 = -377.1434865925225
E100 = -377.04555364288296
E50  = -376.69933448489746
Eg = Single_spin_energy(0,[0,90,180,90,180,90,0,90])
print('E100 = ' + str(1e3*Eg/k) + ' mK')


# In[18]:


from scipy import optimize

angles0 = np.array([0,90,180,90,180,90,0,90])
result = optimize.minimize(Total_energy,angles0)
result.x

