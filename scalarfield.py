import numpy as np
import matplotlib.pyplot as plt
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 10 # Hz
wmag = 2*np.pi*f # radHz
beta = 1 # rad
w = wmag*np.array([np.sin(beta),0,np.cos(beta)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 #cm
# make sure rpuls is reasonable because otherwise it messes up the contour map
M = 2.784e33 # g
# set up 2-dimensional space
maxlam = np.acos(np.sqrt(rpuls/rm))
lam = np.linspace(-maxlam,maxlam,1001) # this is extra slow with such a low number
phi = np.linspace(0,np.pi*2,1001) # this is actually okay
spacelam,spacephi = np.meshgrid(lam,phi) # create 3-D space for vectorisation 
r = np.array([
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.cos(spacephi),
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.sin(spacephi),
	(rm*np.cos(spacelam)**2)*np.sin(spacelam)
	]) # mathematics convention (latitude)
# centrifugal acceleration = -n dot (w*(w*r))*n
# gravitational acceleration = - n dot G*M*r/(np.linalg.norm(pos)**3)
# mathematics convention (latitude) vector along magnetic field lines for every point (lambda, phi)
# n = [-2*sin(lambda),cos(lambda),0] in [r,lambda,phi] coordinates
# n = [-2*sin(lambda),cos(lambda),0]*[r;theta;phi]
# if you put in other coordinates for r, theta, and phi, they become a matrix
# that yields the below equation for n in cartesian coordinates
n = - np.array([np.sin(spacelam)*np.cos(spacelam)*np.cos(spacephi),
				np.sin(spacelam)*np.cos(spacelam)*np.sin(spacephi),
				2*np.sin(spacelam)**2-np.cos(spacelam)**2
				])
n = n/np.linalg.norm(n,axis=0) # normalise
# now we need to evaluate acceleration at every point
ag = -G*M*np.linalg.vecdot(r,n,axis=0)/np.linalg.norm(r,axis=0)**3
# please find a way to make this agonisingly slow job faster
#for i in range(r.shape[1]):
#    for j in range(r.shape[2]):
#        acen[:,i,j] = -n[:,i,j]@np.cross(w,np.cross(w,r[:,i,j]))
#w*(w*r)
acen = np.linalg.vecdot(-n,np.cross(w,np.cross(w,r,axis=0),axis=0),axis=0) #problem here dotting n with each vector
ax = plt.axes(projection='3d')
ax.plot_surface(spacelam,spacephi,np.sign(ag+acen)*np.log10(np.abs(ag+acen)))
ax.set_xlabel(r"$\lambda$ / rad")
ax.set_ylabel(r"$\varphi$ / rad")
ax.set_zlabel(r"$\mathrm{sgn} \left( a \right) \log_{10} \left| a \right|$")
#plt.colorbar()
plt.show()
#lambda* = arctan(1/(sin(phi)*tan(beta)))
