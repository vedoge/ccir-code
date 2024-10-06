import numpy as np
import matplotlib.pyplot as plt
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 10 # Hz
wmag = 2*np.pi*f # radHz
alpha = 1 # rad; inclination of axis of rotation wrt magnetic moment
beta = 0.5 # rad; inclination of normal to disk wrt magnetic moment
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 # cm
# make sure rpuls is reasonable because otherwise it messes up the contour map
M = 2.784e33 # g
alpha = 1 #rad
beta = 0.5 #rad
phi,wt = np.meshgrid(np.arange(0,2*np.pi,np.pi/1000),np.arange(0,2*np.pi,np.pi/1000))

zimg = np.array([np.sin(alpha)*np.sin(wt),
                 -(np.sin(beta)*np.cos(wt)*np.cos(alpha)+np.sin(alpha)*np.cos(beta)),
                 -np.sin(beta)*np.cos(wt)*np.sin(alpha)+np.cos(alpha)*np.cos(beta)
                 ])

lam0 = np.atan(np.linalg.norm(zimg,axis=0)*np.sin(phi)) # works

p# set up 2-dimensional space
maxlam = np.acos(np.sqrt(rpuls/rm)) # works 
np.arange(-maxlam,maxlam,maxlam/1000) # -maxlam -999maxlam/1000 ... 0 ... maxlam
spacelam,spacephi = np.meshgrid(np.arange(-maxlam,maxlam,maxlam/1000),np.arange(0,2*np,pi,np.pi/1000)) # also works
r = np.array([
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.cos(spacephi),
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.sin(spacephi),
	(rm*np.cos(spacelam)**2)*np.sin(spacelam)
	])
# n = [-2*sin(lambda),cos(lambda),0] in [r,lambda,phi] coordinates
# n = [-2*sin(lambda),cos(lambda),0]*[r;lambda;phi]
# if you put in other coordinates for r, theta, and phi, they become a matrix
# that yields the below equation for n in cartesian coordinates
n =-np.array([np.sin(spacelam)*np.cos(spacelam)*np.cos(spacephi),
			  np.sin(spacelam)*np.cos(spacelam)*np.sin(spacephi),
			  2*np.sin(spacelam)**2-np.cos(spacelam)**2
			  ]) # put through the same matrix transposition
n = n/np.linalg.norm(n,axis=0) # normalise
# now we need to evaluate acceleration at every point
ag = -G*M*np.linalg.vecdot(r,n,axis=0)/np.linalg.norm(r,axis=0)**3
# please find a way to make this agonisingly slow job faster
#for i in range(r.shape[1]):
#    for j in range(r.shape[2]):
#        acen[:,i,j] = -n[:,i,j]@np.cross(w,np.cross(w,r[:,i,j]))
acen = np.linalg.vecdot(-n,np.cross(w,np.cross(w,r,axis=0),axis=0),axis=0)
a = ag + acen
plt.contourf(spacelam,spacephi,np.sign(ag+acen)*np.log10(np.abs(ag+acen)))
plt.xlabel("lambda / rad")
plt.ylabel("phi / rad")
plt.title("Contour plot of sgn(a)log10(|a|) / cm s^-2")
plt.colorbar()
plt.plot(phi[0,:],lam0[0,:],)
plt.plot(phi[24,:],lam0[24,:])
plt.plot(phi[49,:],lam0[49,:])
plt.plot(phi[74,:],lam0[74,:])
plt.legend([r'\omega t=0',r'\omega t=pi/2',r'\omega t=pi',r'\omega t=3pi/2'])
plt.show()
particle_phi = 0 # rad
v0 = +10^5 # cm s^-1
t = 0
dt = 0
v = v0 # positive acceleration is in direction of increasing lambda
#lam0[phi,wt]
dlam = 
ds = rm*np.cos(spacelam)*np.sqrt(1+3*np.sin(spacelam)**2)*dlam
dt = 
