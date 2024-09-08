import numpy as np
import matplotlib.pyplot as plt
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 10 # Hz
wmag = 2*np.pi*f # radHz
alpha = 1 # rad
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 #cm
# make sure rpuls is reasonable because otherwise it messes up the contour map
M = 2.784e33 # g
# set up 2-dimensional space
maxlam = np.acos(np.sqrt(rpuls/rm))
lam = np.linspace(-maxlam,maxlam,1001)
phi = np.linspace(0,np.pi*2,1001)
spacelam,spacephi = np.meshgrid(lam,phi)
r = np.array([
			(rm*np.cos(spacelam)**3)*np.cos(spacephi),
			(rm*np.cos(spacelam)**3)*np.sin(spacephi),
			(rm*np.cos(spacelam)**2)*np.sin(spacelam)
			]) # mathematics convention (latitude)
# centrifugal acceleration = -n@(w*(w*r))*n
# gravitational acceleration = n*G*M*pos@n/(np.linalg.norm(pos)**3)
slope_lat = np.atan(0.5/np.tan(spacelam))-spacelam # signed latitude for normal vector at each point
n = np.array([np.cos(slope_lat)*np.cos(spacephi), np.cos(slope_lat)*np.sin(spacephi), np.sin(slope_lat)]) # mathematics convention (latitude)
# now we need to evaluate acceleration at every point (thank god for SIMD)
ag = -G*M*r@n/(np.linalg.norm(r,axis=0)**3)
acen = np.zeros(r.shape) # create array full of 
#plt.contourf(spacelam,spacephi,np.linalg.norm(ag,axis=0)) # is actually a very reasonable graph - the value explodes close to lambda=np.pi/2 (visible in the 10^50 cm^2 yellow bar on either end)
#plt.colorbar()
#plt.show()
# work out how to do a cross product of r[:,i,j] with w
# find a function lambda(phi,w,alpha,beta)
# please find a way to make this agonisingly slow job faster
#for i in range(r.shape[1]):
#    for j in range(r.shape[2]):
#        acen[:,i,j] = -n[:,i,j]@np.cross(w,np.cross(w,r[:,i,j]))
acen = -n@np.cross(w,np.cross(w,r,axis=0),axis=0)
a = np.linalg.norm(acen, axis=0) + np.linalg.norm(ag,axis=0)
loga = np.log10(a)
plt.contourf(spacelam,spacephi,loga)
plt.xlabel("lambda / rad")
plt.ylabel("phi / rad")
plt.title("Contour plot of log(a) / log(cm s^-2)")
plt.colorbar()
plt.show()
