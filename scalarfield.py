import numpy as np
import matplotlib.pyplot as plt
# define coordinate system (lambda, phi) in region [-pi/2,pi/2]x[0,2pi]
# calculate starting lambda from alpha and phi
# cross product easy to implement in cartesian
# convert between the coordinate systems
def surf2cart(lam,phi,rm):
	r = rm*np.cos(lam)**2
	coords = np.array([r*np.cos(lam)*np.cos(phi),r*np.cos(lam)*np.sin(phi),r*np.sin(lam)])
	return coords
def cart2surf(a):
	if a[0] != 0:
		phi = np.atan(a[1]/a[0])
	else:
		phi = 0
	lam = np.acos(np.sqrt(np.linalg.norm(a)/rm))
	return np.array([lam,phi])
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 10 # Hz
wmag = 2*np.pi*f # radHz
alpha = 0 # rad
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 #cm
M = 2.784e33 # g
# set up 2-dimensional space
lam = np.linspace(-np.pi/2,np.pi/2,1001)
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
print(ag.shape)
#plt.contourf(spacelam,spacephi,np.linalg.norm(ag,axis=0)) # is actually a very reasonable graph - the value explodes close to lambda=np.pi/2 (visible in the 10^50 cm^2 yellow bar on either end)
#plt.colorbar()
#plt.show()
# work out how to do a cross product of r[:,i,j] with w
# find a function lambda(phi,w,alpha,beta)
for i in range(r.shape[1]):
	for j in range(r.shape[2]):
		# loop over all the vectors to do |-n@(w x (w x r))| for each point
		# then add |ag| to the result
		#produce a contour chart of log10(result) against x
		# later, we write a simulation that operates on fixed arclength & variable dt
