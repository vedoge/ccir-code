import numpy as np
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
rm = np.array(1e8) # cm
rpuls = np.array(1e6) #cm
M = 2.784e33 # g
# set up 2-dimensional space
lam = np.linspace(-np.pi/2,np.pi/2,1000)
phi = np.linspace(0,np.pi*2,1000)


ntheta = np.pi + lam - np.atan(0.5/np.tan(lam))
nphi = phi + np.pi
n = np.array(	[np.cos(ntheta)*np.cos(nphi),
				np.cos(ntheta)*np.sin(nphi),
				np.sin(ntheta)]
			)
# centrifugal acceleration = -n@(w*(w*r))*n
# gravitational acceleration = n*G*M*pos@n/(np.linalg.norm(pos)**3)
print(ntheta)
print(n)
