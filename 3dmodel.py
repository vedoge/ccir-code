import numpy as np
rmin = 1e6 # cm
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.784e33 # g
# define coordinate system (lambda, phi) in region [-pi/2,pi/2]x[0,2pi]
# calculate starting lambda from alpha and phi
# cross product easy to implement in cartesian
# convert between the coordinate systems
def surf2cart(a):
	lam = a[0]
	phi = a[1]
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
# following parameters are of significance in the calculation
T = 1 # s
wmag = 2*np.pi/T # radHz
alpha = 0 # rad
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 #cm
# the following represents the starting point of the particle in question
phi = 0
# calculate initial elevation
lam = np.pi/2 - np.atan(1/(np.sin(phi)*np.tan(alpha)))
chi = np.tan(0.5/np.atan(lam))
# integrate the equation of motion
while rm*np.cos(lam)**2 >= rpuls:

