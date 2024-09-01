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
f = 1 # Hz
wmag = 2*np.pi*f # radHz
alpha = 0 # rad
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = np.array(1e8) # cm
rpuls = np.array(1e6) #cm
M = 2.784e33 # g
# initial conditions of the particle in question
phi = 0
v = 1e5
dx = 0 # arclength (will be converted to delta-lambda)
# fixed timestep 
t = 0
dt = 1e-4
# calculate initial elevation
if phi != 0 and alpha != 0:
	lam = np.pi/2 - np.atan(1/(np.sin(phi)*np.tan(alpha)))
	chi = np.tan(0.5/np.atan(lam))
else: 
	lam = 0
	chi = np.pi/2
r = rm*np.cos(lam)**2
# integrate the equation of motion
while np.any(r >= rpuls):
	# calculate centrifugal using -n dot omega cross (omega cross position)
	ntheta = 3*np.pi/2 + lam - chi
	nphi = phi + np.pi
	n = np.array([np.cos(ntheta)*np.cos(nphi),np.cos(ntheta)*np.sin(nphi),np.sin(ntheta)])
	pos = surf2cart(lam,phi,rm)
	acen = np.dot(n,np.cross(w,np.cross(w,pos)))*n
	print(t,v)
	ag = n*G*M* np.cos(chi) / (r**2)
	a = ag + acen
	dx = v*dt + 0.5*a*(dt**2) # verlet
	v = v + a*dt # update velocity
	lam += np.linalg.norm(dx) / (rm*np.cos(lam) * (1+3*(np.sin(lam)**2))**0.5) # calculate delta-lambda based on eq. given
	chi = np.atan(1/(2*np.tan(lam)))
	r = rm * (np.cos(lam)**2)
	t += dt
