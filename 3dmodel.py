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

G = 6.67e-8 # dyn cm^2 g^-2
# following parameters are of significance in the calculation
T = 1 # s
wmag = 2*np.pi/T # radHz
alpha = 0 # rad
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 #cm
M = 2.784e33 # g
# the following represents the initial conditions of the particle in question
phi = 0
v0 = 1e5
dx = 0 # arclength (will be converted to delta-lambda)
# fixed timestep 
dt = 1e-5
# calculate initial elevation
lam = np.pi/2 - np.atan(1/(np.sin(phi)*np.tan(alpha)))
chi = np.tan(0.5/np.atan(lam))
# integrate the equation of motion
while rm*np.cos(lam)**2 >= rpuls:
	ag = G*m_pulsar * math.cos(chi) / (r**2) 
	# calculate centrifugal using -n dot omega cross (omega cross position)
	coord = surf2cart(lam,phi,rm) # r
	# NEED NORM FIXME
	norm = # n
	dx = v*dt + 0.5*a*(dt**2) # verlet
	v = v + a*dt # update velocity
	lam += dx / (rm*math.cos(lam) * (1+3*(math.sin(lam)**2))**0.5) # calculate delta-lambda based on eq. given
	chi = math.atan(1/(2*math.tan(lam)))
	r = rm * (math.cos(lam)**2)
	t += dt
