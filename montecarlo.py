import numpy as np
pi = np.pi
sin = np.sin
cos = np.cos
tan = np.tan
asin = np.asin                                                                                                                                                       
acos = np.acos
atan = np.atan
dot = np.linalg.vecdot
cross = np.cross
norm = np.linalg.norm
sqrt = np.sqrt
import matplotlib.pyplot as plt
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.784e33 # g
rpuls = 1e6 # cm
rm = 1e8 # cm
alpha 0.3 # rad
beta = 0.5 # rad
wmag = 2*pi*1 # radHz
w = array([wmag*sin(beta),0,wmag*cos(beta)]) # vector for crossing 
t = 0 # s
dt = 1e-4 # s
# one particle per phi
phi = np.arange(0,2*pi,pi/1000)
wt = wmag*t # ie 0
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
# choose 100 phi-cells for speed
# recalculate lam0 every dt - it's less memory intensive than precalculating
v0 = 1e7 # cm s^-1
v = np.sign(lam0)*v0
lam = lam0
stopped = np.full_like(phi, False) # keep a record of the particles that have stopped 
# eventually we will end up with 10^8 particles in the simulation (more than enough)
while t <= 5:
	r = np.array([
	(rm*cos(lam)**2)*cos(lam)*cos(phi),
	(rm*cos(lam)**2)*cos(lam)*sin(phi),
	(rm*cos(lam)**2)*sin(lam)
	])
	n =-np.array([sin(lam)*cos(lam)*cos(phi),
	sin(lam)*cos(lam)*sin(phi),
	2*sin(lam)**2-cos(lam)**2
	])
	n = n / norm(n,axis=0)
	ag = -G*M*dot(r,n,axis=0)/norm(r,axis=0)**3
	acen = dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
	a = ag + acen
	ds = v*dt + a*(dt**2)*0.5
	v += a*dt
	t += dt
	lam += ds/(rm*cos(lam)*sqrt(1 + 3*(sin(lam)**2)))
	# add some n number of particles
	wt = wmag*dt
	lam = np.append(lam,np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt)))
	# start the particles
	phi = np.append(phi,np.arange(0,2*pi,np.pi/1000)
# add 2e3 particles every 1e-4 s; 2e7 particles per second 





