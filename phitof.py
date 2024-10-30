import numpy as np
import matplotlib.pyplot as plt
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
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 1 # Hz
wmag = 2*pi*f # radHz
alpha = 0.3 # rad; inclination of axis of rotation wrt disk
beta = 0.5 # rad; inclination of magnetic moment wrt axis of rotation
w = wmag*np.array([sin(beta),0,cos(beta)]) # radHz (Cartesian, b-vec coords)
rm = 1e8 # cm
rpuls = 1e6 # cm
# make sure rpuls is reasonable because otherwise it messes up the contour map
M = 2.784e33 # g
# set up 2-dimensional space
maxlam = acos(np.sqrt(rpuls/rm)) # works
# n = [-2*sin(lambda),cos(lambda),0] in [r,lambda,phi] coordinates
# n = [-2*sin(lambda),cos(lambda),0]*[r;lambda;phi]
# if you put in other coordinates for r, theta, and phi, they become a matrix
# that yields the below equation for n in cartesian coordinates
dt = 1e-4 # s
v0 = 1e7 # cm s^-1
t = 0
wt,phi = np.meshgrid(np.arange(0,2*pi,wmag*dt),np.arange(0,2*pi,pi/1000))
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
# start a few particles at different phi cells at wt=0
t = 0
phi = phi[:,0]
lam = lam0[:,0]
print(lam)
v = np.full_like(phi,v0)
v[lam < 0] = -v0
stopped = np.full_like(phi,False)
# we need to associate particles with all phi at wt=0
tof = np.full_like(phi,0) # time of flight
while (not np.all(stopped)) and t <= 5:
	# two fates - impact the surface or be stopped by the potential barrier and return to origin
	# sieve out both at the end of the ode solver
	# need to pick a phi
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
	# record time of flight of particles
	idx = (np.abs(lam) >= maxlam) & (stopped == False)
	print(np.sum(stopped),np.all(stopped))
	tof[idx] = t # if the particle has just stopped &
	stopped[idx] = True
	# fate 2 - particle stopped by potential barrier and returned to starting point
	idx = ((v == 0) & (a == 0)) # acceleration directed to 0 lambda & particle has no KE left (return to centre) 
	stopped[idx] = True
print(tof)
plt.plot(phi,tof)
plt.show()
	