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
# set up a stream of random numbers
rng = np.random.Generator(np.random.SFC64())
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.784e33 # g
rpuls = 1e6 # cm
rm = 1e8 # cm
maxlam = acos(sqrt(rpuls/rm))
alpha = 0.3 # rad
beta = 0.5 # rad
wmag = 2*pi*1 # radHz
w = np.array([wmag*sin(beta),0,wmag*cos(beta)]) # vector for crossing
tmax = 2 # s
dt = 1e-5 # s
particles_per_dt = np.intp(1e3)
# I'm thinking - classify the particles into the cells that they're in
i = 1 # indicates which iteration of the ode solver we are on
# now we have calculated a starting chart
phi0,wt0 = np.meshgrid(np.empty(particles_per_dt), np.arange(0,tmax*wmag, wmag*dt))
phi0 = 2*pi*rng.random(phi0.size).reshape(phi0.shape)
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi0) - np.cos(phi0)*np.sin(alpha)*np.sin(wt0) + np.cos(beta)*np.sin(alpha)*np.cos(wt0)*np.sin(phi0))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt0))
v0 = 1e7*np.sign(lam0)
stopped0 = np.full_like(phi0, False, dtype='bool')
# currently coming up with a better, more performant alternative.
# phi = np.full_like (phi0, np.nan)
# lam = np.full_like(lam0, np.nan)
# v = np.full_like(v0,np.nan)
phi = phi0[0:i,:]
lam = lam0[0:i,:]
v = v0[0:i,:]
stopped = stopped0[0:i,:]
dMdt = np.zeros(np.intp((tmax+dt)/dt)) # collect this many results
t = 0 # s
while t <= tmax:
	lam[stopped == True] = np.nan
	phi[stopped == True] = np.nan
	v[stopped == True] = np.nan
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
	i += 1
	lam += ds/(rm*cos(lam)*sqrt(1 + 3*(sin(lam)**2)))
	dMdt[i-2] = np.sum((np.abs(lam) >= maxlam) & (stopped == False))*(1.672e-24)/(dt) # half protons half electrons - plasma with zero metallicity
	phi = phi0[0:i,:]
	lam = lam0[0:i,:]
	v = v0[0:i,:]
	stopped = stopped0[0:i,:]
	stopped[np.abs(lam) >= maxlam] = True
	print(t)
plt.plot(np.arange(0,t,dt),dMdt)
plt.xlabel(r"t / s")
plt.ylabel(r"$\dot{M}$ / $\mathrm{g s}^{-1}$")
plt.show()
