# finalmontecarlo.py
# Main Monte Carlo Simulation code
# full calculation 
# use with montecarloroutines.py
# add variable dt
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
sin = np.sin
cos = np.cos
tan = np.tan
asin = np.asin
acos = np.acos
atan = np.atan
sinh = np.sinh
cosh = np.cosh
tanh = np.tanh
asinh = np.asinh
acosh = np.acosh
atanh = np.atanh
dot = np.linalg.vecdot
cross = np.cross
norm = np.linalg.norm
sqrt = np.sqrt
# for threading (parallelism)
from threading import Thread, Lock
from concurrent.futures import ThreadPoolExecutor,wait 
# set up a stream of random numbers that are fast
rng = np.random.Generator(np.random.SFC64())
rand = rng.random
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.784e33 # g
rpuls = 1e6 # cm
rm = 1e8 # cm
maxlam = acos(sqrt(rpuls/rm))
alpha = 0.3 # rad
beta = 0.5 # rad
wmag = 2*pi*1 # radHz
# computing parameters
nr_compute_threads = 8 # customize to your needs
w = np.array([wmag*sin(-beta),0,wmag*cos(-beta)])
# define the array
def s(lam):
	return (rm/6)*(sqrt(3)*asinh(sqrt(3)*sin(lam)) + 3*sqrt(3)*sin(lam)*sqrt((sin(lam))**2 + 1/3))

def lat0(phi,wt): # FIXME
	return -asin(sin(phi)*(cos(alpha)*sin(beta) + cos(beta)*sin(alpha)*cos(wt)) - cos(phi)*sin(alpha)*sin(wt))
stot = s(maxlam)-s(-maxlam)

def update_region(plam,pphi,pv):
	global dt, lamendpoints,phiendpoints
	lamcells = np.searchsorted(lammidpoints,plam)
	phicells = np.searchsorted(phimidpoints,pphi)
	accel = a[phicells, lamcells]
	ds = pv*dt + 0.5*accel*(dt**2)
	plam += ds/(rm*cos(plam)*sqrt(1+3*(sin(plam)**2)))
	return plam 

def acc(lam,phi):
	# region - variable from 0 to 3
	# determines which phi coordinates are changed
		# region - variable from 0 to 3
	# determines which lambda coordinates are altered
	r = np.array([
		(rm*cos(lam)**2)*cos(lam)*cos(phi),
		(rm*cos(lam)**2)*cos(lam)*sin(phi),
		(rm*cos(lam)**2)*sin(lam)
	]) # mathematics convention (latitude)
	n = - np.array([
		sin(lam)*cos(lam)*cos(phi),
		sin(lam)*cos(lam)*sin(phi),
		2*sin(lam)**2-cos(lam)**2
	])
	n = n/norm(n,axis=0) # normalise
	# ag = -G*M*dot(r,n,axis=0)/norm(r,axis=0)**3
	# acen = dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
	return -G*M*dot(r,n,axis=0)/(norm(r,axis=0)**3) + dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)

lam = np.zeros(int(1e5))
lam[0] = -maxlam
for i in range(int(1e5)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(1e5*lamprime) # approx. 4*10^2 cm / cell; very reasonable computationally

# generate lambdas
phi_width = pi/250
lam,phi = np.meshgrid(lam,np.arange(0,2*pi,phi_width)) # granularity of phi reasonable
lamendpoints = (lam[0,1:] + lam[0,:-1])/2 #releases GIL - Main thread is allocated a quantum - we're cooked
phiendpoints = (phi[1:,0] + phi[:-1,0])/2

# lam[0,:] is increasing lambda
# phi[:,0] is increasing phi
a = np.empty_like(phi)

with ThreadPoolExecutor(max_workers = nr_compute_threads) as executor:
	a = np.array([res for res in executor.map(acc,lam,phi)]).reshape(phi.shape)
	t = 0 # s
	dt = 0 # s (breaks if doesn't work)
	j = 0
	T = None
	dMdt = None
	tfin = 2 # s
	partspercell = 10
	cells = phi[:,0].size
	parts = cells*partspercell
	dMdt_disk = 1e15 # g s^-1
	pphi = rand(size=(cells,partspercell))*phi_width + phi[:,0:(partspercell)]
	plam = lat0(pphi,0)
	pv = np.sign(plam)*1e7
	lamcells = np.searchsorted(lamendpoints,plam)
	phicells = np.searchsorted(phiendpoints,pphi)
	dt = np.min([np.abs(0.2*np.min((s(plam)-s(lam[phicells,lamcells]))/pv)),1e-5])
	mparts = np.full_like(pphi,dMdt_disk*dt/(parts)) # dM = dMdT * deltaT / nr. of particles 
	while t <= tfin:
		T = np.append(T,t)
		print(dt)
		plam = np.array([i for i in executor.map(update_region,plam,pphi,pv)])
		# find particles that land on surface
		pphi[np.where((plam >= maxlam) | (plam <= -maxlam))] = np.nan
		plam[np.where(pphi == np.nan)] = np.nan
		# process them for mass accretion rate
		dMdt = np.append(dMdt,mparts[np.where(plam==np.nan)]/dt) # calculate mass accretion rate
		# calculate dt
		lamcells = np.searchsorted(lamendpoints,plam)
		phicells = np.searchsorted(phiendpoints,pphi)
		t += dt
		j += 1
		# calculate dt for next iteration
		dt = np.min(0.2*np.min((s(plam)-s(lam[phicells,lamcells]))/pv),1e-5)
		# remove particles that have finished
		pphi = np.extract(pphi,np.where(np.isfinite(plam)))
		plam = np.extract(plam,np.where(np.isfinite(plam)))
		# add particles to simulation
		pphi_delta =  rand(size=(cells,partspercell))*phi_width + phi[:,0:(partspercell)]
		plam = np.append(plam, lat0(pphi_delta, wmag*t))
		pphi = np.append(pphi, pphi_delta)
		mparts = np.append(mparts,np.full_like(pphi_delta,dMdt_disk*dt/parts))
	plt.plot(T,dMdt)
	plt.savefig('fig.eps')
	plt.show()
