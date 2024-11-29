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
	global dt
	lammidpoints = (lam[0,1:] + lam[0,:-1])/2 #releases GIL - Main thread is allocated a quantum - we're cooked
	phimidpoints = (phi[1:,0] + phi[:-1,0])/2
	lamcells = np.searchsorted(lammidpoints,plam)
	phicells = np.searchsorted(phimidpoints,pphi)
    # address phi and lam with phicells and lamcells
    # to reduce the dimension of the calculation
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

lam = np.zeros(int(5e5))
lam[0] = -maxlam
for i in range(int(5e5)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(5e5*lamprime) # approx. 4*10^2 cm / cell; very reasonable computationally

# generate lambdas
phi_width = pi/250
lam,phi = np.meshgrid(lam,np.arange(0,2*pi,phi_width)) # granularity of phi reasonable
# lam[0,:] is increasing lambda
# phi[:,0] is increasing phi
a = np.empty_like(phi)

with ThreadPoolExecutor(max_workers = nr_compute_threads) as executor:
	a = np.array([res for res in executor.map(acc,lam,phi)]).reshape(phi.shape)

	t = 0 # s
	dt = 1e-4 # s (initial)
	j = 0
	tfin = 2 # s
	dMdt = np.zeros(int(tfin / dt)+1)
	partspercell = 10
	cells = phi[:,0].size
	dMdt_disk = 1e15 # g s^-1
	Mpart = dMdt*dt/(phi.shape[0]*partspercell) # mass of individual particle
	pphi = rand(size=(cells,partspercell))*phi_width + phi[:,0:(partspercell)]
	plam = lat0(pphi,0)
	pv = np.sign(plam)*1e7
	while t <= tfin:
		plam = np.array([i for i in executor.map(update_region,plam,pphi,pv)])
		print(np.sum(np.where((plam>=maxlam) | (plam <= -maxlam))))
		dMdt[j] = np.sum(np.where((plam >= maxlam) | (plam <= -maxlam)))*Mpart
		if dMdt[j] != dMdt[j-1]:
			print(dMdt[j])
		pphi[np.where((plam >= maxlam) | (plam <= -maxlam))] = np.nan
		plam[np.where(pphi == np.nan)] = np.nan
		#print(np.max(np.abs(plam)))
		# add new particles 
		t += dt
		j += 1
		pphi = np.extract(pphi,np.where(np.isfinite(plam)))
		plam = np.extract(plam,np.where(np.isfinite(plam)))
		pphi_delta =  rand(size=partspercell*cells)*phi_width
		pphi_delta += np.repeat(phi[:,0],partspercell).flatten()
		plam = np.append(plam, lat0(pphi_delta, wmag*t).copy())
	plt.plot(np.arange(0,tfin+dt,dt),dMdt)
	plt.savefig('fig.eps')
	plt.show()
