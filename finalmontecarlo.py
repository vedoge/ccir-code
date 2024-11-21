# finalmontecarlo.py
# Main Monte Carlo Simulation code
# full calculation 
# use with montecarloroutines.py
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
# from concurrent.futures import ThreadPoolExecutor,wait 
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
	return np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
stot = s(maxlam)-s(-maxlam)
def update_region(region, dt):
	global lam, phi,a
	global plam,pphi,pv
	idx = np.arange(int(region*plam.shape[0]/nr_compute_threads),int((region+1)*plam.shape[0]/nr_compute_threads))
	#lam_parts = plam[idx].reshape(1,1,plam[idx].size)
	#phi_parts = pphi[idx].reshape(1,1,pphi[idx].size)
	# there is a better alternative. Find out midpoints. Use them. 
	#lamcells = np.argmin(np.abs(lam_parts - lam[0,:,np.newaxis]),axis=1).squeeze()
	#phicells = np.argmin(np.abs(phi_parts - phi[:,0, np.newaxis]),axis=1).squeeze()
	lammidpoints = (lam[0,1:] + lam[0,:-1])/2 #releases GIL - Main thread is allocated a quantum - we're cooked
	# use mutex
	phimidpoints = (phi[1:,0] + phi[:-1,0])/2
	lamcells = np.searchsorted(lammidpoints,plam[idx])
	phicells = np.searchsorted(phimidpoints,pphi[idx])
	accel = a[phicells, lamcells]
	ds = pv[idx]*dt + 0.5*accel*(dt**2)
	plam[idx] += ds/(rm*cos(plam[idx])*sqrt(1+3*(sin(plam[idx])**2)))
def acc(lam,phi,region):
	# region - variable from 0 to 3
	# determines which phi coordinates are changed
	global a
	idx = np.arange(int(region*lam.shape[1]/nr_compute_threads),int((region+1)*lam.shape[1]/nr_compute_threads))
	lamref = lam[:,idx]
	phiref = phi[:,idx]
		# region - variable from 0 to 3
	# determines which lambda coordinates are altered
	r = np.array([
		(rm*cos(lamref)**2)*cos(lamref)*cos(phiref),
		(rm*cos(lamref)**2)*cos(lamref)*sin(phiref),
		(rm*cos(lamref)**2)*sin(lamref)
	]) # mathematics convention (latitude)
	n = - np.array([
		sin(lamref)*cos(lamref)*cos(phiref),
		sin(lamref)*cos(lamref)*sin(phiref),
		2*sin(lamref)**2-cos(lamref)**2
	])
	n = n/norm(n,axis=0) # normalise
	# ag = -G*M*dot(r,n,axis=0)/norm(r,axis=0)**3
	# acen = dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
	a[:,idx] = -G*M*dot(r,n,axis=0)/(norm(r,axis=0)**3) + dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
'''
running_acc_threads_mutex.acquire()
running_acc_threads -= 1
running_acc_threads.release()
'''
lam = np.zeros(int(1e5))
lam[0] = -maxlam
for i in range(int(1e5)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(1e5*lamprime) # approx. 2*10^3 cm / cell; very reasonable computationally
# generate lambdas
lam,phi = np.meshgrid(lam,np.arange(0,2*pi,pi/50)) # the relatively low granularity in phi (100 cells) is actually relatively reasonable - it's not as necessary as that of lambda
# lam[0,:] is increasing lambda
# phi[:,0] is increasing phi
# now make things
a = np.empty_like(phi)
threads = [Thread(target=acc, args=(lam,phi,i)).run() for i in range(nr_compute_threads)]
t = 0 # s
dt = 1e-4 # s (initial)
j = 0
tfin = 2*dt # s
dMdt = np.zeros(int(tfin / dt))
print(a)
partspercell = 10 # 500 particles per 10^-5 s (5e7 particles per second, or in other words, 2e10 g per particle) 
Mpart = 1e16/(phi.shape[0]*partspercell) # mass of individual particle
pphi = np.arange(0,2*pi,pi/250)
plam = lat0(pphi,0)
pv = np.sign(plam)*1e7
while t <= tfin:
	threads = [Thread(target=update_region, args=(i,dt)) for i in range(nr_compute_threads)]
	for i in threads: 
		i.start()
	for i in threads: 
		i.join()
	lost = (np.abs(plam) >= maxlam)
	#print(t,np.sum(lost))
	#print(np.max(np.abs(plam)))
	dMdt[j] = np.sum(lost)*Mpart
	if dMdt[j]:
		print(t,dMdt[j])
	# add new particles 
	plam = plam[~lost]
	pphi = pphi[~lost] 
	t += dt
	j += 1
	pphi = np.append(pphi, phi[:,0].copy())
	plam = np.append(plam, lat0(phi[:,0], wmag*t))
# add N particles to each cell (add NM particles)
# M = 100, N = 2
# you need to multithread the update function in the same way, and just finish the rest of the simulation
# you're very close.
plt.plot(np.arange(0,tfin,dt),dMdt)
plt.show()
