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
# for multiprocessing
from threading import Thread
# set up a stream of random numbers
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
w = np.array([wmag*sin(alpha),0,wmag*cos(alpha)])
# define the array
def s(lam):
	return (rm/6)*(sqrt(3)*asinh(sqrt(3)*sin(lam)) + 3*sqrt(3)*sin(lam)*sqrt((sin(lam))**2 + 1/3))
def lat0(phi,wt): # FIXME
	return np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
stot = s(maxlam)-s(-maxlam)
def update_region(region, dt):
	global lam, phi
	global a,plam,pphi,pv
	lam_parts = plam.reshape(1,1,plam.size)
	phi_parts = pphi.reshape(1,1,pphi.size)
	# use the 2nd axis
	# use an Event object to signal (releases GIL when thread is waiting)
	# minimize each axis individually
	print(phi.shape)
	# correct the below with axis= args

	lamcell = np.argmin(np.abs(lam_parts - lam[0,:,np.newaxis]))
	phicell = np.argmin(np.abs(phi_parts - phi[:,0, np.newaxis]))
	print(lamcell,phicell)
def acc(lam,phi):
	global a
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
	a = -G*M*dot(r,n,axis=0)/norm(r,axis=0)**3 + dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
lam = np.zeros(int(1e5))
lam[0] = -maxlam
for i in range(int(1e5)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(1e5*lamprime) # approx. 2*10^3 cm / cell; very reasonable computationally
# generate lambdas
phi, lam = np.meshgrid(np.arange(0,2*pi,pi/50),lam)
pphi = np.arange(0, 2*pi, pi/50)
plam = np.zeros_like(pphi)
# now make things
a = np.empty_like(phi)
t1 = Thread(target=acc,args=(phi,lam))
t1.run()
t2 = Thread(target = update_region, args=(1.5,1e-5)) 
t2.run()
t = 0
lam0 = lat0(phi[0,:],wmag*t)
# add N particles to each cell (add NM particles)
# M = 100, N = 2
