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
import concurrent.futures
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
def lat0(phi,wt):
	return np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
stot = s(maxlam)-s(-maxlam)
def acc(lam,phi,arr):
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
	return -G*M*dot(r,n,axis=0)/norm(r,axis=0)**3 + dot(-n,cross(w,cross(w,r,axis=0),axis=0),axis=0)
lam = np.zeros(int(1e5))
lam[0] = -maxlam
for i in range(int(1e5)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(1e5*lamprime) # approx. 2*10^3 cm / cell; very reasonable computationally

# generate lambdas
phi, lam = np.meshgrid(np.arange(0,2*pi,pi/50),lam)
# now make things
a = np.empty_like(phi)
ac = np.empty.like(phi)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
	a = 
# now, we have acceleration for each lambda and phi
t = 0
lam0 = lat0(phi[0,:],wmag*t)
# add N particles to each cell (add NM particles)
# M = 100, N = 2
