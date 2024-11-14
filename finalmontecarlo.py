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
w = np.array([wmag*sin(alpha),0,wmag*cos(alpha)])
# define the array
def s(lam):
	return (rm/6)*(sqrt(3)*asinh(sqrt(3)*sin(lam)) + 3*sqrt(3)*sin(lam)*sqrt((sin(lam))**2 + 1/3))
stot = s(maxlam)-s(0)
print(stot)
# s = np.linspace(0,stot,1e6) #~100 cm per cell 
lam = np.zeros(int(1e6))
for i in range(int(1e6)-1):
	l = lam[i]
	lamprime = rm*cos(l)*sqrt(1+3*(sin(l)**2))
	lam[i+1] = l + stot/(1e6*lamprime) # approx. 10^2 cm / cell; very reasonable computationally
#plt.plot(s(lam))
#plt.show()
# generate lambdas
phi, lam = np.meshgrid(np.arange(0,2*pi,pi/500),lam)
# now make things
r = np.array([
	(rm*np.cos(lam)**2)*np.cos(lam)*np.cos(phi),
	(rm*np.cos(lam)**2)*np.cos(lam)*np.sin(phi),
	(rm*np.cos(lam)**2)*np.sin(lam)
	]) # mathematics convention (latitude)
# centrifugal acceleration = -n dot (w*(w*r))*n
# gravitational acceleration = - n dot G*M*r/(np.linalg.norm(pos)**3)
# mathematics convention (latitude) vector along magnetic field lines for every point (lambda, phi)
# n = [-2*sin(lambda),cos(lambda),0] in [r,lambda,phi] coordinates
# n = [-2*sin(lambda),cos(lambda),0]*[r;theta;phi]
# if you put in other coordinates for r, theta, and phi, they become a matrix
# that yields the below equation for n in cartesian coordinates
n = - np.array([np.sin(lam)*np.cos(lam)*np.cos(phi),
				np.sin(lam)*np.cos(lam)*np.sin(phi),
				2*np.sin(lam)**2-np.cos(lam)**2
				])
n = n/np.linalg.norm(n,axis=0) # normalise




