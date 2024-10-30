import numpy as np
import matplotlib.pyplot as plt
# the below simulation is very accurate (constant displacement sim) but it's agonizingly slow on modern hardware. Maybe in 2040 with an overclocked 1 nm processor...
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.784e33 # g
rm = 1e8 # cm
rpuls = 1e6 # cm
maxlam = np.acos(np.sqrt(rpuls/rm)) 
wmag = 2 * np.pi * 0.1 # radHz
alpha = 0.2 # rad
w = np.array([wmag*np.sin(alpha),0,wmag*np.cos(alpha)])
beta = 0.1 # rad
phi = np.arange(0,2*np.pi,np.pi/1000)
wt = 0
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))

ds = 1 # cm
# dlam = cos(lambda)*sqrt(1+3*sin^2(lambda))*ds
v0 = 1e7 # cm s^-1
v = v0
phi = 0 # rad
lam = 0 # rad
t = 0 # s 
while lam <= maxlam:
	r = np.array([
	(rm*np.cos(lam)**2)*np.cos(lam)*np.cos(phi),
	(rm*np.cos(lam)**2)*np.cos(lam)*np.sin(phi),
	(rm*np.cos(lam)**2)*np.sin(lam)
	])
	n =-np.array([np.sin(lam)*np.cos(lam)*np.cos(phi),
		np.sin(lam)*np.cos(lam)*np.sin(phi),
		2*np.sin(lam)**2-np.cos(lam)**2
		]) # put through the same matrix transposition
	n = n/np.linalg.norm(n,axis=0) # normalise
	a = np.linalg.vecdot(n,-G*M*r/np.linalg.norm(r)**3) + np.linalg.vecdot(-n,np.cross(w,np.cross(w,r)))
	dt = (np.sqrt(v**2 + 2*a*ds) - v)/a
	if np.isnan(dt):
		break
	v += a*dt
	t += dt
	lam += ds/(rm*np.cos(lam)*np.sqrt(1+3*np.sin(lam)**2))
	print (t,lam)
	