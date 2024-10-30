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
alpha = 0.3 # rad
beta = 0.5 # rad
wmag = 2*pi*1 # radHz
w = np.array([wmag*sin(beta),0,wmag*cos(beta)]) # vector for crossing
tmax = 5 # s
dt = 1e-5 # s
i = 0 # indicates which iteration of the ode solver we are on
# now we have calculated a starting chart
phi0,wt0 = np.meshgrid(np.arange(0,2*pi,2*pi/100), np.arange(0,tmax*wmag, wmag*dt))
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi0) - np.cos(phi0)*np.sin(alpha)*np.sin(wt0) + np.cos(beta)*np.sin(alpha)*np.cos(wt0)*np.sin(phi0))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt0))
v0 = 1e7*np.sign(lam0)
phi = np.full_like (phi0, np.nan)
lam = np.full_like(lam0, np.nan)
v = np.full_like(v0,np.nan)
phi[i,:] = phi0[i,:]
lam[i,:] = lam0[i,:]
v[i,:] = v0[i,:]
t = 0 # s
while t <= tmax:
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
	phi[i,:] = phi0[i,:]
	lam[i,:] = lam0[i,:]
	
	v[i,:] = v0[i,:]
	print(t)