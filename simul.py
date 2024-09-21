import numpy as np
import matplotlib.pyplot as plt
# universal gravitational constant
G = 6.67e-8 # dyn cm^2 g^-2
# parameters of significance in the calculation
f = 10 # Hz
wmag = 2*np.pi*f # radHz
alpha = 1 # rad; inclination of axis of rotation wrt magnetic moment
beta = 0.5 # rad; inclination of normal to disk wrt magnetic moment
w = wmag*np.array([np.sin(alpha),0,np.cos(alpha)]) # radHz (Cartesian)
rm = 1e8 # cm
rpuls = 1e6 # cm
# make sure rpuls is reasonable because otherwise it messes up the contour map
M = 2.784e33 # g
# set up 2-dimensional space
maxlam = np.acos(np.sqrt(rpuls/rm))
lam = np.arange(-maxlam,maxlam,maxlam/1000) # -maxlam -999maxlam/1000 ... 0 ... maxlam
phi = np.arange(0,np.pi*2,np.pi/500) # 0 pi/1000 2pi/1000 ... 2pi
spacelam,spacephi = np.meshgrid(lam,phi)
r = np.array([
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.cos(spacephi),
	(rm*np.cos(spacelam)**2)*np.cos(spacelam)*np.sin(spacephi),
	(rm*np.cos(spacelam)**2)*np.sin(spacelam)
	]) # mathematics convention (latitude)
# centrifugal acceleration = -n dot (w*(w*r))*n
# gravitational acceleration = - n dot G*M*r/(np.linalg.norm(pos)**3)
# mathematics convention (latitude) vector along magnetic field lines for every point (lambda, phi)
# n = [-2*sin(lambda),cos(lambda),0] in [r,lambda,phi] coordinates
# n = [-2*sin(lambda),cos(lambda),0]*[r;theta;phi]
# if you put in other coordinates for r, theta, and phi, they become a matrix
# that yields the below equation for n in cartesian coordinates
n =-np.array([np.sin(spacelam)*np.cos(spacelam)*np.cos(spacephi),
			  np.sin(spacelam)*np.cos(spacelam)*np.sin(spacephi),
			  2*np.sin(spacelam)**2-np.cos(spacelam)**2
			  ]) # put through the same matrix transposition
n = n/np.linalg.norm(n,axis=0) # normalise
# now we need to evaluate acceleration at every point
ag = -G*M*np.linalg.vecdot(r,n,axis=0)/np.linalg.norm(r,axis=0)**3
# please find a way to make this agonisingly slow job faster
#for i in range(r.shape[1]):
#    for j in range(r.shape[2]):
#        acen[:,i,j] = -n[:,i,j]@np.cross(w,np.cross(w,r[:,i,j]))
#w*(w*r)
acen = np.linalg.vecdot(-n,np.cross(w,np.cross(w,r,axis=0),axis=0),axis=0)
a = ag + acen
plt.contourf(spacelam,spacephi,np.sign(ag+acen)*np.log10(np.abs(ag+acen)))
plt.xlabel("lambda / rad")
plt.ylabel("phi / rad")
plt.title("Contour plot of sgn(a)log10(|a|) / cm s^-2")
plt.colorbar()
lam0 = np.pi/2-np.atan(1/(np.sin(phi)*np.tan(beta))) # no problem
lam0[lam0>=np.pi/2] -= np.pi # form a continuous line
plt.plot(lam0,phi)
plt.show()
particle_phi = 0 # rad
phi_idx = np.abs(phi-particle_phi).argmin() # index 1
lam_idx = np.abs(lam-lam0[phi_idx]).argmin() # index 0
v0 = +10^5 # cm s^-1
t = 0
dt = 0
v = v0 # positive acceleration is in direction of increasing lambda
dlam = maxlam/1000
ds = rm*np.cos(lam)*np.sqrt(1+3*np.sin(lam)**2)*dlam
while lam_idx > 0 | lam_idx < 1999:
    # phi_idx is constant (no azimuthal motion)
    # use formula rm*cos(lambda)*sqrt(1+3sin^2(lambda))dlambda
    if np.sign(v) > 0:
    dt = (v - np.sqrt(v^2 - 2*a[lam_idx,phi_idx]*ds[lam_idx]))/(2*a[lam_idx,phi_idx])
    t += dt
    v += a[lam_idx,phi_idx]*dt
    
