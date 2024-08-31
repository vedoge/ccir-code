import numpy as np
import math
def sph2cart(a):
    return ([a[0]*math.sin(a[1])*math.cos(a[2]),a[0]*math.sin(a[1])*math.sin(a[2]),a[0]*math.cos(a[1])])
v0 = np.array([0,0,1e5])
r0 = np.array([1e8,0,0])
n = np.array([0,0,1])
r = r0
v = v0
rmin = 1e6
G = 6.67e-8
M = 2.784e33
f = 0
lam = 0
t = 0
chi = math.pi/2
omega = np.array([2*math.pi*f,0,0])
dt = 1e-5
phi = 0
while np.linalg.norm(r) >= rmin:
    a = np.dot(-G*M*r/(np.linalg.norm(r)**3) + np.dot(n,np.cross(omega,np.cross(r,omega))),n)*n # projection in dir of n    
    r = r + v*dt + 0.5*(dt**2)*a # integrator
    v = v + a*dt
    lam = np.linalg.norm(v*dt + 0.5*(dt**2)*a)/(np.linalg.norm(r0)*math.cos(chi)*math.sqrt(1+3*math.sin(lam)**2))
    chi = math.atan(0.5/math.tan(lam))
    n = (np.array([math.sin(lam)*math.cos(phi),math.sin(lam)*math.sin(phi),math.cos(lam)]) - r/np.linalg.norm(r))
    n = n / np.linalg.norm(n)
    t += dt
    print(t,n) 
