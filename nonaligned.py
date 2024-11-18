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
alpha = 0.3 #rad
#10 KeV
beta =  0.5 #rad
phi = np.arange(0,2*np.pi,np.pi/100)
wt = 0
lam0 = -asin(sin(phi)*(cos(alpha)*sin(beta) + cos(beta)*sin(alpha)*cos(wt)) - cos(phi)*sin(alpha)*sin(wt))
plt.plot(phi,lam0)
wt = np.pi/2
lam0 = -asin(sin(phi)*(cos(alpha)*sin(beta) + cos(beta)*sin(alpha)*cos(wt)) - cos(phi)*sin(alpha)*sin(wt))
plt.plot(phi,lam0)
wt = np.pi
lam0 = -asin(sin(phi)*(cos(alpha)*sin(beta) + cos(beta)*sin(alpha)*cos(wt)) - cos(phi)*sin(alpha)*sin(wt))
plt.plot(phi,lam0)
wt = 3*np.pi/2
lam0 = -asin(sin(phi)*(cos(alpha)*sin(beta) + cos(beta)*sin(alpha)*cos(wt)) - cos(phi)*sin(alpha)*sin(wt))
plt.plot(phi,lam0)
plt.xlabel(r'$\varphi$ / rad')
plt.ylabel(r'$\lambda_0$ / rad')
plt.legend({r'$\omega t = 0$',r'$\omega t = \frac{\pi}{2}$',r'$\omega t = \pi$',r'$\omega t = \frac{3\pi}{2}$'})
plt.show()
