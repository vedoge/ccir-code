import numpy as np
import matplotlib.pyplot as plt
alpha = 0.3 #rad
#10 KeV
beta =  0.5 #rad
#chi,wt = np.meshgrid(np.arange(0,2*np.pi,np.pi/1000),np.arange(0,2*np.pi,np.pi/100))
phi = np.arange(0,2*np.pi,np.pi/100)
#phi = np.acos((np.cos(chi)*np.cos(wt) - np.cos(alpha)*np.sin(chi)*np.sin(wt))/((np.cos(chi)*np.cos(wt) - np.cos(alpha)*np.sin(chi)*np.sin(wt))**2 + (np.sin(chi)*(np.sin(alpha)*np.sin(beta) - np.cos(alpha)*np.cos(beta)*np.cos(wt)) - np.cos(beta)*np.cos(chi)*np.sin(wt))**2)**(1/2))
#lam0 = np.acos(((np.cos(chi)*np.cos(wt) - np.cos(alpha)*np.sin(chi)*np.sin(wt))**2 + (np.sin(chi)*(np.sin(alpha)*np.sin(beta) - np.cos(alpha)*np.cos(beta)*np.cos(wt)) - np.cos(beta)*np.cos(chi)*np.sin(wt))**2)**(1/2))
wt = 0
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
plt.plot(phi,lam0)
wt = np.pi/2
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
plt.plot(phi,lam0)
wt = np.pi
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
plt.plot(phi,lam0)
wt = 3*np.pi/2
lam0 = np.asin(np.cos(alpha)*np.sin(beta)*np.sin(phi) - np.cos(phi)*np.sin(alpha)*np.sin(wt) + np.cos(beta)*np.sin(alpha)*np.cos(wt)*np.sin(phi))/(np.cos(alpha)*np.cos(beta) - np.sin(alpha)*np.sin(beta)*np.cos(wt))
plt.plot(phi,lam0)
plt.xlabel(r'$\varphi$ / rad')
plt.ylabel(r'$\lambda_0$ / rad')
plt.legend({r'$\omega t = 0$',r'$\omega t = \frac{\pi}{2}$',r'$\omega t = \pi$',r'$\omega t = \frac{3\pi}{2}$'})
plt.show()
