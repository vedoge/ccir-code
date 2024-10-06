import numpy as np
import matplotlib.pyplot as plt
alpha = 1 #rad
#10 KeV
beta =  0.5 #rad
chi,wt = np.meshgrid(np.arange(0,2*np.pi,np.pi/100),np.arange(0,2*np.pi,np.pi/100))
lam0 = np.acos((np.abs(np.cos(chi)*np.cos(wt) - np.cos(alpha)*np.sin(chi)*np.sin(wt))**2 + np.abs(np.sin(chi)*(np.sin(alpha)*np.sin(beta) - np.cos(alpha)*np.cos(beta)*np.cos(wt)) - np.cos(beta)*np.cos(chi)*np.sin(wt))**2)**(1/2))
plt.plot(chi[0,:],lam0[0,:])
plt.legend([r'$\omega t=0$',r'$\omega t=pi/2$',r'$\omega t=pi$',r'$\omega t=3pi/2$'])
plt.xlabel(r'$\chi$ / rad')
plt.ylabel(r'$\lambda_0$ / rad')
plt.show()

