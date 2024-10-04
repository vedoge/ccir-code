import numpy as np
import matplotlib.pyplot as plt
alpha = np.pi/2 #rad
beta =  0 #rad
phi,wt = np.meshgrid(np.arange(0,2*np.pi,np.pi/100),np.arange(0,2*np.pi,np.pi/100))
zimg = np.array([np.sin(alpha)*np.sin(wt),
                 -(np.sin(beta)*np.cos(wt)*np.cos(alpha)+np.sin(alpha)*np.cos(beta)),
                 -np.sin(beta)*np.cos(wt)*np.sin(alpha)+np.cos(alpha)*np.cos(beta)
                 ])
lam0 = np.atan(np.linalg.norm(zimg,axis=0)*np.sin(phi))
plt.plot(phi[0,:],lam0[0,:],)
plt.plot(phi[24,:],lam0[24,:])
plt.plot(phi[49,:],lam0[49,:])
plt.plot(phi[74,:],lam0[74,:])
plt.legend([r'$\omega t=0$',r'$\omega t=pi/2$',r'$\omega t=pi$',r'$\omega t=3pi/2$'])
plt.xlabel(r'$\varphi$ / rad')
plt.ylabel(r'$\lambda$ / rad')
plt.show()

