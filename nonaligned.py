import numpy as np
import matplotlib.pyplot as plt
alpha = 1 #rad
beta = 0.5 #rad
phi,wt = np.meshgrid(np.arange(0,2*np.pi,np.pi/100),np.arange(0,2*np.pi,np.pi/100))
zimg = np.array([np.sin(alpha)*np.sin(wt),
                 -(np.sin(beta)*np.cos(wt)*np.cos(alpha)+np.sin(alpha)*np.cos(beta)),
                 -np.sin(beta)*np.cos(wt)*np.sin(alpha)+np.cos(alpha)*np.cos(beta)
                 ])
print(zimg.shape)
lam0 = np.atan(np.linalg.norm(zimg,axis=0)*np.sin(phi))
plt.plot(phi[0,:],lam0[0,:],)
plt.plot(phi[24,:],lam0[24,:])
plt.plot(phi[49,:],lam0[49,:])
plt.plot(phi[74,:],lam0[74,:])
plt.legend(['wt=0','wt=pi/2','wt=pi','wt=3pi/2'])
plt.xlabel('phi')
plt.ylabel('lambda')
plt.show()

