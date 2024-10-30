import numpy as np
import matplotlib.pyplot as plt
G = 6.67e-8 # dyn cm^2 g^-2
M = 2.864e33 # 1.4 Msun
c = 2.998e10 # cm s^-1
rs = 2*G*M/(c**2)
rm= 1e8 # radius of mag-field
rpuls = 1e6 # pulsar inner radius
maxlam = np.acos(np.sqrt(rpuls/rm))
dphi = -1e-4
# things to do:
# t against phi
# constant ds solver
# find the distribution of the photons (TODO Thurs)
# calculate particle motion from the accretion disk (TODO Thurs PM)
# difference between freefall and final velocity (very precise simulation, no rotation, no precession, no radiation)
# pulsars are very important
# as they provide precise time signals -> testing general relativity
# eqn works (tested with matlab + ode45)
# any further problems r to do w/ code
u = 0.5*rs/rpuls # for all of them
alpha,phi = np.meshgrid(np.arange(0,2*np.pi,np.pi/1000),np.pi/2)
uv = -u/np.tan(alpha - phi)
u = np.full_like(uv,u)
r = np.array([0.5*rs/u])
endlams = np.full(alpha.shape, np.nan)
weights = np.zeros(alpha.shape)
oldu = np.full(alpha.shape,np.nan)
rprime = np.full(alpha.shape,np.nan)
i = 1 # iteration counter
total = 0 # nr of geodesics that have left the pulsar
# running without the magnetosphere blocking the light, it seems that the oblique angles enable light to almost succumb to the curvature, circling multiple times before coming out. 
while np.any(u != 0):
	ua = 3*(u**2) - u
	u += uv*dphi + 0.5*ua*(dphi**2)
	uv += ua*dphi
	phi += dphi
	idx = (u >= 0.5*rs/(rm*(np.cos(phi))**2))
	# record the phi coordinate
	endlams[idx] = phi[idx] # record ending lambda / phi
	# now process the velocities
	rprime[idx] = uv[idx]*(-0.5*rs/(u[idx]**2)) # by chain rule 
	weights[idx] += np.cos(np.atan(rprime[idx]) - np.atan(0.5/np.tan(phi[idx])) + np.pi/2)*np.cos(alpha[idx])  # dot product (cos alpha weight of distribution)	
	# alternate photon fates
	u[idx] = 0
	uv[idx] = 0
	u[u <= 0.5*rs/rm] = 0 # photon exits NS sphere of influence
	u[u >= 0.5*rs/rpuls] = 0 # photon impacts NS surface 
	uv[u==0] = 0 # mark all coordinates which have ended here as processed
	i += 1
'''
plt.contourf(alpha,phi,weights)
plt.ylabel(r"$\lambda_0$ / rad")
plt.xlabel(r"$\alpha$ / rad")
plt.colorbar()
plt.title(r"plot of weight against $\alpha$ and $\lambda_0$")
'''
plt.plot(endlams[np.isnan(endlams) == False],weights[np.isnan(endlams) == False])
plt.show()

