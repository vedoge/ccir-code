import numpy as np
import matplotlib.pyplot as plt
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

G = 6.67e-8 # dyn cm^2 g^-2
M = 2.864e33 # 1.4 Msun
c = 2.998e10 # cm s^-1
rs = 2*G*M/(c**2)
rm= 1e8 # radius of mag-field
rpuls = 1e6 # pulsar inner radius
maxlam = np.acos(np.sqrt(rpuls/rm))
dphi = -1e-3
u = 0.5*rs/rpuls # for all of them
alpha, phi0 = np.meshgrid(np.arange(0,2*pi,pi/100),np.arange(0,2*pi,pi/100))
phi = phi0
uv = -u/tan(alpha-phi)
u = np.full_like(phi,u)
uf = np.zeros_Like(u)
ufprime = np.zeros_like(uv)
fphi = np.zeros_like(phi)
ul = np.expand_dims(u.copy(),axis=2)
stopped = np.full_like(phi,False,dtype=bool)
while not np.all(stopped):
	ua = 3*(u**2) - u
	u += uv*dphi + 0.5*ua*(dphi**2)
	uv += ua*dphi
	phi += dphi # check
	ul = np.append(ul[:,:,:],u[:,:,np.newaxis],axis=2)
	print(ul.shape[2],np.sum(stopped))
	# correct fate: 
	# idx = [condition]
	ulim = 0.5*rs/(rm*(sin(phi)**2)) # u limits
	# whatever has cleared the magnetosphere
	idx = ((ul[:,:,-1] >= ulim) & (ul[:,:,-2] <= ulim)) | ((ul[:,:,-1] <= ulim) & (ul[:,:,-2] >= ulim))
	# process the impact angle
	stopped[idx] = True
	# record du/dphi, record u, record phi for the impact
	uf[idx] = u[idx]
	ufprime[idx] = uv[idx]
	fphi[idx] = phi[idx]
	u[idx] = 0 # photon banished to infinity
	uv[idx] = 0
	idx =  (u <= 0.5*rs/rm) | (u >= 0.5*rs/rpuls)
	stopped[idx] = True
	u[idx] = 0
	uv[idx] = 0
# use the formula in the paper to calculate dy/dx for each of the curves
rf = 0.5*rs/uf
rprime = -0.5*rs/(ufprime**2)
dydx_l = (rprime*cos(phi)-rf*sin(phi))/(rprime*sin(phi)+rf*cos(phi))
rprime_b = rm*sin(2*phi)
dydx_b = (rprime_b*cos(phi)-rf*sin(phi))/(rprime_b*sin(phi)+rf*cos(phi))
dot = cos(atan(dydx_l) + atan(dydx_b)
