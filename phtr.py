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
sign = np.sign

G = 6.67e-8 # dyn cm^2 g^-2
M = 2.864e33 # 1.4 Msun
c = 2.998e10 # cm s^-1
rs = 2*G*M/(c**2)
rm= 1e8 # radius of mag-field
rpuls = 1e6 # pulsar inner radius
maxlam = np.acos(np.sqrt(rpuls/rm))
inc = -1e-4
alpha, lam0 = np.meshgrid(np.arange(-pi/2,pi/2,pi/100),np.array([-pi/2,pi/2]))
lam = lam0
print(lam.shape)
u = np.full_like(lam,0.5*rs/rpuls)
uv = -u/tan(alpha-lam)
# result collection arrays
uf = np.full_like(u,np.nan)
ufprime = np.full_like(uv,np.nan)
flam = np.full_like(lam,np.nan)
uprev = u.copy()
# keep track of which light rays are still moving
stopped = np.full_like(lam,False,dtype=bool)
try:
	while not np.all(stopped):
		ua = 3*(u**2) - u
		dlam = sign(uv)*inc
		#dlam = inc
		uprev = u
		u += uv*dlam + 0.5*ua*(dlam**2)
		uv += ua*dlam
		lam += dlam # check
		print(np.sum(~stopped))
		# correct fate: 
		# idx = [condition]
		# whatever has cleared the magnetosphere
		idx = (u >= 0.5*rs/(rm*(cos(lam)**2))) & (uprev <= 0.5*rs/(rm*(cos(lam-dlam)**2)))
		# process the impact angle
		stopped[idx] = True
		uf[idx] = u[idx]
		ufprime[idx] = uv[idx]
		flam[idx] = lam[idx]
		u[idx] = 0 # photon banished to infinity
		uv[idx] = 0
		idx =  (u <= 0.5*rs/rm) | (u >= 0.5*rs/rpuls)
		stopped[idx] = True
		u[idx] = 0
		uv[idx] = 0
	raise KeyboardInterrupt
except KeyboardInterrupt:
	rf = 0.5*rs/uf
	rprime = -0.5*rs*ufprime/(ufprime**2)
	# checked with Sharada Aunty
	# think about the below formula. Check it: with conversion with dy/dx, it should vanish at 35.3 deg. 
	rprime_b = -rm*sin(2*flam)
	gamma = atan((1/rprime_b)*rf) - atan((1/rprime)*rf)
	print(gamma[0,:])
	print(flam[0,:])
	print(maxlam,pi-maxlam)
	plt.plot(flam[0,:],gamma[0,:])
	plt.show()
	plt.plot(flam[1,:], gamma[1,:])
	plt.show()
# now let's think about how to extract this data and turn it into force
# my method - calculate the solid angle attached to each light ray
# calculate the total luminous exitance at the poles (erg s^-1)
# divide by 4pi (2 hemispheres of radiation - together 4pi steradians - unit area) to find radiance (erg s^-1 cm^-2 sr)
# multiply by cos(alpha) - to account for distribution
# then trace light paths
# find the cells of impact
# apply the formulae, get the analytical solutions

