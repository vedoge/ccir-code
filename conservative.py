import numpy as np
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
M = 2.784e33 # g
rm = 1e8 # cm
rpuls = 1e6  # cm
v0 = 1e7 # cm s^-1

dt = 1e-4 # s
t = 0
lam = 0.2 # rad
chi = 0 # rad
v = v0
phi = 0
# compare specific energies of the particle at the start and at the end 
# now integrate
while rm*cos(lam)**2 >= rpuls:
	print(lam)
	chi = atan(0.5/tan(lam))
	a = G*M*cos(chi)/((rm*(cos(lam)**2))**2)
	ds = v*dt + 0.5*a*(dt**2) 
	v += a*dt
	lam += ds/(rm*cos(lam)*sqrt(1 + 3*(sin(lam)**2)))
rf = rm*(cos(lam)**2)
#v -= a*dt
EF = 0.5*(v**2) - G*M/(rf)
E0 = 0.5*(v0**2) - G*M/(rm*cos(0.2)**2)
print(EF-E0)
