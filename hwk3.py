import math
def dotp (a,b) -> float:
    if len(a) != len(b):
        return float(inf)
# m_pulsar = input("mass of pulsar (g) ")
# v0 = input(" initial velocity (cm/s) ")
# r0 = input ("maximum radius of magnetic field line (cm) ")
# r_pulsar = input("radius of pulsar (cm) ")
# dt = input("smallest timestep (s) ")
m_pulsar = 2.864e33 # 1.4 Msun
r0= 4e6 # radius of mag-field
r_pulsar = 1e6 # pulsar inner radius
v0 = 1e5 # 10^5 cm/s
dt = 1e-6 # s
f = 500 # Hz
omega = 2*math.pi*f #rad Hz
r = r0 # initialise
v = v0 # initialise
t = 0 
lam = 0
chi = math.pi / 2
dx = 0
G = 6.67e-8 # dyn cm^2 g^-2
while r >= r_pulsar:
    a= G*m_pulsar * math.cos(chi) / (r**2) - (omega**2 * r * math.cos(lam) * math.cos(chi-lam))  # calculate gravity and centrifugal force (centrifugal force always acts in opposite direction to gravity)
    dx = v*dt + 0.5*a*(dt**2) # verlet
    v = v + a*dt # update velocity 
    lam += dx / (r0*math.cos(lam) * (1+3*(math.sin(lam)**2))**0.5) # calculate delta-lambda based on eq. given
    chi = math.atan(1/(2*math.tan(lam)))
    r = r0 * (math.cos(lam)**2)
    t += dt
print (t,v,r)
omega**2 * r0 * math.cos(lam)**3 * math.cos(chi-lam)
