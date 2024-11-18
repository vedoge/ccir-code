# montecarloroutines.py
# Monte Carlo Routines for finalmontecarlo.py. 
if __name__ == '__main__':
	print("do not run this file directly")
#	raise EOFError
# operations on python variables (not numpy variables) are serialized, so we can use this for coordination...? Utilize the GIL? 
# We need a unified way of computing the mass accretion rate. 
# Currently, my showerthoughts are: 
# global dt
# each thread does the physics, computes the particles that have hit the surface
# mass accretion rate
# rad rate
# use rad map to update acc map
# region is a 2x2 array indexing 

def update_region(region, dt):
	global lam, phi
	global a,plam,pphi,pv
	lam_parts = plam.reshape(1,1,size(plam))
	phi_parts = pphi.reshape(1,1,size(pphi))
	# use the 2nd axis
	# use an Event object to signal (releases GIL when thread is waiting)
	# minimize each axis individually
	lamcell = np.argmin(np.abs(lam_parts - lam[0,:,np.newaxis]))
	phicell = np.argmin(np.abs(phi_parts - phi[:,0, np.newaxis]))
	print(lamcell,phicell)
