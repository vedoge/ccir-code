# montecarloroutines.py
# Monte Carlo Routines for montecarlo.py. 
if __name__ == '__main__':
    print("do not run this file directly")
    raise EOFError
# we need to record individual times for each thread
# The mass accretion method will be done with a synchronize after each thread step
# **kwargs dt passing?
# operations on python variables (not numpy variables) are serialized, so we can use this for coordination...? Utilize the GIL? 
# We need a unified way of computing the mass accretion rate. 
# Currently, my showerthoughts are: 
# global dt
# each thread does the physics, computes the particles that have hit the surface
# mass accretion rate
# rad rate
# use rad map to update acc map
# profit
def update_region(region, particle_lams, particle_phis):

