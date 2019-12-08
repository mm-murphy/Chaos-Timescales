### This is a program for tracking the secular evolution of a system of planets
##     It tracks eccentricity and the AMD held by each planet over time
##  Author: Matthew M Murphy - undergraduate - Stony Brook University
##  Last updated: Dec. 2019
######################################################################################

import rebound                  # I use the REBOUND package, see its online documentation
import numpy as np
import matplotlib.pyplot as plt
import sys

# must first input parameters from command line
sysid = sys.argv[1]              # input system number
mass = float(sys.argv[2])         # planet mass
initial_e = float(sys.argv[3])    # initial eccentricity
alpha = float(sys.argv[4])        # alpha spacing parameter
a1set = float(sys.argv[5])        # semi-major axis
a2set = float(sys.argv[6])
a3set = float(sys.argv[7])
f1set = float(sys.argv[8])        # mean-longitude
f2set = float(sys.argv[9])
f3set = float(sys.argv[10])

# Setting initial system values
# Currently is a 3-planet system
m_1 = mass       # mass of innermost planet
r_1 = 0.000477   # radius of innermost planet
m_2 = m_1        # 2,3 are the other planets moving outward
r_2 = r_1
m_3 = m_2
r_3 = r_2

e_i = initial_e
f_1 = f1set  
f_2 = f2set  
f_3 = f3set    

def Sim():      # Simulation setup
  sim = rebound.Simulation()
  sim.units = ('yr','AU','Msun')
  sim.add(m=1.,r=0.005)                   # solar-like host star
  a1 = a1set                              # separations based on  
  a2 = a2set                              #   alpha = a_inner / a_outer
  a3 = a3set
  sim.add(m=m_1,a=a1,e=e_i,f=f_1,r=r_1)   # zero inclination orbits
  sim.add(m=m_2,a=a2,e=e_i,f=f_2,r=r_2) 
  sim.add(m=m_3,a=a3,e=e_i,f=f_3,r=r_3)
  sim.integrator = 'whfast'                # Using the IAS15 integrator
  sim.dt = 0.03
  sim.ri_whfast.safe_mode=0
  sim.ri_whfast.corrector = 11
  sim.move_to_com()
  return sim

def AMD(m,a,e):          # Canonical AMD definition
  Lambda = m*np.sqrt(sim.G*a*(m+1))/(m+1)
  return Lambda*(1-np.sqrt(1-e**2))

def EncTest(p,N,time):   # Check for orbit crossing, close encounter, or collision
  for j,planet in enumerate(p):
          if j>1:
                  sep = np.sqrt((p[j].x-p[j-1].x)**2 + (p[j].y-p[j-1].y)**2 + (p[j].z-p[j-1].z)**2)
                  Rhillmut = pow((p[j].m+p[j-1].m)/3.,1./3)*((p[j].a + p[j-1].a)/2.)
                  if (p[j].a*(1-p[j].e) < p[j-1].a*(1+p[j-1].e)):
                          return 'true'
                  elif abs(sep) < Rhillmut:
                          return 'true'
  if N != 3:
          return 'true'
  else:
          return 'none'
  
# Now, to finally run the simulation
################################################################
sim = Sim()              # Sets up the system
p_init = sim.particles   # Array of initial condition
inner_period = p_init[1].P

tmax = 1.e4*inner_period                        # Maximum integration time 

Nintsteps = int(tmax / (1*inner_period))
int_times = np.linspace(0,tmax,Nintsteps)

################################################################

sim.collision = 'direct'          # Simulation treats any collisions between planets as a direct
sim.collision_resolve = 'merge'   #   collision and will merge them, conserving momentum and energy

# Creating a data file for output
data_name = str(sysid)+'data'
data_title = data_name+'.txt'
dfile = open(data_title,'w')


# Now running the integration

for i,time in enumerate(int_times):
	sim.integrate(time, exact_finish_time=0)
	p = sim.particles
	check = EncTest(p, sim.N - 1, time)
	if check == 'true':   
		break          
	A1 = AMD(p[1].m,p[1].a,p[1].e)  
	A2 = AMD(p[2].m,p[2].a,p[2].e) 
	A3 = AMD(p[3].m,p[3].a,p[3].e)
	Atot = A1 + A2 + A3
	if (i%1000 == 0):
		wo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{} '.format(time,p[1].e,p[2].e,p[3].e,A1,A2,A3,Atot)
		dfile.write(wo)
		dfile.write('\n')

dfile.close()
