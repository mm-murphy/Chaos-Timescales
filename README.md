# Chaos-Timescales

csim.py is a REBOUND simulation code that I wrote up as part of my Stony Brook Honors thesis (2019-2020). It is an N-body that simulates a planetary system out to a specified amount of orbits. The result is tracking the evolution of eccentricity and angular momentum deficit (AMD) of each planet. We use this simulation to identify secular planetary systems and identify secular timescales, whether they be short-term periodicity in the eccentricity or long-term chaos. 

The python program takes command line arguments as the planetary system's initial parameters. It is run as:

python csim.py sysid mass initial_e alpha a1set a2set a3set f1set f2set f3set

Here, 
sysid = a string identifier for the particular system (ex: sys0)
mass = the mass of each planet in M_sol units (this code treats an equal mass system)
initial_e = the initial eccentricity (treats all planets having the same value)
alpha = separation factor ( 0 < alpha < 1 )
a1,2,3set = initial semi-major axes of each planet, as computed from the separation factor
f1,2,3set = initial orbital angles of each planet

The program will output a .txt file that is named from the passes sysid argument. This file contains the system parameters over time in the order:

time; planet 1 eccentricity; planet 2 eccentricity; planet 3 eccentricity; planet 1 AMD, planet 2 AMD, planet 3 AMD
