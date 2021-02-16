#!/usr/bin/env python
# -*- coding: utf-8 -*-

# cms for gravML project
# contact ch.schaefer@uni-tuebingen.de for questions....

# using gravitational_harmonics https://reboundx.readthedocs.io/en/latest/effects.html 
# to model non-spherical central body which is orbited by smaller object
# paper for the implementation https://academic.oup.com/mnras/article/491/2/2885/5594029


# ChangeLog:
# 2021-02-16 first implementation

import numpy as np
import rebound
import reboundx
import matplotlib.pyplot as plt

# standard units, G=1
# a = 1.0
# -> T = 2\pi
# total number of orbits
Norbits = 10000
tmax = Norbits * 2*np.pi
# samples per orbit
Nsamples = 1
# output number
Nout = Norbits*Nsamples

# one big object, one orbiting object
sim = rebound.Simulation()
sim.add(m=1.0, hash='central_object')
sim.add(m=1e-5, a=1, e=0.2, hash='orbiting_object')
sim.move_to_com()

# now add J2 and J4 to central object
# load reboundx gravitational_harmonnics force
rebx = reboundx.Extras(sim)
gh = rebx.load_force("gravitational_harmonics")
rebx.add_force(gh)

# set J2
sim.particles['central_object'].params['J2'] = 0.1
# set equatorial radius
sim.particles['central_object'].params['R_eq'] = 0.1
# we set J4 to 0.0 for the time being
sim.particles['central_object'].params['J4'] = 0.0

# integration, plotting, and output

times = np.linspace(0, tmax, Nout)
x = np.zeros((2, Nout))
y = np.zeros_like(x)
for i, t in enumerate(times):
    sim.integrate(t)
    print("Current time is %-05.03f, %-02.02f %% done \t\t\t\t\t\t\t\t" % (t, t/tmax*1e2),  end='\r')
    for j, p in enumerate(sim.particles):
        x[j,i] = p.x
        y[j,i] = p.y
        

# create nice pic of orbit
fig, ax = plt.subplots()
ax.scatter(x[0,:], y[0,:], s=1, c='y')
ax.scatter(x[1,:], y[1,:], s=1, c='b')
ax.set_aspect('equal')
ax.grid()
ax.set_xlabel(r'x-coordinate')
ax.set_ylabel(r'y-coordinate')
fig.savefig("orbit.pdf")
plt.close()

# write down coordinates of orbiting object with time tag
np.savetxt("gravML_orbits.txt", np.transpose([times, x[1,:], y[1,:]]))

# in case you want to restart the sim
sim.save("save.bin")




