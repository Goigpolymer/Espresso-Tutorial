
# Copyright (C) 2013-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
This sample simulates a Lennard-Jones fluid maintained at a fixed temperature
by a Langevin thermostat:
1.-SET UP A SIMPLELENNARD JONES SIMULATION 2D :SAVING FILES THERMODYNAMIC PROPERTIES AND FINAL CONFIGURATION
2.-CARRY OUT SIMPLE VISUALIZATION USING THE OPENGL LIBRARY
3.-STORE A TRAJECTORY FILE USING VSF AND VTF FILES READEABLE WITH VMD VISUALIZATION PACKAGE 
4.-Store a final coordinate file to carry out post simulation analysis using Freud-Analysis package.
"""
#Import python libraries
from __future__ import print_function
from espressomd.visualization_opengl import *

import numpy as np
import espressomd
from espressomd import visualization
from espressomd.io.writer import vtf
import MDAnalysis as mda


#Print Espresso features installed
required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

#Import thermodynamic thermostats
from espressomd import thermostat
from random import uniform


print("""
=======================================================
=                    lj_liquid2d.py: basic set up                    =
=======================================================

Program Information:""")
print(espressomd.features())

dev = "cpu"

#Define simulation parameters
# System parameters
#############################################################

#box_l = 30
#density = 0.7
#box=30

n_part=1000
density=0.20
box_l=np.sqrt(n_part/density)
box=box_l;
#
# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5 * lj_sig
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)
system.time_step = 0.001
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
visualizer=espressomd.visualization.openGLLive(system,window_size=[500,500])


# Warmup parameters
#############################################################

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min_dist
min_dist = 0.9

# Integration parameters
int_steps = 2000
int_n_times = 100

#1) SET UP SIMULATION PARAMETERS
#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup (Non-bonded Interactions)
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig,cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())


# Particle setup: Set up N particles randomly across the simulaton domain
#########################################################################

area = box * box
n_part = int(area * density)


a1=np.zeros(n_part)
a2=np.zeros(n_part)
a3=np.zeros(n_part)


for i in range(n_part):
    a1[i]=0.49*box*(uniform(0,1))
    a2[i]=0.49*box*(uniform(0,1))
    a3[i]=0.5*box


for i in range(n_part):
    system.part.add(pos=[a1[i],a2[i],a3[i]],id=i,type=0,fix=[0,0,1])


system.analysis.dist_to(0)


#Printing: the number of particles, box size and the density
print("Simulate {} particles in a cubic simulation box of length {} at density {}."
      .format(n_part, box_l, density).strip())
      
#Printing: warming up initial start      
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print("Start with minimal distance {}".format(act_min_dist))

#Number of Linked cells
system.cell_system.max_num_cells = 2744

#############################################################
#  Warmup Integration   stage                                     #
#############################################################

# open Observable file: this file store basic thermodynamic information about the system
#Simulation time,total Energy, kinetic Energy and Potential Energy.
obs_file = open("pylj_liquid.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\n")

#Printing warming up parameters.
print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# 2)WARMING UP STAGE
# Warmup Integration Loop: rutine warmingup Lj particles.
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(steps=warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.min_dist()
    i += 1
#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap = lj_cap

# Just to see what else we may get from the c code
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
# pprint.pprint(system.part.__getstate__(), width=1)
state = system.__getstate__()
pprint.pprint(state)

# write parameter file

# polyBlockWrite "$name$ident.set" {box_l time_step skin} ""
set_file = open("pylj_liquid.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" %
               (box, system.time_step, system.cell_system.skin))

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# print(initial energies)
energies = system.analysis.energy()
print(energies)

fp=open('LJ.vtf',mode='w+t')
vtf.writevsf(system,fp)
vtf.writevsf(system,fp)


#Main integration Loop
j = 0
for i in range(int_n_times):
    print("run %d at time=%f " % (i, system.time))
    system.integrator.run(steps=int_steps)
    energies = system.analysis.energy()
    print(energies)
    obs_file.write('{ time %s } %s\n' % (system.time, energies))
    linear_momentum = system.analysis.linear_momentum()
    print(linear_momentum)
    visualizer.screenshot(f'screenshot{i:0>5}.png')
    vtf.writevcf(system,fp) 
    sys.stdout.flush()
fp.close( )


# write end configuration
end_file = open("pylj_liquid.end", "w")
end_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
end_file.write("{ particles {id pos type} }")

#Save finalconfiguration
for i in range(n_part):
    end_file.write("%s\n" % system.part[i].pos)
    

obs_file.close()
set_file.close()
end_file.close()

#Store post analysis.file

box2=-0.5+0.5*box
box3=0.5+0.5*box

A=[[0,box],[0,box],[box2,box3]]
print(A)
#A=np.array(Boxsize)
np.savetxt('Box.txt',A)

C=system.part[:].pos       
B=np.array(C)
np.savetxt('LJ.txt',B)




#
# ========================================================="
# Optional MDAnalysis moudle        "
# quantities from MDAnalysis                               "
# ========================================================="
#

#eos = MDA_ESP.Stream(system)

#u = mda.Universe(eos.topology, eos.trajectory)


# let's have a look at the universe
#print(u)

# Inspect atoms
#print(u.atoms)

#print("Positions:")
#print(u.atoms.positions)
#print("Velocities:")
#print(u.atoms.velocities)
#print("Forces:")
#print(u.atoms.forces)
#print("Names:")
#print(u.atoms.names)
#print("IDs:")
#print(u.atoms.ids)
#print("Types:")
#print(u.atoms.types)
#print("Charges:")
#print(u.atoms.charges)

#Writting a PDB file
#u.atoms.write("system.pdb")
#print("===> The initial configuration has been written on system.pdb ")

#from MDAnalysis.analysis.rdf import InterRDF

#charged = u.select_atoms("prop charge  > 0")
#rdf = InterRDF(charged, charged, nbins=7, range=(0, 10))

#RDFcalculation
# This runs so far only over the single frame we have loaded.
# Multiframe averaging must be done by hand

#rdf.run()



# terminate program
print("\nFinished.")
visualizer = openGLLive(system, background_color=[1, 1, 1])
visualizer.run(1)

