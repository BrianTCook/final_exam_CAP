from __future__ import division
import numpy as np
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk

def initial_conditions(Q, N, Rmax, diskmassfrac, seed_disk):

	'''
	returns Sun + 8 planets, as well as a disk around Neptune
	'''

	#get the solar system objects

	sun_and_planets = new_solar_system()

	#throw out Pluto, get Neptune object from new_solar_system
	N_particles = len(sun_and_planets)
	sun_and_planets = sun_and_planets[:N_particles-1]
	Nep_ind = len(sun_and_planets)-1
	Sun = sun_and_planets[0]
	Neptune = sun_and_planets[Nep_ind]

	#hill radius of Neptune
	a_Neptune = np.sqrt(Neptune.x**2 + Neptune.y**2 + Neptune.z**2)
	Neptune_Hill = a_Neptune * (Neptune.mass/Sun.mass)**(1/3.)

	#get the disk with appropriate Salpeter mass distribution

	#Salpeter mass function with a range of 1 order of magnitude
	mmin = 0.0005*Neptune.mass #mmax goes to 0.005
	masses = new_salpeter_mass_distribution(N, mass_min=mmin, mass_max=(10*mmin))

	#set up converter, protoplanetary disk gas particles
	np.random.seed(seed_disk) #random seed for the disk
	converter_gas = nbody_system.nbody_to_si(Neptune.mass, Rmax|units.AU)
	gas = ProtoPlanetaryDisk(N, convert_nbody=converter_gas, Rmin=converter_gas.to_nbody(Neptune_Hill).number, Rmax=converter_gas.to_nbody(Rmax|units.AU).number, q_out=Q, discfraction=diskmassfrac).result

	#attribute Salpeter masses to the gas
	gas.mass = masses

	#move disk to Neptune's phase space coordinates
	gas.x += Neptune.x
	gas.y += Neptune.y
	gas.z += Neptune.z
	gas.vx += Neptune.vx
	gas.vy += Neptune.vy
	gas.vz += Neptune.vz

	return sun_and_planets, gas
