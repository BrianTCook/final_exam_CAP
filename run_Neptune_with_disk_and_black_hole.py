from __future__ import division, print_function
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import operator
import time
from amuse.lab import *
from amuse.couple import bridge
from amuse.units.optparse import OptionParser

from make_initial_conditions import initial_conditions
from add_black_hole import BH_maker

def iterate_through_simulations(Q, N, seed_disk, Rmax, diskmassfrac, BH_mass, x0_BH, y0_BH, vx0_BH, t_end, dt, BH_bool):

	'''
	was having an indentation error, so this seemed
	like the easiest way to toggle the BH being in the simulation
	'''

	if BH_bool == True:

		sun_and_planets, gas = initial_conditions(Q, N, Rmax, diskmassfrac, seed_disk)
		BH = BH_maker(BH_mass, x0_BH, y0_BH, vx0_BH)

		N_objects = len(sun_and_planets)
		sun_and_planets_minus_Neptune = sun_and_planets[:N_objects-1]
		Neptune = sun_and_planets - sun_and_planets_minus_Neptune

		#set up two gravity solvers, external to Neptune system and neptune

		external_bodies = sun_and_planets_minus_Neptune
		converter_external = nbody_system.nbody_to_si(external_bodies.mass.sum(), 500.|units.AU)
		gravity_external = ph4(converter_external)
		gravity_external.particles.add_particles(sun_and_planets_minus_Neptune)
		gravity_external.particles.add_particles(BH)

		mass_internal = Neptune.mass.sum() + gas.mass.sum()
		converter_internal = nbody_system.nbody_to_si(mass_internal, Rmax|units.AU)
		gravity_internal = BHTree(converter_internal)
		gravity_internal.particles.add_particles(Neptune)
		gravity_internal.particles.add_particles(gas)

		gravity = bridge.Bridge()
		gravity.add_system(gravity_internal, (gravity_external,))
		gravity.add_system(gravity_external, (gravity_internal,))
		gravity.timestep = dt|units.yr

		channel_from_gravity_to_sun_and_planets = gravity.particles.new_channel_to(sun_and_planets)
		channel_from_gravity_to_gas = gravity.particles.new_channel_to(gas) 
		channel_from_external_to_BH = gravity.particles.new_channel_to(BH)

		sun = sun_and_planets[0]
		planets = sun_and_planets[1:]

		plt.figure()
		plt.scatter(sun.x.value_in(units.AU), sun.y.value_in(units.AU), s=4, c='y')
		plt.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU), s=2, c='b')
		plt.scatter(gas.x.value_in(units.AU), gas.y.value_in(units.AU), s=1, c='k')
		plt.scatter(BH.x.value_in(units.AU), BH.y.value_in(units.AU), s=2, c='k')

		plt.xlim(-600, 600)
		plt.ylim(-600, 600)
		plt.xlabel('AU')
		plt.ylabel('AU')
		plt.title('Zero time configuration')
		plt.savefig('zero_time_with_BH.png')
		plt.close()

		gravity.evolve_model(t_end|units.yr)
		channel_from_gravity_to_sun_and_planets.copy()
		channel_from_gravity_to_gas.copy()
		channel_from_external_to_BH.copy()

		sun = sun_and_planets[0]
		planets = sun_and_planets[1:]

		plt.scatter(sun.x.value_in(units.AU), sun.y.value_in(units.AU), s=4, c='y')
		plt.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU), s=2, c='b')
		plt.scatter(gas.x.value_in(units.AU), gas.y.value_in(units.AU), s=1, c='k')
		plt.scatter(BH.x.value_in(units.AU), BH.y.value_in(units.AU), s=2, c='k')
		plt.xlim(-600, 600)
		plt.ylim(-600, 600)
		plt.xlabel('AU')
		plt.ylabel('AU')
		plt.title('Final Time: %.02f yr'%(t_end))
		plt.savefig('evolved_with_BH.png')
		plt.close()

		gravity.stop()

	if BH_bool == False:

		sun_and_planets, gas = initial_conditions(Q, N, Rmax, diskmassfrac, seed_disk)

		N_objects = len(sun_and_planets)
		sun_and_planets_minus_Neptune = sun_and_planets[:N_objects-1]
		Neptune = sun_and_planets - sun_and_planets_minus_Neptune

		#set up two gravity solvers, external to Neptune system and neptune

		external_bodies = sun_and_planets_minus_Neptune
		converter_external = nbody_system.nbody_to_si(external_bodies.mass.sum(), 500.|units.AU)
		gravity_external = ph4(converter_external)
		gravity_external.particles.add_particles(sun_and_planets_minus_Neptune)

		mass_internal = Neptune.mass.sum() + gas.mass.sum()
		converter_internal = nbody_system.nbody_to_si(mass_internal, Rmax|units.AU)
		gravity_internal = BHTree(converter_internal)
		gravity_internal.particles.add_particles(Neptune)
		gravity_internal.particles.add_particles(gas)

		gravity = bridge.Bridge()
		gravity.add_system(gravity_internal, (gravity_external,))
		gravity.add_system(gravity_external, (gravity_internal,))
		gravity.timestep = dt|units.yr

		channel_from_gravity_to_sun_and_planets = gravity.particles.new_channel_to(sun_and_planets)
		channel_from_gravity_to_gas = gravity.particles.new_channel_to(gas) 

		sun = sun_and_planets[0]
		planets = sun_and_planets[1:]

		plt.figure()
		plt.scatter(sun.x.value_in(units.AU), sun.y.value_in(units.AU), s=4, c='y')
		plt.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU), s=2, c='b')
		plt.scatter(gas.x.value_in(units.AU), gas.y.value_in(units.AU), s=1, c='k')

		plt.xlim(-50, 50)
		plt.ylim(-50, 50)
		plt.xlabel('AU')
		plt.ylabel('AU')
		plt.title('Zero time configuration')
		plt.savefig('zero_time_without_BH.png')
		plt.close()

		gravity.evolve_model(t_end|units.yr)
		channel_from_gravity_to_sun_and_planets.copy()
		channel_from_gravity_to_gas.copy()

		sun = sun_and_planets[0]
		planets = sun_and_planets[1:]

		plt.scatter(sun.x.value_in(units.AU), sun.y.value_in(units.AU), s=4, c='y')
		plt.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU), s=2, c='b')
		plt.scatter(gas.x.value_in(units.AU), gas.y.value_in(units.AU), s=1, c='k')
		plt.xlim(-50, 50)
		plt.ylim(-50, 50)
		plt.xlabel('AU')
		plt.ylabel('AU')
		plt.title('Final Time: %.02f yr'%(t_end))
		plt.savefig('evolved_without_BH.png')
		plt.close()

		gravity.stop()

def main(Q, N, seed_disk, Rmax, diskmassfrac, BH_mass, x0_BH, y0_BH, vx0_BH, t_end, dt):

	iterate_through_simulations(Q, N, seed_disk, Rmax, diskmassfrac, BH_mass, x0_BH, y0_BH, vx0_BH, t_end, dt, True)
	iterate_through_simulations(Q, N, seed_disk, Rmax, diskmassfrac, BH_mass, x0_BH, y0_BH, vx0_BH, t_end, dt, False)


def new_option_parser():
	
	'''
	define an option parser for initial conditions
	'''

	optparser = OptionParser()
	optparser.add_option('--Q', dest='Q', type='float', default=1.)
	optparser.add_option('--N', dest='N', type='int', default=1000)
	optparser.add_option('--seed_disk', dest='seed_disk', type='int', default=1)
	optparser.add_option('--Rmax', dest='Rmax', type='float', default=0.03)
	optparser.add_option('--diskmassfrac', dest='diskmassfrac', type='float', default=0.01)
	optparser.add_option('--BH_mass', dest='BH_mass', type='float', default=100.)
	optparser.add_option('--x0_BH', dest='x0_BH', type='float', default=-500.)
	optparser.add_option('--y0_BH', dest='y0_BH', type='float', default=200.)
	optparser.add_option('--vx0_BH', dest='vx0_BH', type='float', default=30.)
	optparser.add_option('--t_end', dest='t_end', type='float', default=10.)
	optparser.add_option('--dt', dest='dt', type='float', default=0.1)

	return optparser

if __name__ in '__main__':
	o, arguments = new_option_parser().parse_args()
	main(**o.__dict__)
