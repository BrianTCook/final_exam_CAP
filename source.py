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

class GalacticCenterGravityCode(object):
	def __init__(self, R, M, alpha):
		self.radius = R
		self.mass = M
		self.alpha = alpha

	def get_gravity_at_point(self, eps, x, y, z):
		
		r2 = x**2 + y**2 + z**2
		r = np.sqrt(r2)

		scaled_radius = r/self.radius

		m = self.mass*(scaled_radius)**self.alpha

		fr = constants.G*m/r2
		ax = -fr*x/r
		ay = -fr*y/r
		az = -fr*z/r

		return ax, ay, az

	def circular_velocity(self, r):
		
		scaled_radius = r/self.radius	
		m = self.mass*(scaled_radius)**self.alpha
		vc = np.sqrt(constants.G*m/r)

		return vc

#original wasn't working properly
def LagrangianRadii(stars, massf):

	com = stars.center_of_mass()
	vcom = stars.center_of_mass_velocity()
	n_stars = len(stars)

	d2 = (stars.position - com).lengths_squared()
	m = stars.mass / stars.mass.sum()
	d2m = zip(d2, m)
	d2m.sort(key=operator.itemgetter(0))

	iL = 0
	Lagrad = []

	while len(Lagrad) < len(massf):

		mt = 0 | units.none

		for d2i, mi in d2m:

			if mt < massf[iL]:
				mt += mi
			
			if mt >= massf[iL]:
		    		Lagrad.append(d2i.sqrt())
				break
		    		
		iL += 1

	return Lagrad

def main(N, seed_cluster, seed_salpeter, mmin, mmax, t_end, dt, r_cluster, r_gal, m_gal, alpha, x0, vyfrac):
	
	'''
	initialization of cluster, galaxy
	'''

	seed = None #initialize the random seed, will change value

	#set up the background galaxy
	galaxy_code = GalacticCenterGravityCode(r_gal|units.parsec, m_gal|units.MSun, alpha)

	#set up Salpeter mass distribution for cluster
	seed = seed_salpeter
	masses = new_salpeter_mass_distribution(N, mass_min=mmin|units.MSun, mass_max=mmax|units.MSun)

	#initialize coordinates of stars in the cluster
	converter_cluster = nbody_system.nbody_to_si(masses.sum(), r_cluster|units.parsec)

	seed = seed_cluster
	bodies = new_plummer_model(N, convert_nbody=converter_cluster)

	#endow each star in the Plummer sphere with the appropriate mass
	for i, body in enumerate(bodies):
		body.mass = masses[i]

	#move cluster to initial phase space coordinates:
	for body in bodies:
		body.x += x0|units.parsec
		body.vy += vyfrac*galaxy_code.circular_velocity(x0|units.parsec)
	
	'''
	set up appropriate solver (gravity + stellar evolution)
	'''

	#direct N-body for the cluster
	cluster_code = ph4(converter_cluster)
	cluster_code.particles.add_particles(bodies)

	#stellar evolution code to capture mass loss
	stellar_code = SeBa()
	stellar_code.particles.add_particles(bodies)

	#bridge the cluster with the background galaxy, not the other way around
	gravity_code = bridge.Bridge()
	gravity_code.add_system(cluster_code, (galaxy_code,))

	sim_times_unitless = np.arange(0., t_end, dt)
	sim_times = [ t|units.Myr for t in sim_times_unitless ]

	i = 0 #for naming of .png files for visual

	LR_one, LR_two, LR_three, LR_four = [], [], [], []
	system_mass = []

	t0 = time.time()

	for i, t in enumerate(sim_times):

		if i%10 == 0:
			print('sim_time is %.02f Myr'%(t.value_in(units.Myr)))
			print('clock time is %.02f s'%(time.time()-t0))
			print('')

		MassFraction = [0.1, 0.25, 0.5, 0.75]

		'''
		evolve stellar and gravity separately, see textbook
		'''

		LR = LagrangianRadii(total_system.particles, massf=MassFraction)

		LR_one.append(LR[0].value_in(units.parsec))
		LR_two.append(LR[1].value_in(units.parsec))
		LR_three.append(LR[2].value_in(units.parsec))
		LR_four.append(LR[3].value_in(units.parsec))

		sys_mass = sum([ particle.mass.value_in(units.MSun) for particle in bodies ])

		system_mass.append(sys_mass)

		x = bodies.x.value_in(units.parsec)
		y = bodies.y.value_in(units.parsec)

		plt.figure()
		plt.scatter(x,y,s=1)
		plt.xlim(-200, 200)
		plt.ylim(-200, 200)
		plt.title('Cluster, time = %.02f Myr'%(t.value_in(units.Myr)))
		plt.xlabel(r'$x$ (pc)', fontsize=16)
		plt.ylabel(r'$y$ (pc)', fontsize=16)
		plt.tight_layout()
		plt.savefig('frame_%s.png'%(str(i).rjust(4, '0')))
		plt.close()

	#final snapshot of the simulation
	x = bodies.x.value_in(units.parsec)
	y = bodies.particles.y.value_in(units.parsec)

	plt.figure()
	plt.scatter(x,y,s=1)
	plt.xlim(-200, 200)
	plt.ylim(-200, 200)
	plt.title('Final Snapshot'%(t.value_in(units.Myr)))
	plt.xlabel(r'$x$ (pc)', fontsize=16)
	plt.ylabel(r'$y$ (pc)', fontsize=16)
	plt.tight_layout()
	plt.savefig('final_snapshot.png')
	plt.close()

	#Lagrangian radius plot

	plt.figure()
	plt.semilogy(sim_times_unitless, LR_one, label = r'$10 \%$')
	plt.semilogy(sim_times_unitless, LR_two, label = r'$25 \%$')
	plt.semilogy(sim_times_unitless, LR_three, label = r'$50 \%$')
	plt.semilogy(sim_times_unitless, LR_four, label = r'$75 \%$')
	plt.title('Lagrangian Radii')
	plt.xlabel(r'$t$ (Myr)', fontsize=16)
	plt.ylabel(r'Lagrangian Radius (pc)', fontsize=16)
	plt.legend(loc='best')
	plt.tight_layout()
	plt.savefig('lagrangian_radii.png')
	plt.close()

	#Cluster mass plot

	plt.figure()
	plt.plot(sim_times_unitless, system_mass)
	plt.title('Cluster Mass')
	plt.xlabel(r'$t$ (Myr)', fontsize=16)
	plt.ylabel(r'Mass ($M_{\odot}$)', fontsize=16)
	plt.legend(loc='best')
	plt.tight_layout()
	plt.savefig('lagrangian_radii.png')
	plt.close()

	cluster_code.stop()		

	print('testing has worked!')
	return 1

def new_option_parser():
	
	'''
	define an option parser for initial conditions
	'''

	optparser = OptionParser()
	optparser.add_option('--N', dest='N', type='int', default=20)
	optparser.add_option('--seed_cluster', dest='seed_cluster', type='int', default=121)
	optparser.add_option('--seed_salpeter', dest='seed_salpeter', type='int', default=242)
	optparser.add_option('--mmin', dest='mmin', type='float', default=1.)
	optparser.add_option('--mmax', dest='mmax', type='float', default=100.)
	optparser.add_option('--t_end', dest='t_end', type='float', default=10.)
	optparser.add_option('--dt', dest='dt', type='float', default=0.1)		
	optparser.add_option('--r_cluster', dest='r_cluster', type='float', default=0.1)
	optparser.add_option('--r_gal', dest='r_gal', type='float', default=100.)
	optparser.add_option('--m_gal', dest='m_gal', type='float', default=7.22e9)
	
	#the textbook uses alpha=1.2 in their galacticcentergravitycode but does not indicate that this is a Plummer model
	optparser.add_option('--alpha', dest='alpha', type='float', default=1.0)
	optparser.add_option('--x0', dest='x0', type='float', default=200.)
	optparser.add_option('--vyfrac', dest='vyfrac', type='float', default=0.2)

	return optparser

if __name__ in '__main__':
	o, arguments = new_option_parser().parse_args()
	main(**o.__dict__)

