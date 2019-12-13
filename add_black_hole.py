from amuse.lab import *

def BH_maker(BH_mass, x0_BH, y0_BH, vx0_BH):

	BH = Particles(1)
	BH.mass = BH_mass|units.MSun
	BH.radius = 2*units.constants.G*BH.mass / (units.constants.c)**2
	BH.x = x0_BH|units.AU
	BH.y = y0_BH|units.AU
	BH.z = 0.|units.AU
	BH.vx = vx0_BH|(units.km/units.s)
	BH.vy = 0.|(units.km/units.s)
	BH.vz = 0.|(units.km/units.s)

	return BH
