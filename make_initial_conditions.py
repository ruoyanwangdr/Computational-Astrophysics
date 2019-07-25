import matplotlib
matplotlib.use('Agg')
import numpy
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.lab import Particles, units
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.lab import nbody_system
from amuse.lab import new_powerlaw_mass_distribution
from matplotlib import pyplot

def planetary_system():

    # Generate the sun, Saturn, and Titan. Positions and velocities acquired from textbook.

    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (0.005717, -0.00538, -0.0000213) | units.AU
    sun.velocity = (7.893, 11.894, 0.20642) | (units.m/units.s)

    saturn = particles[1]
    saturn.mass = 95.16 | units.MEarth
    saturn.radius = 763.59 | units.REarth
    saturn.position = (-2.075, 8.7812, 0.3273) | units.AU
    saturn.velocity = (-9.9109, 2.236, -0.2398) | units.kms

    titan = particles[2]
    titan.mass = 0.0225 | units.MEarth
    titan.radius = 0.404 | units.REarth
    titan.position = (-2.0668, 8.7264, 0.3355) | units.AU # Saturn's postion plus 1221870 km in AU.
    titan.velocity = (-4.3409, 3.334, 5.3302) | units.kms # Saturn's velocity plus 5.57.

    particles.move_to_center()

    numpy.random.seed(1)

    # Generate the disk.

    disk_a = 1221870 | units.km # titan's semi-major axis
    disk_e = 0.0288 # eccentricity
    hill_radius = disk_a.value_in(units.AU) * (1 - disk_e) * (titan.mass / (3 * saturn.mass))**(1/3)
    disk_rmin = 0.03 * hill_radius
    disk_rmax = 0.5 * hill_radius

    disk_mass = new_powerlaw_mass_distribution(1000,
                                            2.2e-5*titan.mass,
                                            2.2e-3*titan.mass,
                                            alpha=-2.0)
    Mtot = disk_mass.sum()
    converter = nbody_system.nbody_to_si(titan.mass, 1. | units.AU)

    disk = ProtoPlanetaryDisk(100,
                            convert_nbody=converter,
                            densitypower=1.5,
                            Rmin=disk_rmin,
                            Rmax=disk_rmax,
                            q_out=10.,
                            discfraction=0.1).result

    disk.move_to_center()

    disk.x += titan.x
    disk.y += titan.y
    disk.vx += titan.vx
    disk.vy += titan.vy


    planet_attributes=["x", "y", "z", "vx", "vy", "vz",
                       "mass", "semimajor_axis", "eccentricity"]
    disk_attributes=["x", "y", "z", "vx", "vy", "vz",
                     "mass", "u", "rho", "h_smooth"]

    write_set_to_file(particles, 'SunSaturnTitan_i0000.amuse', 'amuse', attribute_names=planet_attributes)
    write_set_to_file(disk, 'disk_i0000.amuse', 'amuse', attribute_names=disk_attributes)

    #test = read_set_from_file('SunSaturnTitan_i0000.amuse', 'amuse')
    #print test


    return particles, disk

#if __name__ in ('__main__'):
#    planetary_system()
