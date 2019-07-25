from make_initial_conditions import planetary_system
from amuse.lab import units
from amuse.couple import bridge
from amuse.ext.composition_methods import *
from amuse.lab import nbody_system
from amuse.lab import Huayno
from amuse.units import quantities
from matplotlib import pyplot
import subprocess

def integrate_solar_system():

    # Get initial conditions.
    solar = planetary_system()[0]
    debris = planetary_system()[1]

    star = solar[0]
    planet = solar[1]
    moon = solar [2]

    # Integrators set to Huayno.
    converter = nbody_system.nbody_to_si(star.mass, 1. | units.AU)
    system_gravity = Huayno(converter)
    system_gravity.particles.add_particle(star)
    system_gravity.particles.add_particle(planet)
    system_gravity.particles.add_particle(moon)

    disk_gravity = Huayno(converter)
    disk_gravity.particles.add_particle(debris)

    # Set up channels.
    channel_from_system_to_framework = system_gravity.particles.new_channel_to(solar)
    channel_from_disk_to_framework = disk_gravity.particles.new_channel_to(debris)

    # Set up bridge.
    gravity = bridge.Bridge(use_threading=False, method=SPLIT_4TH_S_M4)
    gravity.add_system(system_gravity, (disk_gravity,))
    gravity.add_system(disk_gravity, (system_gravity,))

    # Record positions.
    x_planet = quantities.AdaptingVectorQuantity()
    y_planet = quantities.AdaptingVectorQuantity()
    x_moon = quantities.AdaptingVectorQuantity()
    y_moon = quantities.AdaptingVectorQuantity()
    x_debris = quantities.AdaptingVectorQuantity()
    y_debris = quantities.AdaptingVectorQuantity()

    # Evolving time.
    t = 0 | units.day
    dt = 1 | units.day # to make sure getting 250 snapshots
    t_end = 1 | units.yr # 29.4571 yr Saturn's orbital period

    # Record time and energies.
    time = quantities.AdaptingVectorQuantity()
    system_E = quantities.AdaptingVectorQuantity()
    disk_E = quantities.AdaptingVectorQuantity()
    bridge_E = quantities.AdaptingVectorQuantity()


    # Evolve model.
    print 'Begin evolving'
    while t < t_end:
        t += dt
        gravity.evolve_model(t)

        channel_from_system_to_framework.copy()
        channel_from_disk_to_framework.copy()

        x_planet.append(planet.x)
        y_planet.append(planet.y)
        x_moon.append(moon.x)
        y_moon.append(moon.y)
        x_debris.append(debris.x)
        y_debris.append(debris.y)

        time.append(t)
        system_E.append(system_gravity.kinetic_energy + system_gravity.potential_energy)
        disk_E.append( disk_gravity.kinetic_energy + disk_gravity.potential_energy)
        bridge_E.append(gravity.kinetic_energy + gravity.potential_energy)

        pyplot.clf()
        pyplot.xlabel('x [au]')
        pyplot.ylabel('y [au]')
        #pyplot.xlim(-2.2, -1.9)
        #pyplot.ylim(8.7, 8.8)
        pyplot.scatter(star.x.value_in(units.AU), star.y.value_in(units.AU), color='black', label='star')
        pyplot.scatter(planet.x.value_in(units.AU), planet.y.value_in(units.AU), color='blue', label='planet')
        pyplot.scatter(moon.x.value_in(units.AU), moon.y.value_in(units.AU), color='red', label='moon')
        pyplot.scatter(debris.x.value_in(units.AU), debris.y.value_in(units.AU), color='green',s=1, label='disk')
        #pyplot.legend()
        pyplot.savefig('/home/rywang/capexam/plots/animation/positions_{0:03}.png'.format(t.number))

        print 'Generating plot{}'.format(t.number)

    gravity.stop()
    print 'End evolving'

    # Making animation.
    subprocess.call("mencoder 'mf:///home/rywang/capexam/plots/animation/positions_*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o /home/rywang/capexam/plots/animation/animation.mpg", shell=True)

    return time, system_E, disk_E, bridge_E, x_planet, y_planet, x_moon, y_moon, x_debris, y_debris

#if __name__ in ('__main__'):

#    integrate_solar_system()
