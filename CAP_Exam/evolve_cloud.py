#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:40:11 2019

@author: rywang
"""

import numpy as np
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot

import math
from amuse.lab import *
from optparse import OptionParser
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.ic.gasplummer import new_plummer_gas_model



class GalacticCenterGravityCode(object):
    def __init__(self,R, M, alpha):
        self.radius=R
        self.mass=M
        self.alpha=alpha

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.mass*(r/self.radius)**self.alpha  
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def circular_velocity(self,r):  
        m=self.mass*(r/self.radius)**self.alpha  
        vc=(constants.G*m/r)**0.5
        return vc


    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2)**0.5
        c=constant.G*self.mass/self.radius**self.alpha    
        phi=c/(alpha-1)*(r**(self.alpha-1)-R**(self.alpha-1))
        return phi    
        
def plummer_model_A(nbody, N, Mcloud,
                            Rcloud, parameters = []):
    np.random.seed(353)

    converter=nbody_system.nbody_to_si(Mcloud,Rcloud)
    bodies=new_plummer_gas_model(N,convert_nbody=converter)

    code=nbody(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

def plummer_model_B(nbody, N, Mcloud,
                            Rcloud, parameters = []):
    np.random.seed(424)

    converter=nbody_system.nbody_to_si(Mcloud,Rcloud)
    bodies=new_plummer_gas_model(N,convert_nbody=converter)

    code=nbody(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

def merge_two_clouds(bodies, particles_in_encounter, tcoll):
    common_position = particles_in_encounter.center_of_mass()
    common_velosity = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = common_position
    new_particle.velocity = common_velosity
    new_particle.radius = \
        particles_in_encounter.radius.sum()/len(particles_in_encounter.radius)
    bodies.add_particles(new_particle)
    bodies.remove_particles(particles_in_encounter)
    return new_particle


def plot_cloud(x, y):

    from prepare_figure import single_frame, get_distinct
    colors = get_distinct(1)     
    frame = single_frame('X [pc]', 'Y [pc]')
    
    pyplot.title('Cloud particle positions')
    pyplot.xlim(-200, 200)
    pyplot.ylim(-200, 200)
    pyplot.scatter(x,y, c=colors[0], s=50, lw=0)

    save_file = 'positions_M1e5.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file, '\n'
    pyplot.show()

def evolve_cloud_in_galaxy(N, RinitA, RinitB, tend, timestep, MA, RA, MB, RB):


    Rgal = 1 | units.kpc
    Mgal = 7.22e9 | units.MSun
    alpha = 1.2
    galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)

    cloud_codeA = plummer_model_A(BHTree, N, MA, RA,
                                           parameters=[("epsilon_squared",
                                                        (0.01 | units.parsec)**2)])
    
    cloud_codeB = plummer_model_B(BHTree, N, MB, RB,
                                           parameters=[("epsilon_squared",
                                                        (0.01 | units.parsec)**2)])
    
    cloudA = cloud_codeA.particles.copy()    
    cloudB = cloud_codeA.particles.copy() 
    cloudA.x += RinitA
    cloudB.x += RinitB
    cloudA.vy = 0.01*galaxy_code.circular_velocity(RinitA)
    cloudB.vy = 0.05*galaxy_code.circular_velocity(RinitB)
    channelA = cloudA.new_channel_to(cloud_codeA.particles)
    channelA.copy_attributes(["x","y","z","vx","vy","vz"])
    channelB = cloudB.new_channel_to(cloud_codeB.particles)
    channelB.copy_attributes(["x","y","z","vx","vy","vz"])

    system = bridge(verbose=False)
    system.add_system(cloud_codeA, (galaxy_code,))
    system.add_system(cloud_codeB, (galaxy_code,))

    times = quantities.arange(0|units.Myr, tend, timestep)
    for i,t in enumerate(times):
        system.evolve_model(t,timestep=timestep)
          
    x = system.particles.x.value_in(units.parsec)
    y = system.particles.y.value_in(units.parsec)
    cloud_codeA.stop()
    cloud_codeB.stop()

    return x, y

if __name__ == "__main__":
    N=100
    RinitA=200 | units.parsec
    RinitB=194 | units.parsec
    timestep=0.1 | units.Myr
    endtime = 3 | units.Myr
    McloudA = 1.e5 | units.MSun
    RcloudA = 0.1 | units.parsec
    McloudB = 1.e5 | units.MSun
    RcloudB = 1.0 | units.parsec
    x, y = evolve_cloud_in_galaxy(N, RinitA, RinitB, endtime, timestep,
                                    McloudA, RcloudA, McloudB, RcloudB)
    
    t_end = 30
    nsteps = 10
    Mmax=100 
    Qvir=0.5
    
    t_end = t_end | nbody_system.time
    dt = t_end/float(nsteps)

    bodies = new_plummer_model(N)
    masses = new_powerlaw_mass_distribution(N, mass_min=0.1|units.MSun,
                                            mass_max=10|units.MSun, alpha=-2.35)
    masses /= masses.sum()
    bodies.mass = masses | nbody_system.mass

    bodies.radius = 0.001 | nbody_system.length
    Mtot_init = 1 | nbody_system.mass

    gravity = ph4(number_of_workers=4)
    gravity.particles.add_particles(bodies)
    bodies.scale_to_standard(virial_ratio=Qvir)

    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    Nenc = 0
    dEk_enc = zero    
    dEp_enc = zero

    t = []
    MassFraction = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    r10 = []
    r25 = []
    r50 = []
    r75 = []
    tcoll = []
    rcoll = []
    mcoll = []
    
    while time < t_end:

        RL = LagrangianRadii(gravity.particles, massf=MassFraction)
        pos,coreradius,coredens = \
            gravity.particles.densitycentre_coreradius_coredens()
        print "Cluster at time=", time, "core_radius=", coreradius, \
              "L_radii =",
        for rl in RL:
            print rl.number, " ",
        print "length"

        t.append(time.number)
        r10.append(RL[4].number)
        r25.append(RL[5].number)
        r50.append(RL[6].number)
        r75.append(RL[7].number)

        gravity.evolve_model(time+dt)
        time = gravity.model_time

        if stopping_condition.is_set():
            Ek_enc = gravity.kinetic_energy 
            Ep_enc = gravity.potential_energy
            print "At time=", time, "number of encounters=", \
                len(stopping_condition.particles(0))
            for ci in range(len(stopping_condition.particles(0))): 
                particles_in_encounter = Particles(
                    particles=[stopping_condition.particles(0)[ci],
                               stopping_condition.particles(1)[ci]])
                particles_in_encounter = \
                    particles_in_encounter.get_intersecting_subset_in(bodies)

                new_particle = merge_two_clouds(bodies,
                                               particles_in_encounter,
                                               gravity.model_time)
                bodies.synchronize_to(gravity.particles)
                Nenc+=1
                RL = LagrangianRadii(gravity.particles, massf=MassFraction)

                print "Resolve encounter number", Nenc
                pos,coreradius,coredens = \
                    gravity.particles.densitycentre_coreradius_coredens()
                print "Collision at time=", time, new_particle.mass.sum(), \
                    new_particle.position.length(), "Nstars= ", len(bodies), \
                    "Ncoll=", Nenc, "core_radius=", coreradius, "L_radii=", 
                for rl in RL:
                    print rl.number, " ",
                print "length"

                pos = new_particle[0].position.number
                rc = math.sqrt(np.sum(pos**2))
                mc = new_particle[0].mass.number
                tcoll.append(time.number)
                rcoll.append(rc)
                mcoll.append(mc)

            dEk_enc += Ek_enc - gravity.kinetic_energy 
            dEp_enc += Ep_enc - gravity.potential_energy

        bodies.move_to_center()
        channel_from_gd_to_framework.copy()

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E=", Etot, "Q=", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
        print "dE(enc)=", dEk_enc, dEp_enc
        Etot_init -= dE
        Etot_prev = Etot

    gravity.stop()
    
    plot_cloud(x, y)

    pyplot.figure(figsize=(14,12))
    pyplot.title('Lagrangian radii')
    pyplot.plot(t, r10, 'r', label='10%')
    pyplot.plot(t, r25, 'y', label='25%')
    pyplot.plot(t, r50, 'b', label='50%')
    pyplot.plot(t, r75, 'g', label='75%')
    pyplot.legend()
    size = np.array(mcoll)
    mmin = np.amin(size)
    size = 2 + 15*np.log10(size/mmin)
    pyplot.xlabel('t [N-body]')
    pyplot.ylabel('R [N-body]')
    pyplot.yscale('log')

    save_file = 'gravity_collision.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
    
