import matplotlib
matplotlib.use('Agg')
from make_initial_conditions import planetary_system
from run_debris_disk_around_moon import integrate_solar_system
from amuse.lab import units
from matplotlib import pyplot

def plot_track():

    #figure = pyplot.figure(figsize=(10, 10))
    #pyplot.rcParams.update({'font.size': 10})
    #plot = figure.add_subplot(1,1,1)
    #ax = pyplot.gca()
    #ax.minorticks_on()
    #ax.locator_params(nbins=3)

    t = integrate_solar_system()[0]
    sys_e = integrate_solar_system()[1]
    disc_e = integrate_solar_system()[2]
    bridge_e = integrate_solar_system()[3]

    print 'energy imported'

    x_planet = integrate_solar_system()[4]
    y_planet = integrate_solar_system()[5]
    x_moon = integrate_solar_system()[6]
    y_moon = integrate_solar_system()[7]
    x_disc = integrate_solar_system()[8]
    y_disc = integrate_solar_system()[9]

    print 'positions imported'

    pyplot.figure()
    pyplot.xlabel('x [au]')
    pyplot.ylabel('y [au]')
    pyplot.scatter(x_disc[0].value_in(units.AU), y_disc[0].value_in(units.AU), color='green',s=1, label='disk')
    pyplot.scatter(x_moon[0].value_in(units.AU), y_moon[0].value_in(units.AU), color='red', label='moon')
    pyplot.scatter(x_planet[0].value_in(units.AU), y_planet[0].value_in(units.AU), color='blue', label='planet')
    pyplot.legend()
    pyplot.savefig('/home/rywang/capexam/plots/initial_pos.png')

    pyplot.figure()
    pyplot.xlabel('x [au]')
    pyplot.ylabel('y [au]')
    pyplot.scatter(x_disc[49].value_in(units.AU), y_disc[49].value_in(units.AU), color='green',s=1, label='disk')
    pyplot.scatter(x_moon[49].value_in(units.AU), y_moon[49].value_in(units.AU), color='red', label='moon')
    pyplot.scatter(x_planet[49].value_in(units.AU), y_planet[49].value_in(units.AU), color='blue', label='planet')
    pyplot.legend()
    pyplot.savefig('/home/rywang/capexam/plots/equal_distance_1_pos.png')

    pyplot.figure()
    pyplot.xlabel('x [au]')
    pyplot.ylabel('y [au]')
    pyplot.scatter(x_disc[99].value_in(units.AU), y_disc[99].value_in(units.AU), color='green',s=1, label='disk')
    pyplot.scatter(x_moon[99].value_in(units.AU), y_moon[99].value_in(units.AU), color='red', label='moon')
    pyplot.scatter(x_planet[99].value_in(units.AU), y_planet[99].value_in(units.AU), color='blue', label='planet')
    pyplot.legend()
    pyplot.savefig('/home/rywang/capexam/plots/equal_distance_2_pos.png')

    pyplot.figure()
    pyplot.xlabel('x [au]')
    pyplot.ylabel('y [au]')
    pyplot.scatter(x_disc[-1].value_in(units.AU), y_disc[-1].value_in(units.AU), color='green',s=1, label='disk')
    pyplot.scatter(x_moon[-1].value_in(units.AU), y_moon[-1].value_in(units.AU), color='red', label='moon')
    pyplot.scatter(x_planet[-1].value_in(units.AU), y_planet[-1].value_in(units.AU), color='blue', label='planet')
    pyplot.legend()
    pyplot.savefig('/home/rywang/capexam/plots/final_pos.png')

    pyplot.xlabel('time')
    pyplot.ylabel('energy')
    #pyplot.yscale('log')
    pyplot.plot(t.number, sys_e.number, label='planets energy')
    pyplot.plot(t.number, disc_e.number, label='disk energy')
    pyplot.plot(t.number, bridge_e.number, label='bridge energy')
    pyplot.legend()
    pyplot.savefig('/home/rywang/capexam/plots/energy_time.png')

    #pyplot.show()

if __name__ in ('__main__'):

    plot_track()
