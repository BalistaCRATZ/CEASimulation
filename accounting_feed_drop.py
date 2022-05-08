
import numpy as np
import scipy.optimize
from octopus import Fluid, Manifold, Orifice, PropertySource, utils
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend

# Plotting libraries
import matplotlib.pyplot as plt

# our utilities
from utils import *

#Some variables
initial_ullage_pressure = 19e5
initial_ullage_fraction = 0.05 # fraction of total volume replaced by ullage
initial_propellant_mass = 15
fuel_density = 841.9 #Density of IPA and H20 mixture at 15C, from ic.gc.ca and NIST
oxidiser_density = 1001.2 #Density of N20 at -20C, from NIST
initial_rho = ((fuel_density + 3.5 * oxidiser_density)/4.5)*(1-initial_ullage_fraction)
tank_volume = initial_propellant_mass/initial_rho
print(tank_volume)
dt = 0.01
t_max = 0.8
feed_pressure_drop = 3e5

def get_density(ullage):
    return ((fuel_density + 3.5 * oxidiser_density)/4.5)*(1-ullage)

## Code for the injector pressure drop calculator

data = json.load(open("variables.json"))


def get_chamber_pressure_and_OF(p_tank, ox_area):
    # setup octopus sim for oxidiser
    oxidiser = Fluid(name="NitrousOxide", eos="HEOS")
    ox_property_source = PropertySource(p=p_tank, T=data["tank_temp"])
    oxidiser.set_state(P=p_tank, T=data["tank_temp"])
    ox_manifold = Manifold(fluid=oxidiser, parent=ox_property_source, A=1)
    ox_orifice = Orifice(manifold=ox_manifold, A=ox_area, Cd=0.7)
    # END INIT #
    chamber_pressure = scipy.optimize.fsolve(func=mass_flow_diff, x0=np.array(10e5), args=(p_tank, ox_orifice))[0]
    OF = mdot_ox_calc(p_tank, chamber_pressure, ox_orifice)/mdot_fuel_calc(p_tank, chamber_pressure)
    return [chamber_pressure, OF]

def main():
    # DO NOT EDIT -- MUST INCLUDE THIS INIT #
    # precalculating fuel flow area from datafile
    ox_area = data["slot_width"] * data["num_slots"] * (
            data["slot_height"] - data["pintle_sleeve_thickness"] * np.sin(data["alpha"])
    )
    t_list = np.arange(0,t_max, dt)
    p_tank = initial_ullage_pressure
    ullage = initial_ullage_fraction
    density = initial_rho

    p_tank_list = []
    p_chamber_list = []
    OF_ratio_list = []
    mdot_list = []

    for t in t_list:
        p_tank_list.append(p_tank)
        p_feed = p_tank - feed_pressure_drop
        p_chamber, OF = get_chamber_pressure_and_OF(p_feed, ox_area)
        OF_ratio_list.append(OF)
        p_chamber_list.append(p_chamber)
        mdot = mdot_nozzle_calc(p_chamber, OF)
        mdot_list.append(mdot)
        prev_density = density
        density = get_density(ullage)
        Vdot = mdot/density
        propellant_volume = (1-ullage)*tank_volume - Vdot*dt
        ullage = (tank_volume - propellant_volume)/tank_volume


        ## Modelling Isentropic expansion of gas to determine new tank_pressure
        (MOLECULAR_MASS, gamma) = C.get_Throat_MolWt_gamma(p_chamber, OF)
        cp = C.get_HeatCapacities(p_chamber, OF)[1] #in kJ/kg-K
        R = cp - (cp/gamma) # in KJ/kg-K
        R = 1000*R  # in J/kg-K
        p_tank = ((density/prev_density)**(gamma))*p_tank #Using isentropic formula

    for i in range(len(t_list)):
        print(f"Time : {t_list[i]}")
        print(f"Tank Pressure: {p_tank_list[i]}")
        print(f"Chamber Pressure: {p_chamber_list[i]}")
    print(OF)

    plt.xlabel("Time(s)")
    plt.ylabel("Pressure (Pa)")
    plt.plot (t_list, p_tank_list, label = "Tank Pressure (Pa)")
    plt.plot(t_list, p_chamber_list, label = "Chamber Pressure (Pa)")
    # plt.plot(t_list, OF_ratio_list, label = "OF Ratio")
    plt.legend()
    plt.show()
    plt.savefig("figs/tfig")

    # plt.xlabel("Time(s)")
    # plt.ylabel("OF Ratio")
    # plt.plot(t_list, OF_ratio_list)
    # plt.show()

    # plt.xlabel("Time(s)")
    # plt.ylabel("mdot (kgs)")
    # plt.plot(t_list, mdot_list)
    # plt.show()


if __name__ == "__main__":
    main()