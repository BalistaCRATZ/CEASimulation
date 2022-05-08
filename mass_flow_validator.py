import json
import math

import numpy as np
import scipy.optimize
from octopus import Fluid, Manifold, Orifice, PropertySource, utils
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend

# Plotting libraries
import matplotlib.pyplot as plt

#Some variables
initial_ullage_pressure = 19e5
initial_ullage_fraction = 0.05 # fraction of total volume replaced by ullage
initial_propellant_mass = 150
fuel_density = 841.9 #Density of IPA and H20 mixture at 15C, from ic.gc.ca and NIST
oxidiser_density = 1001.2 #Density of N20 at -20C, from NIST
initial_rho = ((fuel_density + 3.5 * oxidiser_density)/4.5)*(1-initial_ullage_fraction)
tank_volume = initial_propellant_mass/initial_rho
dt = 0.1
t_max = 10

def get_density(ullage):
    return ((fuel_density + 3.5 * oxidiser_density)/4.5)*(1-ullage)

## Code for the injector pressure drop calculator

data = json.load(open("variables.json"))

Isoblend = newFuelBlend(fuelL = ["Isopropanol", "H2O"], fuelPcentL = [80, 20])
C = CEA_Obj(oxName = "N2O", fuelName= Isoblend, pressure_units= 'Pa', specific_heat_units='kJ/kg-K', density_units = 'kg/m^3', temperature_units = 'K')

# header function for the solver (can copy this format for other solver functions)
def annular_gap_func(mdot, dp):
    Ao = data["area_reduction_ratio"] * 0.25 * np.pi * (data["pintle_gap_OD"] ** 2 - data["pintle_sleeve_OD"] ** 2)
    return dp - utils.dp_annular_gap(
        D_outer=data["pintle_gap_OD"],
        D_inner=data["pintle_sleeve_OD"],
        mdot=mdot,
        L=data["pintle_gap_length"],
        rho=data["fuel_density"],
        mu=data["fuel_viscosity"]
    ) - (0.5 * mdot ** 2) / (data["fuel_density"] * Ao ** 2)


# wrapped functions for mass flow rate as a function of upstream and downstream pressure
def mdot_fuel_calc(p_tank, p_chamber):
    return scipy.optimize.fsolve(func=annular_gap_func, x0=np.array(0.167), args=(p_tank - p_chamber))[0]  # 0.178


def mdot_ox_calc(p_tank, p_chamber, orifice: Orifice):
    orifice.manifold.parent._p = p_tank
    return orifice.m_dot_dyer(p_chamber)


def mdot_nozzle_calc(p_chamber, OF):
    CHAMBER_TEMP = C.get_Temperatures(p_chamber, OF)[1] #expansion ratio is set to default of 40
    (MOLECULAR_MASS, gamma) = C.get_Throat_MolWt_gamma(p_chamber, OF)
    cp = C.get_HeatCapacities(p_chamber, OF)[1] #in kJ/kg-K
    R = cp - (cp/gamma)
    R = 1000*R

    return ((data["throat_area"] * p_chamber) / np.sqrt(CHAMBER_TEMP)) * np.sqrt(gamma / R) * math.pow(
        (gamma + 1) / 2, -(gamma + 1) / (2 * (gamma - 1)))

def mass_flow_diff(p_chamber, p_tank, orifice: Orifice):
    m_f = mdot_fuel_calc(p_tank, p_chamber)
    m_o = mdot_ox_calc(p_tank, p_chamber, orifice)
    OF = m_o / m_f
    m_n = mdot_nozzle_calc(p_chamber, OF)
    return m_n - m_f - m_o


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

    p_tank_max = 25e5
    p_tank_list = np.arange(10e5, p_tank_max, 0.25e5)
    p_chamber_list = []
    mdot_ox_list = []
    mdot_fuel_list = []

    for p_tank in p_tank_list:
        p_chamber, OF = get_chamber_pressure_and_OF(p_tank, ox_area)
        p_chamber_list.append(p_chamber)
        mdot_fuel = mdot_fuel_calc(p_tank, p_chamber)
        oxidiser = Fluid(name="NitrousOxide", eos="HEOS")
        ox_property_source = PropertySource(p=p_tank, T=data["tank_temp"])
        oxidiser.set_state(P=p_tank, T=data["tank_temp"])
        ox_manifold = Manifold(fluid=oxidiser, parent=ox_property_source, A=1)
        ox_orifice = Orifice(manifold=ox_manifold, A=ox_area, Cd=0.7)
        mdot_ox = mdot_ox_calc(p_tank, p_chamber, ox_orifice)
        mdot_ox_list.append(mdot_ox)
        mdot_fuel_list.append(mdot_fuel)
        # mdot = mdot_nozzle_calc(p_chamber, OF)
        # prev_density = density
        # density = get_density(ullage)
        # Vdot = mdot/density
        # propellant_volume = (1-ullage)*tank_volume - Vdot*dt
        # ullage = (tank_volume - propellant_volume)/tank_volume


        ## Modelling Isentropic expansion of gas to determine new tank_pressure
        # (MOLECULAR_MASS, gamma) = C.get_Throat_MolWt_gamma(p_chamber, OF)
        # cp = C.get_HeatCapacities(p_chamber, OF)[1] #in kJ/kg-K
        # R = cp - (cp/gamma) # in KJ/kg-K
        # R = 1000*R  # in J/kg-K
        # p_tank = ((density/prev_density)**(gamma))*p_tank #Using isentropic formula

    # for i in range(len(t_list)):
    #     print(f"Time : {t_list[i]}")
    #     print(f"Tank Pressure: {p_tank_list[i]}")
    #     print(f"Chamber Pressure: {p_chamber_list[i]}")
    # print(OF)

    plt.xlabel("Tank Pressure(Pa)")
    plt.ylabel("Mass Flow Rate (kg s-1)")
    plt.plot (p_tank_list, mdot_fuel_list, label = "Fuel Mass Flow Rate (kg s-1)")
    plt.plot(p_tank_list, mdot_ox_list, label = "Oxidiser Mass Flow Rate (kg s-1)")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()