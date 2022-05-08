import json
import math
import numpy as np
import scipy.optimize
from octopus import Fluid, Manifold, Orifice, PropertySource, utils
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend

data = json.load(open("variables.json"))

Isoblend = newFuelBlend(fuelL = ["Isopropanol", "H2O"], fuelPcentL = [80, 20])
C = CEA_Obj(oxName = "N2O", fuelName= Isoblend, pressure_units= 'Pa', specific_heat_units='kJ/kg-K', density_units = 'kg/m^3', temperature_units = 'K')

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
    CHAMBER_TEMP = C.get_Temperatures(p_chamber, OF)[1]  # expansion ratio is set to default of 40
    (MOLECULAR_MASS, gamma) = C.get_Throat_MolWt_gamma(p_chamber, OF)
    cp = C.get_HeatCapacities(p_chamber, OF)[1]  # in kJ/kg-K
    R = cp - (cp / gamma)
    R = 1000 * R

    return ((data["throat_area"] * p_chamber) / np.sqrt(CHAMBER_TEMP)) * np.sqrt(gamma / R) * math.pow(
        (gamma + 1) / 2, -(gamma + 1) / (2 * (gamma - 1)))


def mass_flow_diff(p_chamber, p_tank, orifice: Orifice):
    m_f = mdot_fuel_calc(p_tank, p_chamber)
    m_o = mdot_ox_calc(p_tank, p_chamber, orifice)
    OF = m_o / m_f
    m_n = mdot_nozzle_calc(p_chamber, OF)
    return m_n - m_f - m_o
