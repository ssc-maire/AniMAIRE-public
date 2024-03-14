import numpy as np
import pandas as pd
import pkg_resources
import ParticleRigidityCalculationTools as PRCT
from numba import njit
from metpy.constants import g
from metpy.calc import height_to_pressure_std
from metpy.units import units


NM64_proton_response = pd.read_csv(pkg_resources.resource_filename(__name__,'NM64_proton_response.csv'), 
                                    header=None)
NM64_proton_response.columns = ["proton_rigidity_GV","cts_per_primary_per_cm2"]

NM64_proton_filepath_tabulated_2013 = pkg_resources.resource_filename(__name__,'NM64_proton_response_tabulated_2013.csv')
NM64_proton_response_tabulated_2013_epn = pd.read_csv(NM64_proton_filepath_tabulated_2013,header=None)
NM64_proton_response_tabulated_2013_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_proton_response_tabulated_2013_epn["Yield_per_m2_per_sr"] = NM64_proton_response_tabulated_2013_epn["Yield_arb_units"] * 0.109
NM64_proton_response_tabulated_2013 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_proton_response_tabulated_2013_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                             fluxInEnergyMeVform=NM64_proton_response_tabulated_2013_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=1,
                                             particleChargeAU=1)

NM64_alpha_filepath_tabulated_2013 = pkg_resources.resource_filename(__name__,'NM64_alpha_response_tabulated_2013.csv')
NM64_alpha_response_tabulated_2013_epn = pd.read_csv(NM64_alpha_filepath_tabulated_2013,header=None)
NM64_alpha_response_tabulated_2013_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_alpha_response_tabulated_2013_epn["Yield_per_m2_per_sr"] = NM64_alpha_response_tabulated_2013_epn["Yield_arb_units"] * 8.37e-2 * 4 #need to multiply by 4 because at the moment its in per nucleon
NM64_alpha_response_tabulated_2013 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_alpha_response_tabulated_2013_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_tabulated_2013_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)

# 2020 response values and the response function below are taken from https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA027433

NM64_alpha_filepath_tabulated_2020 = pkg_resources.resource_filename(__name__,'NM64_alpha_response_tabulated_2020.csv')
NM64_alpha_response_tabulated_2020_epn = pd.read_csv(NM64_alpha_filepath_tabulated_2020,header=None)
NM64_alpha_response_tabulated_2020_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_alpha_response_tabulated_2020_epn["Yield_per_m2_per_sr"] = NM64_alpha_response_tabulated_2020_epn["Yield_arb_units"] * 4 #need to multiply by 4 because at the moment its in per nucleon
NM64_alpha_response_tabulated_2020 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_alpha_response_tabulated_2020_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_tabulated_2020_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)

@njit
def convert_particle_energy_to_rigidity(atomic_mass, atomic_charge, particle_name, energy_GeV_per_nucleon):

    rest_mass_GeV_dict = {"proton":0.938,"alpha":3.727}

    R_value = (atomic_mass / atomic_charge) * np.sqrt(energy_GeV_per_nucleon * (energy_GeV_per_nucleon + (2 * rest_mass_GeV_dict[particle_name])))

    return R_value

@np.vectorize
def get_NM64_response_value(particle_name:str,energy_GeV_per_nucleon:float):
    
    dict_of_parameters={"proton":{
                                  0.0: {"a0":-12.104,"a1":13.879,"a2":-8.6616,"a3":0.0,},
                                  1.28:{"a0":-8.76,"a1":2.831,"a2":0.428,"a3":-0.186,},
                                  10.0:  {"a0":-4.763,"a1":1.206,"a2":-0.0365,"a3":0.0,}
                                  },
                        "alpha":{
                                  0.0: {"a0":-9.5396,"a1":7.2827,"a2":-3.0659,"a3":0.5082,},
                                  1.7:{"a0":-8.65,"a1":4.9329,"a2":-1.2022,"a3":0.1179,},
                                  15.0:  {"a0":-4.763,"a1":1.206,"a2":-0.0365,"a3":0.0,}
                                  }
                        }
    
    energy_limits_keys = list(dict_of_parameters[particle_name].keys())

    if energy_GeV_per_nucleon <= energy_limits_keys[1]:
        params_to_use = dict_of_parameters[particle_name][energy_limits_keys[0]]
    elif energy_GeV_per_nucleon < energy_limits_keys[2]:
        params_to_use = dict_of_parameters[particle_name][energy_limits_keys[1]]
    else:
        params_to_use = dict_of_parameters[particle_name][energy_limits_keys[2]]

    if particle_name == "proton":
        R_value = convert_particle_energy_to_rigidity(1,1,"proton",energy_GeV_per_nucleon)
    elif particle_name == "alpha":
        R_value = convert_particle_energy_to_rigidity(4,2,"alpha",energy_GeV_per_nucleon)

    ln_yield_value = sum_yield_values(list(params_to_use.values()), R_value)

    return np.exp(ln_yield_value)

@njit
def sum_yield_values(params_to_use, R_value):
    ln_yield_value = 0.0
    for l in [0,1,2,3]:
        ln_yield_value += params_to_use[l] * (np.log(R_value)**l)
    return ln_yield_value

def get_parameter_value(particle_name:str,parameter_label:str,energy_GeV_per_nucleon:float):

    dict_of_parameters={
                        "proton":{
                                  "A":{"b5":6.945e-9,"b4":-1.461e-7,"b3":1.115e-6,"b2":-3.402e-6,"b1":3.355e-6,"b0":-9.823e-7},
                                  "B":{"b5":-3.963e-6,"b4":8.091e-5,"b3":-6.394e-4,"b2":2.348e-3,"b1":-4.713e-3,"b0":1.186e-2}
                                 },
                        "alpha":{
                                  "A":{"b5":9.422e-9,"b4":-2.284e-7,"b3":2.037e-6,"b2":-7.828e-6,"b1":1.203e-5,"b0":-5.545e-6},
                                  "B":{"b5":-5.351e-6,"b4":1.316e-4,"b3":-1.226e-3,"b2":5.176e-3,"b1":-1.017e-2,"b0":1.458e-2}
                                 },
                        }

    params_to_use = dict_of_parameters[particle_name][parameter_label]
    
    if particle_name == "proton":
        #R_value = PRCT.convertParticleEnergyToRigidity(energy_GeV_per_nucleon * 1000)
        R_value = convert_particle_energy_to_rigidity(1,1,"proton",energy_GeV_per_nucleon)
    elif particle_name == "alpha":
        #R_value = PRCT.convertParticleEnergyToRigidity(energy_GeV_per_nucleon * 1000 * 4, particleMassAU=4, particleChargeAU=2)
        R_value = convert_particle_energy_to_rigidity(4,2,"alpha",energy_GeV_per_nucleon)
    #R_value = np.sqrt(energy_GeV_per_nucleon * (energy_GeV_per_nucleon + 1.876))

    output_parameter = 0.0
    for l in [0,1,2,3,4,5]:
        output_parameter += params_to_use[f"b{l}"] * (np.log(R_value)**l)

    return output_parameter

#@njit
def get_NM64_response_value_atmospheric_depth(particle_name:str,energy_GeV_per_nucleon:float,atmospheric_depth_g_per_cm2):

    ground_level_values = get_NM64_response_value(particle_name,energy_GeV_per_nucleon)

    A = get_parameter_value(particle_name, "A", energy_GeV_per_nucleon)
    B = get_parameter_value(particle_name, "B", energy_GeV_per_nucleon)

    values_at_depth = get_values_at_depth(atmospheric_depth_g_per_cm2, ground_level_values, A, B)

    return values_at_depth

@njit
def get_values_at_depth(atmospheric_depth_g_per_cm2, ground_level_values, A, B):
    values_at_depth = ground_level_values * np.exp((A * (1000-atmospheric_depth_g_per_cm2)**2) + (B * (1000-atmospheric_depth_g_per_cm2)) )
    return values_at_depth

def get_NM64_response_value_altitude(particle_name:str,energy_GeV_per_nucleon:float,altitude_in_km:float):

    pressure_at_altitude_hectopascal = height_to_pressure_std(altitude_in_km * units("km"))
    atmospheric_depth_g_per_cm2 = (pressure_at_altitude_hectopascal / g).to("g / (cm**2)").magnitude

    values_at_altitude = get_NM64_response_value_atmospheric_depth(particle_name, energy_GeV_per_nucleon, atmospheric_depth_g_per_cm2)

    return values_at_altitude

#response_value_epns = list(range(0.1,1000,(1000-0.1)/1000))
#response_value_epns = np.geomspace(0.1,1000,100000)
global response_value_epns
response_value_epns = np.geomspace(0.1,1000,10000)

proton_response_values = get_NM64_response_value("proton",response_value_epns)

NM64_proton_response_tabulated_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
                                                                "Yield_per_m2_per_sr":proton_response_values})
NM64_proton_response_tabulated_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_proton_response_tabulated_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                             fluxInEnergyMeVform=NM64_proton_response_tabulated_2020_functional_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=1,
                                             particleChargeAU=1)

alpha_response_values = get_NM64_response_value("alpha",response_value_epns)

NM64_alpha_response_tabulated_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
              "Yield_per_m2_per_sr":alpha_response_values * 4})
NM64_alpha_response_tabulated_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_alpha_response_tabulated_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_tabulated_2020_functional_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)