import numpy as np
import pandas as pd

from pandarallel import pandarallel
from .utils.AsymptoticDirectionDataframe import acquireWeightingFactors, get_apply_method
pandarallel.initialize(progress_bar=True)

from .spectralCalculations.particleDistribution import particleDistribution

from atmosphericRadiationDoseAndFlux import doseAndFluxCalculator as DAFcalc
from scipy.interpolate import interp1d
from scipy.integrate import trapz

import pkg_resources

import ParticleRigidityCalculationTools as PRCT

from numba import njit, vectorize, typed, types, typeof

NM64_proton_response = pd.read_csv(pkg_resources.resource_filename(__name__,'data/NM64_proton_response.csv'), 
                                    header=None)
NM64_proton_response.columns = ["proton_rigidity_GV","cts_per_primary_per_cm2"]

NM64_proton_filepath_Mishev_2013 = pkg_resources.resource_filename(__name__,'data/NM64_proton_response_Mishev_2013.csv')
NM64_proton_response_Mishev_2013_epn = pd.read_csv(NM64_proton_filepath_Mishev_2013,header=None)
NM64_proton_response_Mishev_2013_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_proton_response_Mishev_2013_epn["Yield_per_m2_per_sr"] = NM64_proton_response_Mishev_2013_epn["Yield_arb_units"] * 0.109
NM64_proton_response_Mishev_2013 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_proton_response_Mishev_2013_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                             fluxInEnergyMeVform=NM64_proton_response_Mishev_2013_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=1,
                                             particleChargeAU=1)

NM64_alpha_filepath_Mishev_2013 = pkg_resources.resource_filename(__name__,'data/NM64_alpha_response_Mishev_2013.csv')
NM64_alpha_response_Mishev_2013_epn = pd.read_csv(NM64_alpha_filepath_Mishev_2013,header=None)
NM64_alpha_response_Mishev_2013_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_alpha_response_Mishev_2013_epn["Yield_per_m2_per_sr"] = NM64_alpha_response_Mishev_2013_epn["Yield_arb_units"] * 8.37e-2 * 4 #need to multiply by 4 because at the moment its in per nucleon
NM64_alpha_response_Mishev_2013 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_alpha_response_Mishev_2013_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_Mishev_2013_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)

# 2020 response values and the response function below are taken from https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA027433

NM64_alpha_filepath_Mishev_2020 = pkg_resources.resource_filename(__name__,'data/NM64_alpha_response_Mishev_2020.csv')
NM64_alpha_response_Mishev_2020_epn = pd.read_csv(NM64_alpha_filepath_Mishev_2020,header=None)
NM64_alpha_response_Mishev_2020_epn.columns = ["Energy_per_nucleon_GeV_per_n","Yield_arb_units"]
NM64_alpha_response_Mishev_2020_epn["Yield_per_m2_per_sr"] = NM64_alpha_response_Mishev_2020_epn["Yield_arb_units"] * 4 #need to multiply by 4 because at the moment its in per nucleon
NM64_alpha_response_Mishev_2020 = PRCT.convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV=NM64_alpha_response_Mishev_2020_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_Mishev_2020_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)

@njit
def convert_particle_energy_to_rigidity(atomic_mass, atomic_charge, particle_name, energy_GeV_per_nucleon):

    rest_mass_GeV_dict = {"proton":0.938,"alpha":3.727}

    R_value = (atomic_mass / atomic_charge) * np.sqrt(energy_GeV_per_nucleon * (energy_GeV_per_nucleon + (2 * rest_mass_GeV_dict[particle_name])))

    return R_value

# dict_of_parameter_values = typed.Dict.empty(
#                             key_type=types.unicode_type,
#                             value_type=types.float64[:],
#                         )
# dict_for_proton_limits = typed.Dict.empty(
#                         key_type=types.float64[:],
#                         value_type=typeof(dict_of_parameter_values),
#                     )

# dict_for_alpha_limits = typed.Dict.empty(
#                         key_type=types.float64[:],
#                         value_type=typeof(dict_of_parameter_values),
#                     )

# # dict_of_parameters = typed.Dict.empty(
# #                         key_type=types.unicode_type,
# #                         value_type=typeof(dict_for_limits),
# #                     )

# dict_for_proton_limits[0.0] = {"a0":-12.104,"a1":13.879,"a2":-8.6616,"a3":0.0,}

@np.vectorize
#@vectorize #(nopython=True)
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
        #R_value = PRCT.convertParticleEnergyToRigidity(energy_GeV_per_nucleon * 1000)
        R_value = convert_particle_energy_to_rigidity(1,1,"proton",energy_GeV_per_nucleon)
    elif particle_name == "alpha":
        #R_value = PRCT.convertParticleEnergyToRigidity(energy_GeV_per_nucleon * 1000 * 4, particleMassAU=4, particleChargeAU=2)
        R_value = convert_particle_energy_to_rigidity(4,2,"alpha",energy_GeV_per_nucleon)
    #R_value = np.sqrt(energy_GeV_per_nucleon * (energy_GeV_per_nucleon + 1.876))

    ln_yield_value = sum_yield_values(list(params_to_use.values()), R_value)

    return np.exp(ln_yield_value)

@njit
def sum_yield_values(params_to_use, R_value):
    ln_yield_value = 0.0
    for l in [0,1,2,3]:
        #ln_yield_value += params_to_use[f"a{l}"] * (np.log(R_value)**l)
        ln_yield_value += params_to_use[l] * (np.log(R_value)**l)
    return ln_yield_value

#@njit
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

from metpy.constants import g
from metpy.calc import height_to_pressure_std
from metpy.units import units

def get_NM64_response_value_altitude(particle_name:str,energy_GeV_per_nucleon:float,altitude_in_km:float):

    pressure_at_altitude_hectopascal = height_to_pressure_std(altitude_in_km * units("km"))
    atmospheric_depth_g_per_cm2 = (pressure_at_altitude_hectopascal / g).to("g / (cm**2)").magnitude

    values_at_altitude = get_NM64_response_value_atmospheric_depth(particle_name, energy_GeV_per_nucleon, atmospheric_depth_g_per_cm2)

    return values_at_altitude

#response_value_epns = list(range(0.1,1000,(1000-0.1)/1000))
#response_value_epns = np.geomspace(0.1,1000,100000)
response_value_epns = np.geomspace(0.1,1000,10000)

proton_response_values = get_NM64_response_value("proton",response_value_epns)

NM64_proton_response_Mishev_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
                                                                "Yield_per_m2_per_sr":proton_response_values})
NM64_proton_response_Mishev_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_proton_response_Mishev_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                             fluxInEnergyMeVform=NM64_proton_response_Mishev_2020_functional_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=1,
                                             particleChargeAU=1)

alpha_response_values = get_NM64_response_value("alpha",response_value_epns)

NM64_alpha_response_Mishev_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
              "Yield_per_m2_per_sr":alpha_response_values * 4})
NM64_alpha_response_Mishev_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_alpha_response_Mishev_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000*4,
                                             fluxInEnergyMeVform=NM64_alpha_response_Mishev_2020_functional_epn["Yield_per_m2_per_sr"],
                                             particleMassAU=4,
                                             particleChargeAU=2)

class singleParticleEngineInstance():

    def __init__(self, 
                 particle_distribution:particleDistribution, 
                 dfofAsymptoticDirections:pd.DataFrame, 
                 list_of_altitudes_in_km:list,
                 generate_NM_count_rates:bool):
        
        self.particle_distribution = particle_distribution
        self.list_of_altitudes_in_km = list_of_altitudes_in_km
        self.dfOfAllAsymptoticDirections = dfofAsymptoticDirections
        self.generate_NM_count_rates=generate_NM_count_rates
        self.rigiditySpectrumParamDict = {}
        self.pitchAngleDistributionParamDict = {}

    def getThePitchAngleDistribution(self):
        print("acquiring pitch angle distributions...")
        pitchAngleDistToUse = self.particle_distribution.momentum_distribution.pitch_angle_distribution
        try:
            print("IMF latitude set to",self.IMFlatitude)
            print("IMF longitude set to",self.IMFlongitude)
        except AttributeError:
            self.IMFlatitude = -1000.0 # distribution is probably isotropic: while a better solution is probably necessary, set to dummy variables for now.
            self.IMFlongitude = -1000.0
        pitchAngleDistToUse.setInterplanetaryMagFieldDirection(self.IMFlatitude,self.IMFlongitude)
        return pitchAngleDistToUse

    def getAsymptoticDirsAndRun(self):

        self.acquireDFofAllAsymptoticDirections()
        sortedOutputDoseRates = self.runOverSpecifiedAltitudes()

        return sortedOutputDoseRates

    def runOverSpecifiedAltitudes(self):
        
        #######################################################
        # assign each asymptotic direction a weighting factor in accordance with its pitch angle
        print("assigning pitch angle weighting factors...")
        df_with_weighting_factors_full_angles = acquireWeightingFactors(self.dfOfAllAsymptoticDirections,
                                                            self.particle_distribution.momentum_distribution)
        #######################################################

        # print("df_with_weighting_factors:")
        # print(df_with_weighting_factors.to_string())
        # import pickle as pkl
        # with open("weightingFactorsDF_2.pkl","wb") as weightingFactorsFile:
        #     pkl.dump(df_with_weighting_factors, weightingFactorsFile)

        #######################################################
        df_with_weighting_factors_full_angles.to_pickle("df_with_weighting_factors_full_angles.pkl")
        # convert the dataframe such that the maxmimum weighting factor for the initial azimuths and zeniths is used.
        #df_with_weighting_factors = get_max_weighting_factors_for_multi_angle_magcos_runs(df_with_weighting_factors_full_angles)

        df_with_weighting_factors = get_mean_weighting_factors_for_multi_angle_magcos_runs(df_with_weighting_factors_full_angles)
        #######################################################

        #######################################################
        # calculate dose rates at each point based on spectra and weighting factors
        print("converting spectra and asymptotic directions to particle fluxes and dose rates...")

        sortedOutputDoseRates = self.calc_output_dose_flux(df_with_weighting_factors, 
                                                           self.list_of_altitudes_in_km, 
                                                           self.particle_distribution.particle_species.particleName)
        
        df_with_weighting_factors.to_pickle("weighting_factor_DF.pkl")
        sortedOutputDoseRates.weighting_factor_input_DF = df_with_weighting_factors
        
        print("output dose rates calculated successfully!")

        #exit()

        #raise Exception("mdsiondsidosdids")

        # 5 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        return sortedOutputDoseRates

    def calc_output_dose_flux(self, asymp_dir_DF_with_weighting_factors, list_of_altitudes_in_km, particle_name):

        spectrum_to_function_conversion_function = lambda row: interp1d(x=np.array(row["Rigidity"]),
                                                                        y=np.array(row["fullRigidityPitchWeightingFactor"]),
                                                                        bounds_error=False,
                                                                        fill_value=0.0)

        DFofSpectraForEachCoord = asymp_dir_DF_with_weighting_factors.groupby(["initialLatitude","initialLongitude"]).apply(spectrum_to_function_conversion_function)
        #print(DFofSpectraForEachCoord)
        outputDoseRatesForAltitudeRange = get_apply_method(DFofSpectraForEachCoord)(lambda spectrum:DAFcalc.calculate_from_rigidity_spec(
                                                                                                inputRigidityDistributionFunctionGV=lambda x:float(spectrum(x)),
                                                                                                altitudesInkm=list_of_altitudes_in_km,
                                                                                                particleName = particle_name,
                                                                                                )
                                                                            )
        
        #exit()
        outputDoseRatesOnlyDF = pd.concat(outputDoseRatesForAltitudeRange.tolist(), ignore_index=True)

        index_list = [-2,-1,0,1,2,3,4,5,6,7,8]

        if self.generate_NM_count_rates == True:
            print("calculating neutron monitor count rates...")
            outputDoseRatesOnlyDF["NM64_cr_unnorm"] = self.get_neutron_monitor_count_rates(list_of_altitudes_in_km, particle_name, DFofSpectraForEachCoord)
            print("neutron monitor count rates successfully determined.")

            index_list = [-2,-1,0,1,2,3,4,5,6,7,8,9]
        
        outputDoseRatesOnlyDF_with_lats_and_longs = self.reinsert_original_lats_and_longs(outputDoseRatesOnlyDF,outputDoseRatesForAltitudeRange)
        #######################################################

        sortedOutputDoseRates = outputDoseRatesOnlyDF_with_lats_and_longs  \
                                        .iloc[:,index_list] \
                                        .sort_values(["latitude","longitude","altitude (km)"], ignore_index=True)
        return sortedOutputDoseRates

    def get_neutron_monitor_count_rates(self, list_of_altitudes_in_km, particle_name, DFofSpectraForEachCoord):
        NM64_values_list = []

        for altitude_in_km in list_of_altitudes_in_km:
            NM64_vals_DF = get_apply_method(DFofSpectraForEachCoord)(lambda x:self.calculate_unnormed_NM64_cr(x,particle_name = particle_name, altitude_in_km=altitude_in_km)).values
            NM64_values_list.append(NM64_vals_DF)

        #outputDoseRatesOnlyDF["NM64_cr_unnorm"] = np.repeat(NM64_vals_DF,len(list_of_altitudes_in_km))
        new_cr_unnorm_list = np.array(NM64_values_list).T.flatten()
        return new_cr_unnorm_list

    def calculate_unnormed_NM64_cr(self, interpolated_spectra_row, particle_name, altitude_in_km = 0.0):
        #spectra_weighting_factors_across_NM = NM64_proton_response["proton_rigidity_GV"].apply(interpolated_spectra_row)

        #response_values = get_NM64_response_value(particle_name,response_value_epns)
        response_values = get_NM64_response_value_altitude(particle_name, response_value_epns, altitude_in_km)

        NM64_response_Mishev_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
                                                                        "Yield_per_m2_per_sr":response_values})

        if particle_name == "proton":
            # #NM64_response_to_use = NM64_proton_response_Mishev_2013
            # NM64_response_to_use = NM64_response_Mishev_2020_functional

            NM64_response_Mishev_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_response_Mishev_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                                    fluxInEnergyMeVform=NM64_response_Mishev_2020_functional_epn["Yield_per_m2_per_sr"],
                                                    particleMassAU=1,
                                                    particleChargeAU=1)
            
        elif particle_name == "alpha":
            # #NM64_response_to_use = NM64_alpha_response_Mishev_2013
            # #NM64_response_to_use = NM64_alpha_response_Mishev_2020
            # NM64_response_to_use = NM64_alpha_response_Mishev_2020_functional
            # #NM64_response_to_use = NM64_proton_response_Mishev_2013
            # # NM64_response_to_use = pd.DataFrame({"Rigidity":[0.0,1.0],
            # #                                      "Rigidity distribution values":[0.0,0.0]})

            NM64_response_Mishev_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_response_Mishev_2020_functional_epn["Energy_per_nucleon_GeV_per_n"] * 1000 * 4,
                                                    fluxInEnergyMeVform=NM64_response_Mishev_2020_functional_epn["Yield_per_m2_per_sr"],
                                                    particleMassAU=4,
                                                    particleChargeAU=2)

        else:
            raise Exception("ERROR: particle name did not match either proton or alpha!")
        
        NM64_response_to_use = NM64_response_Mishev_2020_functional

        spectra_weighting_factors_across_NM = NM64_response_to_use["Rigidity"].apply(interpolated_spectra_row)
        
        full_NM64_spectra_weighting_factor = spectra_weighting_factors_across_NM * NM64_response_to_use["Rigidity distribution values"]
        output_unnormalised_NM64_cr = trapz(full_NM64_spectra_weighting_factor, NM64_response_to_use["Rigidity"])
        return output_unnormalised_NM64_cr

    def reinsert_original_lats_and_longs(self, new_doserates_DF:pd.DataFrame,original_data_frame:pd.Series):
        original_lats = original_data_frame.reset_index(level=[0,1])["initialLatitude"]
        original_longs = original_data_frame.reset_index(level=[0,1])["initialLongitude"]

        list_of_altitude_dfs = []
        for altitude_in_km in new_doserates_DF["altitude (km)"].unique():
            df_for_specific_altitude = new_doserates_DF[new_doserates_DF["altitude (km)"] == altitude_in_km].reset_index(drop=True).copy()
            df_for_specific_altitude["latitude"] = original_lats
            df_for_specific_altitude["longitude"] = original_longs
            list_of_altitude_dfs.append(df_for_specific_altitude)

        return pd.concat(list_of_altitude_dfs, ignore_index=True)
    
def get_max_weighting_factors_for_multi_angle_magcos_runs(multi_angle_DF:pd.DataFrame):

    max_weighting_factor_DF = multi_angle_DF[multi_angle_DF[["Rigidity","initialLatitude","initialLongitude","fullRigidityPitchWeightingFactor"]].groupby(["Rigidity","initialLatitude","initialLongitude"], group_keys=False).apply(lambda x:x == x.max())["fullRigidityPitchWeightingFactor"]]

    return max_weighting_factor_DF.reset_index(drop=True)

def get_mean_weighting_factors_for_multi_angle_magcos_runs(multi_angle_DF:pd.DataFrame):

    max_weighting_factor_DF = multi_angle_DF.groupby(["Rigidity","initialLatitude","initialLongitude"], group_keys=False).mean()

    return max_weighting_factor_DF.reset_index()