import numpy as np
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)
from atmosphericRadiationDoseAndFlux import doseAndFluxCalculator as DAFcalc
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import ParticleRigidityCalculationTools as PRCT

from .data.NM64_responses import get_NM64_response_value_altitude, response_value_epns
from .AsymptoticDirectionProcessing import acquireWeightingFactors, get_apply_method
from .spectralCalculations.particleDistribution import particleDistribution

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
                                                            self.particle_distribution)
        #######################################################

        # print("df_with_weighting_factors:")
        # print(df_with_weighting_factors.to_string())
        # import pickle as pkl
        # with open("weightingFactorsDF_2.pkl","wb") as weightingFactorsFile:
        #     pkl.dump(df_with_weighting_factors, weightingFactorsFile)

        #######################################################
        df_with_weighting_factors_full_angles.to_pickle("df_with_weighting_factors_full_angles.pkl")

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

        return sortedOutputDoseRates

    def calc_output_dose_flux(self, asymp_dir_DF_with_weighting_factors, list_of_altitudes_in_km, particle_name):

        spectrum_to_function_conversion_function = lambda row: interp1d(x=np.array(row["Rigidity"]),
                                                                        y=np.array(row["fullRigidityPitchWeightingFactor"]),
                                                                        bounds_error=False,
                                                                        fill_value=0.0)

        DFofSpectraForEachCoord = asymp_dir_DF_with_weighting_factors.groupby(["initialLatitude","initialLongitude"]).apply(spectrum_to_function_conversion_function)
        outputDoseRatesForAltitudeRange = get_apply_method(DFofSpectraForEachCoord)(lambda spectrum:DAFcalc.calculate_from_rigidity_spec(
                                                                                                inputRigidityDistributionFunctionGV=lambda x:float(spectrum(x)),
                                                                                                altitudesInkm=list_of_altitudes_in_km,
                                                                                                particleName = particle_name,
                                                                                                )
                                                                            )
        
        outputDoseRatesOnlyDF = pd.concat(outputDoseRatesForAltitudeRange.tolist(), ignore_index=True)

        index_list = [-2,-1,0,1,2,3,4,5,6,7,8]

        if self.generate_NM_count_rates == True:
            print("calculating neutron monitor count rates...")
            outputDoseRatesOnlyDF["NM64_cr_unnorm"] = self.get_neutron_monitor_count_rates(list_of_altitudes_in_km, particle_name, DFofSpectraForEachCoord)
            print("neutron monitor count rates successfully determined.")

            index_list = [-2,-1,0,1,2,3,4,5,6,7,8,9]
        
        outputDoseRatesOnlyDF_with_lats_and_longs = self.reinsert_original_lats_and_longs(outputDoseRatesOnlyDF,outputDoseRatesForAltitudeRange)

        sortedOutputDoseRates = outputDoseRatesOnlyDF_with_lats_and_longs  \
                                        .iloc[:,index_list] \
                                        .sort_values(["latitude","longitude","altitude (km)"], ignore_index=True)
        return sortedOutputDoseRates

    def get_neutron_monitor_count_rates(self, list_of_altitudes_in_km, particle_name, DFofSpectraForEachCoord):
        NM64_values_list = []

        for altitude_in_km in list_of_altitudes_in_km:
            NM64_vals_DF = get_apply_method(DFofSpectraForEachCoord)(lambda x:self.calculate_unnormed_NM64_cr(x,particle_name = particle_name, altitude_in_km=altitude_in_km)).values
            NM64_values_list.append(NM64_vals_DF)

        new_cr_unnorm_list = np.array(NM64_values_list).T.flatten()
        return new_cr_unnorm_list

    def calculate_unnormed_NM64_cr(self, interpolated_spectra_row, particle_name, altitude_in_km = 0.0):
        response_values = get_NM64_response_value_altitude(particle_name, response_value_epns, altitude_in_km)

        NM64_response_tabulated_2020_functional_epn = pd.DataFrame({"Energy_per_nucleon_GeV_per_n":response_value_epns,
                                                                        "Yield_per_m2_per_sr":response_values})

        if particle_name == "proton":

            NM64_response_tabulated_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_response_tabulated_2020_functional_epn["Energy_per_nucleon_GeV_per_n"]*1000,
                                                    fluxInEnergyMeVform=NM64_response_tabulated_2020_functional_epn["Yield_per_m2_per_sr"],
                                                    particleMassAU=1,
                                                    particleChargeAU=1)
            
        elif particle_name == "alpha":

            NM64_response_tabulated_2020_functional = PRCT.convertParticleEnergySpecToRigiditySpec(NM64_response_tabulated_2020_functional_epn["Energy_per_nucleon_GeV_per_n"] * 1000 * 4,
                                                    fluxInEnergyMeVform=NM64_response_tabulated_2020_functional_epn["Yield_per_m2_per_sr"],
                                                    particleMassAU=4,
                                                    particleChargeAU=2)

        else:
            raise Exception("ERROR: particle name did not match either proton or alpha!")
        
        NM64_response_to_use = NM64_response_tabulated_2020_functional

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