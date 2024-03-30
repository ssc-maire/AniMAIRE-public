import numpy as np
import datetime as dt
import spaceweather as sw
import pandas as pd

from .utils import get_correctly_formatted_particle_dist_list
from .anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import DLRmodelSpectrum, CommonModifiedPowerLawSpectrum, CommonModifiedPowerLawSpectrumSplit
from .anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import gaussianBeeckPitchAngleDistribution, isotropicPitchAngleDistribution, gaussianPitchAngleDistribution
from .anisotropic_MAIRE_engine.generalEngineInstance import generalEngineInstance, default_array_of_lats_and_longs

default_altitudes_in_kft = [0,10,20] + [i for i in range(25, 61 + 1, 3)]

def run_from_spectra(
        proton_rigidity_spectrum=None,
        alpha_rigidity_spectrum=None,
        reference_pitch_angle_latitude=0.0, #None,
        reference_pitch_angle_longitude=45.0, #None,
        proton_pitch_angle_distribution=isotropicPitchAngleDistribution(),
        alpha_pitch_angle_distribution=isotropicPitchAngleDistribution(),
        altitudes_in_kft=default_altitudes_in_kft,
        altitudes_in_km=None,
        Kp_index=None,
        date_and_time=dt.datetime.utcnow(),
        array_of_lats_and_longs=default_array_of_lats_and_longs,
        cache_magnetocosmics_run=True,
        generate_NM_count_rates=False,
        use_default_9_zeniths_azimuths=False,
        **mag_cos_kwargs,
):
    
    if Kp_index is None:
        Kp_index = sw.ap_kp_3h()[pd.to_datetime(sw.ap_kp_3h().index) < date_and_time].iloc[-1]["Kp"]
        #raise Exception("Error: no Kp index specified!")
    
    if (altitudes_in_km is not None) and (altitudes_in_kft is not default_altitudes_in_kft):
        raise Exception("Error: only one of altitudes_in_km and altitudes_in_kft should be supplied!")
    
    if altitudes_in_km is None:
        altitudes_in_km = np.array(altitudes_in_kft) * 0.3048

    if (proton_rigidity_spectrum is None) and (alpha_rigidity_spectrum is None):
        raise Exception("Error: either a proton rigidity spectrum or an alpha rigidity spectrum must be specified!")
    
    list_of_particle_distributions = get_correctly_formatted_particle_dist_list(proton_rigidity_spectrum, 
                                                                                alpha_rigidity_spectrum, 
                                                                                reference_pitch_angle_latitude, 
                                                                                reference_pitch_angle_longitude, 
                                                                                proton_pitch_angle_distribution, 
                                                                                alpha_pitch_angle_distribution)

    
    engine_to_run = generalEngineInstance(list_of_particle_distributions,
                                          list_of_altitudes_km=altitudes_in_km,
                                          Kp_index=Kp_index,
                                          date_and_time=date_and_time,
                                          reference_latitude=reference_pitch_angle_latitude,
                                          reference_longitude=reference_pitch_angle_longitude,
                                          array_of_lats_and_longs=array_of_lats_and_longs,
                                          cache_magnetocosmics_runs=cache_magnetocosmics_run,
                                          generate_NM_count_rates=generate_NM_count_rates)
    
    output_dose_rate_DF = engine_to_run.getAsymptoticDirsAndRun(use_default_9_zeniths_azimuths,**mag_cos_kwargs)

    print("Success!")

    return output_dose_rate_DF

def run_from_power_law_gaussian_distribution(J0, gamma, deltaGamma, sigma, 
                                             reference_pitch_angle_latitude, reference_pitch_angle_longitude, 
                                             Kp_index,date_and_time,
                                             use_split_spectrum=False,
                                             **kwargs):
    
    if use_split_spectrum == True:
        spec_to_use = CommonModifiedPowerLawSpectrumSplit
    else:
        spec_to_use = CommonModifiedPowerLawSpectrum
        

    output_dose_rate_DF = run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianPitchAngleDistribution(normFactor=1,sigma=sigma),
        Kp_index=Kp_index,date_and_time=date_and_time,
        **kwargs,
    )

    return output_dose_rate_DF

def run_from_double_power_law_gaussian_distribution(J0, gamma, deltaGamma, sigma_1, sigma_2,
                                                    B, alpha_prime,
                                             reference_pitch_angle_latitude, reference_pitch_angle_longitude, 
                                             Kp_index,date_and_time,
                                             use_split_spectrum=False,
                                             **kwargs):
    
    if use_split_spectrum == True:
        spec_to_use = CommonModifiedPowerLawSpectrumSplit
    else:
        spec_to_use = lambda J0,gamma,deltaGamma:CommonModifiedPowerLawSpectrum(J0,gamma,deltaGamma, lowerLimit=0.814529,upperLimit=21.084584)
        

    output_dose_rate_DF = run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianPitchAngleDistribution(normFactor=1,sigma=sigma_1) + (B * gaussianPitchAngleDistribution(normFactor=1,sigma=sigma_2,alpha=alpha_prime)),
        Kp_index=Kp_index,date_and_time=date_and_time,
        **kwargs,
    )

    return output_dose_rate_DF

def run_from_power_law_Beeck_gaussian_distribution(J0, gamma, deltaGamma, A, B, 
                                             reference_pitch_angle_latitude, reference_pitch_angle_longitude, 
                                             Kp_index,date_and_time,
                                             use_split_spectrum=False,
                                             **kwargs):
    
    if use_split_spectrum == True:
        spec_to_use = CommonModifiedPowerLawSpectrum
    else:
        spec_to_use = CommonModifiedPowerLawSpectrumSplit

    output_dose_rate_DF = run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianBeeckPitchAngleDistribution(normFactor=1,A=A,B=B),
        Kp_index=Kp_index,date_and_time=date_and_time,
        **kwargs,
    )

    return output_dose_rate_DF

def run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=None,
                                      W_parameter=None,
                                      Kp_index=None,date_and_time=None,
                                      **kwargs):
    
    if (W_parameter is None) and (OULU_count_rate_in_seconds is None):
        print("As no OULU count rates or W parameters were specified, the count rate of OULU in the past will be determined using the supplied date and time for the purposes of calculating the incoming spectra.")
        DLR_model_date_and_time = date_and_time
    else:
        DLR_model_date_and_time = None

    output_dose_rate_DF = run_from_spectra(
        proton_rigidity_spectrum=DLRmodelSpectrum(atomicNumber=1, date_and_time=DLR_model_date_and_time, OULUcountRateInSeconds=OULU_count_rate_in_seconds, W_parameter=W_parameter),
        alpha_rigidity_spectrum=DLRmodelSpectrum(atomicNumber=2, date_and_time=DLR_model_date_and_time, OULUcountRateInSeconds=OULU_count_rate_in_seconds, W_parameter=W_parameter),
        Kp_index=Kp_index,date_and_time=date_and_time,
        **kwargs,
    )

    return output_dose_rate_DF




    
