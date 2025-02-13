import numpy as np
import datetime as dt
import spaceweather as sw
import pandas as pd
from typing import Callable, List, Optional

from .utils import get_correctly_formatted_particle_dist_list, get_kp_index, validate_altitudes
from .anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import DLRmodelSpectrum, CommonModifiedPowerLawSpectrum, CommonModifiedPowerLawSpectrumSplit
from .anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import gaussianBeeckPitchAngleDistribution, isotropicPitchAngleDistribution, gaussianPitchAngleDistribution
from .anisotropic_MAIRE_engine.generalEngineInstance import generalEngineInstance, default_array_of_lats_and_longs

def run_from_spectra(
        proton_rigidity_spectrum: Optional[Callable[[float], float]] = None,
        alpha_rigidity_spectrum: Optional[Callable[[float], float]] = None,
        reference_pitch_angle_latitude: float = 0.0,
        reference_pitch_angle_longitude: float = 45.0,
        proton_pitch_angle_distribution: Callable[[float, float], float] = isotropicPitchAngleDistribution(),
        alpha_pitch_angle_distribution: Callable[[float, float], float] = isotropicPitchAngleDistribution(),
        altitudes_in_kft: Optional[List[float]] = None,
        altitudes_in_km: Optional[List[float]] = None,
        Kp_index: Optional[int] = None,
        date_and_time: dt.datetime = dt.datetime.utcnow(),
        array_of_lats_and_longs: np.ndarray = default_array_of_lats_and_longs,
        cache_magnetocosmics_run: bool = True,
        generate_NM_count_rates: bool = False,
        use_default_9_zeniths_azimuths: bool = False,
        asymp_dir_file: Optional[str] = None,
        **mag_cos_kwargs,
) -> pd.DataFrame:
    """
    Perform a run to calculate dose rates across Earth's atmosphere based on proton, alpha particle + heavier ions, or proton + alpha particle + heavier ions spectra.

    Parameters:
    - proton_rigidity_spectrum: callable, optional
        Function describing the proton rigidity spectrum in units of cm-2 s-1 sr-1 (GV/n)-1.
    - alpha_rigidity_spectrum: callable, optional
        Function describing the alpha particle rigidity spectrum in units of cm-2 s-1 sr-1 (GV/n)-1.
    - reference_pitch_angle_latitude: float, optional
        Reference latitude in GEO coordinates representing a pitch angle of 0 in the supplied pitch angle distribution.
    - reference_pitch_angle_longitude: float, optional
        Reference longitude in GEO coordinates representing a pitch angle of 0 in the supplied pitch angle distribution.
    - proton_pitch_angle_distribution: callable, optional
        Function describing the proton pitch angle distribution.
    - alpha_pitch_angle_distribution: callable, optional
        Function describing the alpha particle pitch angle distribution.
    - altitudes_in_kft: list, optional
        List of altitudes in kilofeet to perform calculations for.
    - altitudes_in_km: list, optional
        List of altitudes in kilometers to perform calculations for.
    - Kp_index: int, optional
        Kp index representing geomagnetic conditions.
    - date_and_time: datetime, optional
        Date and time for the simulation.
    - array_of_lats_and_longs: array, optional
        Array of latitudes and longitudes to perform calculations for.
    - cache_magnetocosmics_run: bool, optional
        Whether to cache the results of MAGNETOCOSMICS simulations.
    - generate_NM_count_rates: bool, optional
        Whether to generate neutron monitor count rates.
    - use_default_9_zeniths_azimuths: bool, optional
        Whether to use the mean of nine different asymptotic directions to calculate dose rates.
    - **mag_cos_kwargs: additional keyword arguments
        Additional arguments to pass to AsympDirsCalculator.

    Returns:
    - output_dose_rate_DF: DataFrame
        DataFrame containing the calculated dose rates.
    """
    
    if Kp_index is None:
        Kp_index = get_kp_index(date_and_time)
    
    altitudes_in_km = validate_altitudes(altitudes_in_km, altitudes_in_kft)

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
                                          generate_NM_count_rates=generate_NM_count_rates,
                                          asymp_dir_file=asymp_dir_file)
    
    output_dose_rate_DF = engine_to_run.getAsymptoticDirsAndRun(use_default_9_zeniths_azimuths,**mag_cos_kwargs)

    print("Success!")

    return output_dose_rate_DF

def run_from_power_law_gaussian_distribution(
        J0: float, gamma: float, deltaGamma: float, sigma: float, 
        reference_pitch_angle_latitude: float, reference_pitch_angle_longitude: float, 
        Kp_index: int, date_and_time: dt.datetime,
        use_split_spectrum: bool = False,
        asymp_dir_file: Optional[str] = None,
        **kwargs
) -> pd.DataFrame:
    """
    Perform a run to calculate dose rates using a combined power law rigidity spectrum and Gaussian pitch angle distribution.

    Parameters:
    - J0: float
        Normalization factor for the rigidity spectrum.
    - gamma: float
        Spectral index for the rigidity spectrum.
    - deltaGamma: float
        Modification factor for the spectral index.
    - sigma: float
        Standard deviation for the Gaussian pitch angle distribution.
    - reference_pitch_angle_latitude: float
        Reference latitude in GEO coordinates representing a pitch angle of 0.
    - reference_pitch_angle_longitude: float
        Reference longitude in GEO coordinates representing a pitch angle of 0.
    - Kp_index: int
        Kp index representing geomagnetic conditions.
    - date_and_time: datetime
        Date and time for the simulation.
    - use_split_spectrum: bool, optional
        Whether to use a split spectrum.
    - **kwargs: additional keyword arguments
        Additional arguments to pass to run_from_spectra.

    Returns:
    - output_dose_rate_DF: DataFrame
        DataFrame containing the calculated dose rates.
    """
    spec_to_use = CommonModifiedPowerLawSpectrumSplit if use_split_spectrum else CommonModifiedPowerLawSpectrum

    return run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianPitchAngleDistribution(normFactor=1,sigma=sigma),
        Kp_index=Kp_index,date_and_time=date_and_time,
        asymp_dir_file=asymp_dir_file,
        **kwargs,
    )

def run_from_double_power_law_gaussian_distribution(
        J0: float, gamma: float, deltaGamma: float, sigma_1: float, sigma_2: float,
        B: float, alpha_prime: float,
        reference_pitch_angle_latitude: float, reference_pitch_angle_longitude: float, 
        Kp_index: int, date_and_time: dt.datetime,
        use_split_spectrum: bool = False,
        asymp_dir_file: Optional[str] = None,
        **kwargs
) -> pd.DataFrame:
    """
    Perform a run to calculate dose rates using a double power law rigidity spectrum and Gaussian pitch angle distribution.

    Parameters:
    - J0: float
        Normalization factor for the rigidity spectrum.
    - gamma: float
        Spectral index for the rigidity spectrum.
    - deltaGamma: float
        Modification factor for the spectral index.
    - sigma_1: float
        Standard deviation for the first Gaussian pitch angle distribution.
    - sigma_2: float
        Standard deviation for the second Gaussian pitch angle distribution.
    - B: float
        Scaling factor for the second Gaussian pitch angle distribution.
    - alpha_prime: float
        Shift factor for the second Gaussian pitch angle distribution.
    - reference_pitch_angle_latitude: float
        Reference latitude in GEO coordinates representing a pitch angle of 0.
    - reference_pitch_angle_longitude: float
        Reference longitude in GEO coordinates representing a pitch angle of 0.
    - Kp_index: int
        Kp index representing geomagnetic conditions.
    - date_and_time: datetime
        Date and time for the simulation.
    - use_split_spectrum: bool, optional
        Whether to use a split spectrum.
    - **kwargs: additional keyword arguments
        Additional arguments to pass to run_from_spectra.

    Returns:
    - output_dose_rate_DF: DataFrame
        DataFrame containing the calculated dose rates.
    """
    spec_to_use = CommonModifiedPowerLawSpectrumSplit if use_split_spectrum else lambda J0,gamma,deltaGamma:CommonModifiedPowerLawSpectrum(J0,gamma,deltaGamma, lowerLimit=0.814529,upperLimit=21.084584)

    return run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianPitchAngleDistribution(normFactor=1,sigma=sigma_1) + (B * gaussianPitchAngleDistribution(normFactor=1,sigma=sigma_2,alpha=alpha_prime)),
        Kp_index=Kp_index,date_and_time=date_and_time,
        asymp_dir_file=asymp_dir_file,
        **kwargs,
    )

def run_from_power_law_Beeck_gaussian_distribution(
        J0: float, gamma: float, deltaGamma: float, A: float, B: float, 
        reference_pitch_angle_latitude: float, reference_pitch_angle_longitude: float, 
        Kp_index: int, date_and_time: dt.datetime,
        use_split_spectrum: bool = False,
        asymp_dir_file: Optional[str] = None,
        **kwargs
) -> pd.DataFrame:
    """
    Perform a run to calculate dose rates using a power law rigidity spectrum and Beeck Gaussian pitch angle distribution.

    Parameters:
    - J0: float
        Normalization factor for the rigidity spectrum.
    - gamma: float
        Spectral index for the rigidity spectrum.
    - deltaGamma: float
        Modification factor for the spectral index.
    - A: float
        Parameter A for the Beeck Gaussian pitch angle distribution.
    - B: float
        Parameter B for the Beeck Gaussian pitch angle distribution.
    - reference_pitch_angle_latitude: float
        Reference latitude in GEO coordinates representing a pitch angle of 0.
    - reference_pitch_angle_longitude: float
        Reference longitude in GEO coordinates representing a pitch angle of 0.
    - Kp_index: int
        Kp index representing geomagnetic conditions.
    - date_and_time: datetime
        Date and time for the simulation.
    - use_split_spectrum: bool, optional
        Whether to use a split spectrum.
    - **kwargs: additional keyword arguments
        Additional arguments to pass to run_from_spectra.

    Returns:
    - output_dose_rate_DF: DataFrame
        DataFrame containing the calculated dose rates.
    """
    spec_to_use = CommonModifiedPowerLawSpectrum if use_split_spectrum else CommonModifiedPowerLawSpectrumSplit

    return run_from_spectra(
        proton_rigidity_spectrum=spec_to_use(J0, gamma, deltaGamma),
        reference_pitch_angle_latitude=reference_pitch_angle_latitude,
        reference_pitch_angle_longitude=reference_pitch_angle_longitude,
        proton_pitch_angle_distribution=gaussianBeeckPitchAngleDistribution(normFactor=1,A=A,B=B),
        Kp_index=Kp_index,date_and_time=date_and_time,
        asymp_dir_file=asymp_dir_file,
        **kwargs,
    )

def run_from_DLR_cosmic_ray_model(
        OULU_count_rate_in_seconds: Optional[float] = None,
        W_parameter: Optional[float] = None,
        Kp_index: Optional[int] = None,
        date_and_time: Optional[dt.datetime] = None,
        asymp_dir_file: Optional[str] = None,
        **kwargs
) -> pd.DataFrame:
    """
    Perform a run to calculate dose rates using the DLR cosmic ray model.

    Parameters:
    - OULU_count_rate_in_seconds: float, optional
        OULU count rate in seconds.
    - W_parameter: float, optional
        W parameter for the DLR model.
    - Kp_index: int, optional
        Kp index representing geomagnetic conditions.
    - date_and_time: datetime, optional
        Date and time for the simulation.
    - **kwargs: additional keyword arguments
        Additional arguments to pass to run_from_spectra.

    Returns:
    - output_dose_rate_DF: DataFrame
        DataFrame containing the calculated dose rates.
    """
    if (W_parameter is None) and (OULU_count_rate_in_seconds is None):
        print("As no OULU count rates or W parameters were specified, the count rate of OULU in the past will be determined using the supplied date and time for the purposes of calculating the incoming spectra.")
        DLR_model_date_and_time = date_and_time
    else:
        DLR_model_date_and_time = None

    return run_from_spectra(
        proton_rigidity_spectrum=DLRmodelSpectrum(atomicNumber=1, date_and_time=DLR_model_date_and_time, OULUcountRateInSeconds=OULU_count_rate_in_seconds, W_parameter=W_parameter),
        alpha_rigidity_spectrum=DLRmodelSpectrum(atomicNumber=2, date_and_time=DLR_model_date_and_time, OULUcountRateInSeconds=OULU_count_rate_in_seconds, W_parameter=W_parameter),
        Kp_index=Kp_index,date_and_time=date_and_time,
        asymp_dir_file=asymp_dir_file,
        **kwargs,
    )





