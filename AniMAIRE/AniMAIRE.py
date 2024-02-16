import numpy as np
import datetime as dt
import plotly.express as px

from .utils import get_correctly_formatted_particle_dist_list
from .anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import DLRmodelSpectrum, MishevModifiedPowerLawSpectrum, MishevModifiedPowerLawSpectrumSplit
from .anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import gaussianBeeckPitchAngleDistribution, isotropicPitchAngleDistribution, gaussianPitchAngleDistribution

from .anisotropic_MAIRE_engine.maireSengine import generalEngineInstance, default_array_of_lats_and_longs

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
        date_and_time=dt.datetime.now(),
        array_of_lats_and_longs=default_array_of_lats_and_longs,
        cache_magnetocosmics_run=True,
        generate_NM_count_rates=True,
        **mag_cos_kwargs,
):
    
    if Kp_index is None:
        raise Exception("Error: no Kp index specified!")
    
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
    
    output_dose_rate_DF = engine_to_run.getAsymptoticDirsAndRun(**mag_cos_kwargs)

    print("Success!")

    return output_dose_rate_DF

def run_from_power_law_gaussian_distribution(J0, gamma, deltaGamma, sigma, 
                                             reference_pitch_angle_latitude, reference_pitch_angle_longitude, 
                                             Kp_index,date_and_time,
                                             use_split_spectrum=False,
                                             **kwargs):
    
    if use_split_spectrum == True:
        spec_to_use = MishevModifiedPowerLawSpectrumSplit
    else:
        spec_to_use = MishevModifiedPowerLawSpectrum
        

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
        spec_to_use = MishevModifiedPowerLawSpectrumSplit
    else:
        spec_to_use = lambda J0,gamma,deltaGamma:MishevModifiedPowerLawSpectrum(J0,gamma,deltaGamma, lowerLimit=0.814529,upperLimit=21.084584)
        

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
        spec_to_use = MishevModifiedPowerLawSpectrum
    else:
        spec_to_use = MishevModifiedPowerLawSpectrumSplit

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

def create_single_dose_map_plot(DF_to_use,
                              selected_altitude_in_km):

    if selected_altitude_in_km is not None:
        DF_to_use = DF_to_use[round(DF_to_use["altitude (km)"],4) == selected_altitude_in_km]

    if len(DF_to_use) == 0:
        raise Exception("Error: specified altitude in kilometers did not match any of the altitudes in kilometers in the inputted DataFrame.")

    doseRateMap = px.scatter(DF_to_use, x="longitude",y="latitude",color="adose",
                            symbol_sequence=["square"],
                            range_y=[-90,90],
                            range_x=[0,360])
    
    doseRateMap.update_traces(marker={'size': 10})
    doseRateMap.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1, range=[-90,90]))
    doseRateMap.update_xaxes(range=[0,360])
    doseRateMap.update_yaxes(range=[-90,90])
    doseRateMap.update_layout(xaxis_scaleanchor="y")

    doseRateMap.update_layout(autosize=False,
                            width=800,
                            height=600)

    doseRateMap.show()

    return doseRateMap

import seaborn as sns
import geopandas
import matplotlib.pyplot as plt

def create_single_dose_map_plot_plt(heatmap_DF_to_Plot, hue_range = (2,7), heatmap_s = 63,edgecolor='white',
                    dose_type = "adose"):

    ############################ creating background world map and dose image
    currentFigure = plt.gcf()

    currentFigure.set_figheight(10)
    currentFigure.set_figwidth(10)

    #heatmap_DF_to_Plot = pd.read_csv(file_path_to_read, delimiter=',')
    heatmap_DF_to_Plot["SEU (Upsets/hr/Gb)"] = heatmap_DF_to_Plot["SEU"] * (60.0 * 60.0) * 1e9
    heatmap_DF_to_Plot["SEL (Latch-ups/hr/device)"] = heatmap_DF_to_Plot["SEL"] * (60.0 * 60.0)
    heatmap_DF_to_Plot["longitudeTranslated"] = heatmap_DF_to_Plot["longitude"].apply(lambda x:x-360.0 if x > 180.0 else x)

    scatterPlotAxis = sns.scatterplot(data=heatmap_DF_to_Plot,x="longitudeTranslated",y="latitude",
                    hue=dose_type, hue_norm=hue_range, palette="crest",
                    zorder=10,
                    marker="s",s=heatmap_s,edgecolor=edgecolor,
                    legend="brief",
                    )#ax=axToPlotOn)

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world.plot(color="None",edgecolor="black",lw=1.5,ax=scatterPlotAxis,zorder=20)

    ####################################################################

    plt.ylim([-90,90])
    plt.xlim([-175,180])
    plt.grid(True)
    plt.xlabel("Longitude (degrees)")
    plt.ylabel("Latitude (degrees)")

    plt.legend(title=dose_type,loc="center left",bbox_to_anchor=(1.1,0.5))

    return plt.gca()


    
