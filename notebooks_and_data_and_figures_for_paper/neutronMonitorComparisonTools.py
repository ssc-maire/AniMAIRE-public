import plotting_modules
import datetime as dt
import kpindex
import geomagindices as gi
import spaceweather as sw

from AniMAIRE import AniMAIRE
import numpy as np
import pandas as pd

import plotly.graph_objects as go
import matplotlib.pyplot as plt

from AsympDirsCalculator import AsympDirsTools

import scipy
from scipy import optimize as spoptimize

from tabula.io import read_pdf

from joblib import Memory
#cache_location = './GLE71_nm_comparisons'
cache_location = './GLE71_nm_comparisons_new'
cache_location_flat = './GLE71_nm_comparisons_new_flat'
memory = Memory(cache_location, verbose=0)
memory_adjusted_flat_production = Memory(cache_location_flat, verbose=0)

default_list_of_zeniths_and_azimuths = [
                                        [0.0,0.0],
                                        [16.0,0.0],
                                        [16.0,90.0],
                                        [16.0,180.0],
                                        [16.0,270.0],
                                        [32.0,0.0],
                                        [32.0,90.0],
                                        [32.0,180.0],
                                        [32.0,270.0],
                                        ]

def get_Kp_index_for_date_and_time(date_and_time):
    kp_index_data = sw.ap_kp_3h(update=True)

    Kp_index = int(kp_index_data[
        (kp_index_data.index > (date_and_time-dt.timedelta(hours=3)).strftime("%m-%d-%Y %H:%M:%S")) &
        (kp_index_data.index <= date_and_time.strftime("%m-%d-%Y %H:%M:%S"))
        ]["Kp"].iloc[0])
    
    return Kp_index

def get_nm_pre_increase_cr(datetime_for_event, station_code):

    nm_data_for_datetime = plotting_modules.getAllStationDataForDatetime(datetime_for_event)

    pre_increase_cr = nm_data_for_datetime[nm_data_for_datetime["STATION"] == station_code]["preincreaseCountRate(cts/s)"].iloc[0]

    return pre_increase_cr

def collect_relevant_input_data(datetime_for_event):
    nm_data_for_datetime = plotting_modules.getAllStationDataForDatetime(datetime_for_event)

    OULU_pre_increase_cr = nm_data_for_datetime[nm_data_for_datetime["STATION"] == "OULU"]["preincreaseCountRate(cts/s)"].iloc[0]

    Kp_index_to_use = get_Kp_index_for_date_and_time(datetime_for_event)

    array_of_lats_and_longs = nm_data_for_datetime[["latitude","longitude"]].applymap(float).to_numpy()

    array_of_altitudes_in_km = nm_data_for_datetime["altitude(m)"].apply(float).unique() / 1000.0
    return nm_data_for_datetime,OULU_pre_increase_cr,Kp_index_to_use,array_of_lats_and_longs,array_of_altitudes_in_km

def importParamsFromMishevPDF():

    GLE59_date = dt.datetime(year=2000,month=7,day=14,hour=0,minute=0,tzinfo=dt.timezone.utc)

    specDetailsPDF = read_pdf("MishevGLEfittingPaper.pdf",pages=5)
    print(specDetailsPDF)
    parametersTable = specDetailsPDF[0]
    parametersTable.dropna(axis=1,inplace=True)
    parametersTable.columns = ["Time","J0","gamma","deltaGamma","sigmaSquared","Phi","Lambda"]
    parametersTable["J0"] = parametersTable["J0"].apply(lambda s:s.replace(",",""))
    parametersTable["sigmaSquared"] = parametersTable["sigmaSquared"].apply(lambda s:s.replace(" ",""))
    parametersTable["Time"] = parametersTable["Time"].apply(lambda s:s.replace(":","_"))
    parametersTable.iloc[:,1:] = parametersTable.iloc[:,1:].applymap(float)
    parametersTable["Lambda"] = -1 * parametersTable["Lambda"]

    parametersTable["J0"] = parametersTable["J0"]

    parametersTable["datetime"] = parametersTable["Time"].apply(lambda x:GLE59_date + dt.timedelta(seconds=(float(x.split("_")[0])*60*60) + \
                                                                                            (float(x.split("_")[1])*60)))

    return parametersTable.drop("Time",axis=1)

def importParamsFromMishevPDF2():

    GLE59_date = dt.datetime(year=2000,month=7,day=14,hour=0,minute=0,tzinfo=dt.timezone.utc)

    specDetailsPDF = read_pdf("MishevGLEfittingPaper2.pdf",pages=6)
    print(specDetailsPDF)
    parametersTable = specDetailsPDF[0]
    parametersTable.dropna(axis=1,inplace=True)
    parametersTable = pd.DataFrame(np.array(list(map(lambda x:" ".join(x).split(" "), list(parametersTable.iloc[1:].values))))[:,[0,3,4,5,6,7,8]])
    parametersTable.columns = ["Time","J0","gamma","deltaGamma","sigmaSquared","Phi","Lambda"]
    print(parametersTable)
    parametersTable["J0"] = parametersTable["J0"].apply(lambda s:s.replace(",",""))
    parametersTable["sigmaSquared"] = parametersTable["sigmaSquared"].apply(lambda s:s.replace(" ",""))
    parametersTable["Time"] = parametersTable["Time"].apply(lambda s:s.replace(":","_"))
    parametersTable["Lambda"] = parametersTable["Lambda"].apply(lambda s:s.replace('\U00002212',"-"))
    parametersTable.iloc[:,1:] = parametersTable.iloc[:,1:].applymap(float)

    parametersTable["J0"] = parametersTable["J0"]

    parametersTable["datetime"] = parametersTable["Time"].apply(lambda x:GLE59_date + dt.timedelta(seconds=(float(x.split("_")[0])*60*60) + \
                                                                                            (float(x.split("_")[1])*60)))

    return parametersTable.drop("Time",axis=1)

def importParamsFromMishevPDF2_GLE70():

    GLE70_date = dt.datetime(year=2006,month=12,day=13,hour=0,minute=0,tzinfo=dt.timezone.utc)

    specDetailsPDF = read_pdf("MishevGLEfittingPaper2.pdf",pages=10)
    print(specDetailsPDF)
    parametersTable = specDetailsPDF[0]
    parametersTable.dropna(axis=1,inplace=True)
    parametersTable = pd.DataFrame(np.array(list(map(lambda x:" ".join(x).split(" "), list(parametersTable.iloc[1:].values))))[:,[0,3,4,5,6,7,8]])
    parametersTable.columns = ["Time","J0","gamma","deltaGamma","sigmaSquared","Phi","Lambda"]
    print(parametersTable)
    parametersTable["J0"] = parametersTable["J0"].apply(lambda s:s.replace(",",""))
    parametersTable["sigmaSquared"] = parametersTable["sigmaSquared"].apply(lambda s:s.replace(" ",""))
    parametersTable["Time"] = parametersTable["Time"].apply(lambda s:s.replace(":","_"))
    parametersTable["Lambda"] = parametersTable["Lambda"].apply(lambda s:s.replace('\U00002212',"-"))
    parametersTable["Phi"] = parametersTable["Phi"].apply(lambda s:s.replace('\U00002212',"-"))
    parametersTable.iloc[:,1:] = parametersTable.iloc[:,1:].applymap(float)

    parametersTable["J0"] = parametersTable["J0"]

    parametersTable["datetime"] = parametersTable["Time"].apply(lambda x:GLE70_date + dt.timedelta(seconds=(float(x.split("_")[0])*60*60) + \
                                                                                            (float(x.split("_")[1])*60)))

    return parametersTable.drop("Time",axis=1)

def importParamsFromMishevPDF_GLE71():

    GLE71_date = dt.datetime(year=2012,month=5,day=17,hour=0,minute=0,tzinfo=dt.timezone.utc)

    specDetailsPDF = pd.read_html("MishevGLE71fittingPaper.html")
    print(specDetailsPDF)
    parametersTable = specDetailsPDF[1]
    parametersTable.dropna(axis=1,inplace=True)
    #parametersTable = pd.DataFrame(np.array(list(map(lambda x:" ".join(x).split(" "), list(parametersTable.iloc[1:].values))))[:,[0,3,4,5,6,7,8]])
    parametersTable.columns = ["Time","J0","gamma","deltaGamma","sigma_1_Squared","B","sigma_2_Squared","alpha_prime","Phi","Lambda","D"]
    print(parametersTable)
    #parametersTable["J0"] = parametersTable["J0"].apply(lambda s:s.replace(",",""))
    #parametersTable["sigma_1_Squared"] = parametersTable["sigma_1_Squared"].apply(lambda s:s.replace(" ",""))
    #parametersTable["sigma_2_Squared"] = parametersTable["sigma_2_Squared"].apply(lambda s:s.replace(" ",""))
    #parametersTable["Time"] = parametersTable["Time"].apply(lambda s:s.replace(":","_").replace('\U00002212',"-").replace("–","-").split("-")[0])
    parametersTable["Time"] = parametersTable["Time"].apply(lambda s:s.replace(":","_").replace('\U00002212',"-").replace("–","-").split("-")[1])
    #parametersTable["Lambda"] = parametersTable["Lambda"].apply(lambda s:s.replace('\U00002212',"-"))
    parametersTable["Phi"] = parametersTable["Phi"].apply(lambda s:s.replace('\U00002212',"-"))
    parametersTable["gamma"] = parametersTable["gamma"].apply(lambda s:s.replace('\U00002212',"-"))

    parametersTable["gamma"] = parametersTable["gamma"].apply(float) * -1

    parametersTable.iloc[:,1:] = parametersTable.iloc[:,1:].applymap(float)

    parametersTable["J0"] = parametersTable["J0"]

    parametersTable["datetime"] = parametersTable["Time"].apply(lambda x:GLE71_date + dt.timedelta(seconds=(float(x.split("_")[0])*60*60) + \
                                                                                            (float(x.split("_")[1])*60)))

    return parametersTable.drop(["Time","D"],axis=1)

def get_params_GLE73():
    beginning_datetime_for_event = dt.datetime(year=2021,month=10,day=28,hour=15,minute=55,second=0,tzinfo=dt.timezone.utc)
    params_to_use = pd.read_html("https://link.springer.com/article/10.1007/s11207-022-02026-0/tables/2")[0]
    params_to_use["datetime"] = params_to_use["Integration interval UT"].apply(lambda x:x.split("–")[0]) \
                                            .apply(lambda x:dt.timedelta(hours=int(x.split(":")[0]) - beginning_datetime_for_event.hour, 
                                                                        minutes=int(x.split(":")[1])  - beginning_datetime_for_event.minute)) \
                                            .apply(lambda x:x + beginning_datetime_for_event)
    
    params_to_use.drop('Integration interval UT',inplace=True, axis=1)
    params_to_use.rename(columns={'\(J_{0}\) [m−2 s−1 sr−1 GV−1]':'J0'},inplace=True)
    params_to_use.rename(columns={'γ':"gamma"},inplace=True)
    params_to_use.rename(columns={'δγ':"deltaGamma"},inplace=True)
    params_to_use.rename(columns={'\(\sigma ^{2}\) [rad2]':"sigmaSquared"},inplace=True)
    params_to_use.rename(columns={'Ψ [degrees]':"Phi"},inplace=True)
    params_to_use.rename(columns={'Λ [degrees]':"Lambda"},inplace=True)

    #params_to_use["J0"] = params_to_use["J0"].apply(lambda s:s.replace(",",""))
    #params_to_use["sigmaSquared"] = params_to_use["sigmaSquared"].apply(lambda s:s.replace(" ",""))
    params_to_use["Lambda"] = params_to_use["Lambda"].apply(lambda s:s.replace('\U00002212',"-"))
    params_to_use["Phi"] = params_to_use["Phi"].apply(lambda s:s.replace('\U00002212',"-"))
    params_to_use[["J0","gamma","deltaGamma","sigmaSquared","Phi","Lambda"]] = params_to_use[["J0","gamma","deltaGamma","sigmaSquared","Phi","Lambda"]].applymap(float)


    return params_to_use

def mean_value_of_PAD(sigma):
    # numerator = sigma * (1-np.exp(-(np.pi**2)/(sigma**2)))
    # denominator = np.sqrt(np.pi) * scipy.special.erf(np.pi/sigma)

    # mean_pitch_angle = numerator / denominator

    # mean_value = np.exp(-(mean_pitch_angle**2)/(sigma**2))

    integral_of_PAD = 0.5 * np.sqrt(np.pi) * sigma * scipy.special.erf(np.pi/sigma)

    return integral_of_PAD / np.pi

list_of_dose_columns = ["edose","adose","dosee","tn1","tn2","tn3","SEU","SEL","NM64_cr_unnorm"]
list_of_flat_dose_columns = ["flat_edose","flat_adose","flat_dosee","flat_tn1","flat_tn2","flat_tn3","flat_SEU","flat_SEL","flat_NM64_cr_unnorm"]

list_of_dose_columns_no_NM = ["edose","adose","dosee","tn1","tn2","tn3","SEU","SEL"]
list_of_flat_dose_columns_no_NM = ["flat_edose","flat_adose","flat_dosee","flat_tn1","flat_tn2","flat_tn3","flat_SEU","flat_SEL"]

highestMaxRigValue=1010
maxRigValue=20
minRigValue=0.1
#minRigValue=1
nIncrements_high=60
nIncrements_low=200
#nIncrements_low=2000
#nIncrements_low=20000
n_zeniths_azimuths = 1

# highestMaxRigValue=20
# maxRigValue=19.9
# minRigValue=0.1
# nIncrements_high=2
# nIncrements_low=200
# #nIncrements_low=25
# n_zeniths_azimuths = 1

array_of_zeniths_and_azimuths = np.array([(i,j) for i in np.linspace(0,30,n_zeniths_azimuths) for j in np.linspace(0,360,n_zeniths_azimuths)])

def get_percentage_increase_at_NM(date_and_time_to_run_for, 
                                  Mishev_params_to_run_for,station_code, 
                                  single_station_run_mode=False, 
                                  anisotropy_mode = "power_law_gaussian",
                                  location_lat_displacement = 0.0,
                                  location_long_displacement = 0.0,
                                  **kwargs):

    if single_station_run_mode:
        nm_data, nm_DLR_dose_rate, nm_Mishev_dose_rate, flat_nm_Mishev_dose_rate = \
            determine_all_relevant_dose_rates(date_and_time_to_run_for, 
                                              Mishev_params_to_run_for, 
                                              station_to_run_for="all", 
                                              anisotropy_mode=anisotropy_mode,
                                              location_lat_displacement = location_lat_displacement,
                                              location_long_displacement = location_long_displacement,
                                              **kwargs)
    else:
        nm_data, nm_DLR_dose_rate, nm_Mishev_dose_rate, flat_nm_Mishev_dose_rate = \
            determine_all_relevant_dose_rates(date_and_time_to_run_for, 
                                              Mishev_params_to_run_for, 
                                              station_to_run_for=station_code, 
                                              anisotropy_mode=anisotropy_mode,
                                              location_lat_displacement = location_lat_displacement,
                                              location_long_displacement = location_long_displacement,
                                              **kwargs)

    station_data_to_run_for = nm_data.query(f"STATION == '{station_code}'")

    for dose_rate_DF in [nm_DLR_dose_rate, nm_Mishev_dose_rate, flat_nm_Mishev_dose_rate]:
        dose_rate_DF["latitude"] = dose_rate_DF["latitude"] - location_lat_displacement
        dose_rate_DF["longitude"] = dose_rate_DF["longitude"] - location_long_displacement

    nm_DLR_dose_rate_single = find_dose_rate_line(nm_DLR_dose_rate, station_data_to_run_for)
    nm_Mishev_dose_rate_single = find_dose_rate_line(nm_Mishev_dose_rate, station_data_to_run_for)
    flat_nm_Mishev_dose_rate_single = find_dose_rate_line(flat_nm_Mishev_dose_rate, station_data_to_run_for)

    percentage_increase_rates = nm_Mishev_dose_rate_single.copy()
    percentage_increase_rates["datetime"] = date_and_time_to_run_for
    percentage_increase_rates["nm_percent_increase"] = station_data_to_run_for["% INC."].iloc[0]
    try:
        percentage_increase_rates[list_of_dose_columns] = (nm_Mishev_dose_rate_single[list_of_dose_columns] / nm_DLR_dose_rate_single[list_of_dose_columns]) * 100.0
        percentage_increase_rates[list_of_flat_dose_columns] = (flat_nm_Mishev_dose_rate_single[list_of_dose_columns] / nm_DLR_dose_rate_single[list_of_dose_columns]) * 100.0
    except KeyError:
        percentage_increase_rates[list_of_dose_columns_no_NM] = (nm_Mishev_dose_rate_single[list_of_dose_columns_no_NM] / nm_DLR_dose_rate_single[list_of_dose_columns_no_NM]) * 100.0
        percentage_increase_rates[list_of_flat_dose_columns_no_NM] = (flat_nm_Mishev_dose_rate_single[list_of_dose_columns_no_NM] / nm_DLR_dose_rate_single[list_of_dose_columns_no_NM]) * 100.0

    #return percentage_increase_rates, nm_DLR_dose_rate_single, nm_Mishev_dose_rate_single, flat_nm_Mishev_dose_rate_single
    return percentage_increase_rates #, nm_DLR_dose_rate_single


# # #@testmemory.cache()
# # #@memoryjacobian.cache()
# #@memory.cache()
# def determine_all_relevant_dose_rates(date_and_time_to_run_for, 
#                                       Mishev_params_to_run_for, 
#                                       station_to_run_for="all", 
#                                       anisotropy_mode = "power_law_gaussian",
#                                       custom_GCR_proton_spectrum = None,
#                                       custom_GCR_alpha_spectrum = None,
#                                       multiply_sigmasquares_by_factor = 1,
#                                       lats_and_longs_to_run_for=None,
#                                       altitudes_to_run_for_in_km=None,
#                                       **kwargs):
    
#     nm_data, OULU_pre_increase, Kp_index_at_time, array_of_all_lats_and_longs, array_of_alts = collect_relevant_input_data(date_and_time_to_run_for)

#     if station_to_run_for is not "all":
#         nm_data = nm_data[nm_data["STATION"] == station_to_run_for]

#     #Kp_index_at_time = 1

#     #output_nm_pre_increase_cr = get_nm_pre_increase_cr(date_and_time_to_run_for, station_code)

#     if lats_and_longs_to_run_for is None:
#         lats_and_longs_to_run_for = nm_data[["latitude","longitude"]].drop_duplicates().to_numpy()
#     if lats_and_longs_to_run_for=="default_MAIRE":
#         lats_and_longs_to_run_for = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

#     if altitudes_to_run_for_in_km is None:
#         altitudes_to_run_for_in_km = (nm_data["altitude(m)"] / 1000.0).unique()
#     if altitudes_to_run_for_in_km == "default_MAIRE":
#         altitudes_to_run_for_in_km = [0,10,20] + [i for i in range(25, 61 + 1, 3)]

#     if custom_GCR_proton_spectrum is None:
#         nm_DLR_dose_rate = AniMAIRE.run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=OULU_pre_increase,
#                                                     Kp_index=Kp_index_at_time,
#                                                     date_and_time=date_and_time_to_run_for,
#                                                     array_of_lats_and_longs=lats_and_longs_to_run_for,
#                                                     altitudes_in_km=altitudes_to_run_for_in_km,
#                                                     **kwargs)
#     else:
#         nm_DLR_dose_rate = AniMAIRE.run_from_spectra(
#                                                     proton_rigidity_spectrum=custom_GCR_proton_spectrum,
#                                                     alpha_rigidity_spectrum=custom_GCR_alpha_spectrum,
#                                                     Kp_index=Kp_index_at_time,
#                                                     date_and_time=date_and_time_to_run_for,
#                                                     array_of_lats_and_longs=lats_and_longs_to_run_for,
#                                                     altitudes_in_km=altitudes_to_run_for_in_km,
#                                                     **kwargs)

#     if anisotropy_mode == "power_law_gaussian":
#         nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
#             J0=Mishev_params_to_run_for["J0"],
#             gamma=Mishev_params_to_run_for["gamma"],
#             deltaGamma=Mishev_params_to_run_for["deltaGamma"],
#             sigma=np.sqrt(Mishev_params_to_run_for["sigmaSquared"]),
#             reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
#             reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
#             Kp_index=Kp_index_at_time,
#             date_and_time=date_and_time_to_run_for,
#             array_of_lats_and_longs=lats_and_longs_to_run_for,
#             altitudes_in_km=altitudes_to_run_for_in_km,

#             highestMaxRigValue=highestMaxRigValue,
#             maxRigValue=maxRigValue,
#             minRigValue=minRigValue,
#             nIncrements_high=nIncrements_high,
#             nIncrements_low=nIncrements_low,

#             #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
#             **kwargs
#         )

#         flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigmaSquared"]) * 100_000
#         flat_J0_to_use = Mishev_params_to_run_for["J0"] * mean_value_of_PAD(Mishev_params_to_run_for["sigmaSquared"])

#         flat_nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
#             J0=flat_J0_to_use,
#             gamma=Mishev_params_to_run_for["gamma"],
#             deltaGamma=Mishev_params_to_run_for["deltaGamma"],
#             sigma=flat_sigma_to_use,
#             reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
#             reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
#             Kp_index=Kp_index_at_time,
#             date_and_time=date_and_time_to_run_for,
#             array_of_lats_and_longs=lats_and_longs_to_run_for,
#             altitudes_in_km=altitudes_to_run_for_in_km,

#             highestMaxRigValue=highestMaxRigValue,
#             maxRigValue=maxRigValue,
#             minRigValue=minRigValue,
#             nIncrements_high=nIncrements_high,
#             nIncrements_low=nIncrements_low,

#             **kwargs
#             #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
#         )

#     elif anisotropy_mode == "double_power_law_gaussian":
#         nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
#             J0=Mishev_params_to_run_for["J0"],
#             gamma=Mishev_params_to_run_for["gamma"],
#             deltaGamma=Mishev_params_to_run_for["deltaGamma"],
#             sigma_1=np.sqrt(multiply_sigmasquares_by_factor * Mishev_params_to_run_for["sigma_1_Squared"]),
#             sigma_2=np.sqrt(multiply_sigmasquares_by_factor * Mishev_params_to_run_for["sigma_2_Squared"]),
#             B=Mishev_params_to_run_for["B"],
#             alpha_prime=Mishev_params_to_run_for["alpha_prime"],
#             reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
#             reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
#             Kp_index=Kp_index_at_time,
#             date_and_time=date_and_time_to_run_for,
#             array_of_lats_and_longs=lats_and_longs_to_run_for,
#             altitudes_in_km=altitudes_to_run_for_in_km,

#             highestMaxRigValue=highestMaxRigValue,
#             maxRigValue=maxRigValue,
#             minRigValue=minRigValue,
#             nIncrements_high=nIncrements_high,
#             nIncrements_low=nIncrements_low,

#             **kwargs
#             #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
#         )

#         flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigma_1_Squared"]) * 100_000
#         flat_J0_to_use = Mishev_params_to_run_for["J0"] * mean_value_of_PAD(Mishev_params_to_run_for["sigma_1_Squared"])

#         flat_nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
#             J0=flat_J0_to_use,
#             gamma=Mishev_params_to_run_for["gamma"],
#             deltaGamma=Mishev_params_to_run_for["deltaGamma"],
#             sigma_1=flat_sigma_to_use,
#             sigma_2=flat_sigma_to_use,
#             B=Mishev_params_to_run_for["B"],
#             alpha_prime=Mishev_params_to_run_for["alpha_prime"],
#             reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
#             reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
#             Kp_index=Kp_index_at_time,
#             date_and_time=date_and_time_to_run_for,
#             array_of_lats_and_longs=lats_and_longs_to_run_for,
#             altitudes_in_km=altitudes_to_run_for_in_km,

#             highestMaxRigValue=highestMaxRigValue,
#             maxRigValue=maxRigValue,
#             minRigValue=minRigValue,
#             nIncrements_high=nIncrements_high,
#             nIncrements_low=nIncrements_low,

#             **kwargs
#             #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
#         )

#     else:
#         raise Exception("Error: anisotropy mode not among valid options!")
#     #exit()
#     #indfidfbdbudufduf

    
    
#     return nm_data,nm_DLR_dose_rate,nm_Mishev_dose_rate,flat_nm_Mishev_dose_rate

def determine_all_relevant_dose_rates(date_and_time_to_run_for, 
                                      Mishev_params_to_run_for, 
                                      station_to_run_for="all", 
                                      anisotropy_mode = "power_law_gaussian",
                                      custom_GCR_proton_spectrum = None,
                                      custom_GCR_alpha_spectrum = None,
                                      multiply_sigmasquares_by_factor = 1,
                                      lats_and_longs_to_run_for=None,
                                      altitudes_to_run_for_in_km=None,
                                      phi_displacement=0.0,
                                      lambda_displacement=0.0,
                                      location_lat_displacement = 0.0,
                                      location_long_displacement = 0.0,
                                      **kwargs):
    
    nm_data, OULU_pre_increase, Kp_index_at_time, array_of_all_lats_and_longs, array_of_alts = collect_relevant_input_data(date_and_time_to_run_for)

    if station_to_run_for != "all":
        nm_data = nm_data[nm_data["STATION"] == station_to_run_for]

    #Kp_index_at_time = 1

    #output_nm_pre_increase_cr = get_nm_pre_increase_cr(date_and_time_to_run_for, station_code)

    if lats_and_longs_to_run_for is None:
        lats_and_longs_to_run_for = nm_data[["latitude","longitude"]].drop_duplicates().to_numpy()
    if lats_and_longs_to_run_for=="default_MAIRE":
        lats_and_longs_to_run_for = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

    lats_and_longs_to_run_for = lats_and_longs_to_run_for + np.array([[location_lat_displacement,location_long_displacement]])

    if altitudes_to_run_for_in_km is None:
        altitudes_to_run_for_in_km = (nm_data["altitude(m)"] / 1000.0).unique()
    if altitudes_to_run_for_in_km == "default_MAIRE":
        altitudes_to_run_for_in_km = np.array([0,10,20] + [i for i in range(25, 61 + 1, 3)]) * 0.3048

    if custom_GCR_proton_spectrum is None:
        nm_DLR_dose_rate = AniMAIRE.run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=OULU_pre_increase,
                                                    Kp_index=Kp_index_at_time,
                                                    date_and_time=date_and_time_to_run_for,
                                                    array_of_lats_and_longs=lats_and_longs_to_run_for,
                                                    altitudes_in_km=altitudes_to_run_for_in_km,
                                                    **kwargs)
    else:
        nm_DLR_dose_rate = AniMAIRE.run_from_spectra(
                                                    proton_rigidity_spectrum=custom_GCR_proton_spectrum,
                                                    alpha_rigidity_spectrum=custom_GCR_alpha_spectrum,
                                                    Kp_index=Kp_index_at_time,
                                                    date_and_time=date_and_time_to_run_for,
                                                    array_of_lats_and_longs=lats_and_longs_to_run_for,
                                                    altitudes_in_km=altitudes_to_run_for_in_km,
                                                    **kwargs)

    if anisotropy_mode == "power_law_gaussian":

        flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigmaSquared"]) * 100_000
        flat_J0_to_use = Mishev_params_to_run_for["J0"] * mean_value_of_PAD(Mishev_params_to_run_for["sigmaSquared"])

        flat_nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
            J0=flat_J0_to_use,
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma=flat_sigma_to_use,
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

        nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
            J0=Mishev_params_to_run_for["J0"],
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma=np.sqrt(Mishev_params_to_run_for["sigmaSquared"]),
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
            **kwargs
        )

    elif anisotropy_mode == "double_power_law_gaussian":

        flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigma_1_Squared"]) * 100_000
        flat_J0_to_use = Mishev_params_to_run_for["J0"] #* mean_value_of_PAD(Mishev_params_to_run_for["sigma_1_Squared"])

        flat_nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
            J0=flat_J0_to_use,
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma_1=flat_sigma_to_use,
            sigma_2=flat_sigma_to_use,
            B=0.0, #Mishev_params_to_run_for["B"],
            alpha_prime=Mishev_params_to_run_for["alpha_prime"],
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"] + phi_displacement,
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"] + lambda_displacement,
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

        nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
            J0=Mishev_params_to_run_for["J0"],
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma_1=np.sqrt(Mishev_params_to_run_for["sigma_1_Squared"]),
            sigma_2=np.sqrt(multiply_sigmasquares_by_factor * Mishev_params_to_run_for["sigma_2_Squared"]),
            B=Mishev_params_to_run_for["B"],
            alpha_prime=Mishev_params_to_run_for["alpha_prime"],
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"] + phi_displacement,
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"] + lambda_displacement,
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

    else:
        raise Exception("Error: anisotropy mode not among valid options!")
    #exit()
    #indfidfbdbudufduf

    #print(nm_Mishev_dose_rate)
    
    return nm_data,nm_DLR_dose_rate,nm_Mishev_dose_rate,flat_nm_Mishev_dose_rate

@memory.cache()
def determine_all_relevant_dose_rates_GCR_only(date_and_time_to_run_for, 
                                      Mishev_params_to_run_for, 
                                      station_to_run_for="all", 
                                      anisotropy_mode = "power_law_gaussian",
                                      custom_GCR_proton_spectrum = None,
                                      custom_GCR_alpha_spectrum = None,
                                      multiply_sigmasquares_by_factor = 1,
                                      lats_and_longs_to_run_for=None,
                                      altitudes_to_run_for_in_km=None,
                                      **kwargs):
    
    nm_data, OULU_pre_increase, Kp_index_at_time, array_of_all_lats_and_longs, array_of_alts = collect_relevant_input_data(date_and_time_to_run_for)

    if station_to_run_for != "all":
        nm_data = nm_data[nm_data["STATION"] == station_to_run_for]

    #Kp_index_at_time = 1

    #output_nm_pre_increase_cr = get_nm_pre_increase_cr(date_and_time_to_run_for, station_code)

    if lats_and_longs_to_run_for is None:
        lats_and_longs_to_run_for = nm_data[["latitude","longitude"]].drop_duplicates().to_numpy()
    if lats_and_longs_to_run_for=="default_MAIRE":
        lats_and_longs_to_run_for = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

    if altitudes_to_run_for_in_km is None:
        altitudes_to_run_for_in_km = (nm_data["altitude(m)"] / 1000.0).unique()
    if altitudes_to_run_for_in_km == "default_MAIRE":
        altitudes_to_run_for_in_km = np.array([0,10,20] + [i for i in range(25, 61 + 1, 3)]) * 0.3048

    if custom_GCR_proton_spectrum is None:
        nm_DLR_dose_rate = AniMAIRE.run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=OULU_pre_increase,
                                                    Kp_index=Kp_index_at_time,
                                                    date_and_time=date_and_time_to_run_for,
                                                    array_of_lats_and_longs=lats_and_longs_to_run_for,
                                                    altitudes_in_km=altitudes_to_run_for_in_km,
                                                    **kwargs)
    else:
        nm_DLR_dose_rate = AniMAIRE.run_from_spectra(
                                                    proton_rigidity_spectrum=custom_GCR_proton_spectrum,
                                                    alpha_rigidity_spectrum=custom_GCR_alpha_spectrum,
                                                    Kp_index=Kp_index_at_time,
                                                    date_and_time=date_and_time_to_run_for,
                                                    array_of_lats_and_longs=lats_and_longs_to_run_for,
                                                    altitudes_in_km=altitudes_to_run_for_in_km,
                                                    **kwargs)
    
    return nm_DLR_dose_rate

@memory.cache()
def determine_all_relevant_dose_rates_GLE_only(date_and_time_to_run_for, 
                                      Mishev_params_to_run_for, 
                                      station_to_run_for="all", 
                                      anisotropy_mode = "power_law_gaussian",
                                      custom_GCR_proton_spectrum = None,
                                      custom_GCR_alpha_spectrum = None,
                                      multiply_sigmasquares_by_factor = 1,
                                      lats_and_longs_to_run_for=None,
                                      altitudes_to_run_for_in_km=None,
                                      **kwargs):
    
    nm_data, OULU_pre_increase, Kp_index_at_time, array_of_all_lats_and_longs, array_of_alts = collect_relevant_input_data(date_and_time_to_run_for)

    if station_to_run_for != "all":
        nm_data = nm_data[nm_data["STATION"] == station_to_run_for]

    #Kp_index_at_time = 1

    #output_nm_pre_increase_cr = get_nm_pre_increase_cr(date_and_time_to_run_for, station_code)

    if lats_and_longs_to_run_for is None:
        lats_and_longs_to_run_for = nm_data[["latitude","longitude"]].drop_duplicates().to_numpy()
    if lats_and_longs_to_run_for=="default_MAIRE":
        lats_and_longs_to_run_for = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

    if altitudes_to_run_for_in_km is None:
        altitudes_to_run_for_in_km = (nm_data["altitude(m)"] / 1000.0).unique()
    if altitudes_to_run_for_in_km == "default_MAIRE":
        altitudes_to_run_for_in_km = np.array([0,10,20] + [i for i in range(25, 61 + 1, 3)]) * 0.3048

    if anisotropy_mode == "power_law_gaussian":
        nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
            J0=Mishev_params_to_run_for["J0"],
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma=np.sqrt(Mishev_params_to_run_for["sigmaSquared"]),
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
            **kwargs
        )

    elif anisotropy_mode == "double_power_law_gaussian":
        nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
            J0=Mishev_params_to_run_for["J0"],
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma_1=np.sqrt(multiply_sigmasquares_by_factor * Mishev_params_to_run_for["sigma_1_Squared"]),
            sigma_2=np.sqrt(multiply_sigmasquares_by_factor * Mishev_params_to_run_for["sigma_2_Squared"]),
            B=Mishev_params_to_run_for["B"],
            alpha_prime=Mishev_params_to_run_for["alpha_prime"],
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

    else:
        raise Exception("Error: anisotropy mode not among valid options!")
    #exit()
    #indfidfbdbudufduf

    return nm_Mishev_dose_rate

#@memory.cache()
@memory_adjusted_flat_production.cache()
def determine_all_flat_relevant_dose_rates(date_and_time_to_run_for, 
                                      Mishev_params_to_run_for, 
                                      station_to_run_for="all", 
                                      anisotropy_mode = "power_law_gaussian",
                                      custom_GCR_proton_spectrum = None,
                                      custom_GCR_alpha_spectrum = None,
                                      multiply_sigmasquares_by_factor = 1,
                                      lats_and_longs_to_run_for=None,
                                      altitudes_to_run_for_in_km=None,
                                      **kwargs):
    
    nm_data, OULU_pre_increase, Kp_index_at_time, array_of_all_lats_and_longs, array_of_alts = collect_relevant_input_data(date_and_time_to_run_for)

    if station_to_run_for != "all":
        nm_data = nm_data[nm_data["STATION"] == station_to_run_for]

    #Kp_index_at_time = 1

    #output_nm_pre_increase_cr = get_nm_pre_increase_cr(date_and_time_to_run_for, station_code)

    if lats_and_longs_to_run_for is None:
        lats_and_longs_to_run_for = nm_data[["latitude","longitude"]].drop_duplicates().to_numpy()
    if lats_and_longs_to_run_for=="default_MAIRE":
        lats_and_longs_to_run_for = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

    if altitudes_to_run_for_in_km is None:
        altitudes_to_run_for_in_km = (nm_data["altitude(m)"] / 1000.0).unique()
    if altitudes_to_run_for_in_km == "default_MAIRE":
        altitudes_to_run_for_in_km = np.array([0,10,20] + [i for i in range(25, 61 + 1, 3)]) * 0.3048

    if anisotropy_mode == "power_law_gaussian":

        flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigmaSquared"]) * 100_000
        flat_J0_to_use = Mishev_params_to_run_for["J0"] #* mean_value_of_PAD(Mishev_params_to_run_for["sigmaSquared"])

        flat_nm_Mishev_dose_rate = AniMAIRE.run_from_power_law_gaussian_distribution(
            J0=flat_J0_to_use,
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma=flat_sigma_to_use,
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

    elif anisotropy_mode == "double_power_law_gaussian":

        flat_sigma_to_use = np.sqrt(Mishev_params_to_run_for["sigma_1_Squared"]) * 100_000
        flat_J0_to_use = Mishev_params_to_run_for["J0"] #* mean_value_of_PAD(Mishev_params_to_run_for["sigma_1_Squared"])

        flat_nm_Mishev_dose_rate = AniMAIRE.run_from_double_power_law_gaussian_distribution(
            J0=flat_J0_to_use,
            gamma=Mishev_params_to_run_for["gamma"],
            deltaGamma=Mishev_params_to_run_for["deltaGamma"],
            sigma_1=flat_sigma_to_use,
            sigma_2=flat_sigma_to_use,
            B=0.0, #Mishev_params_to_run_for["B"],
            alpha_prime=Mishev_params_to_run_for["alpha_prime"],
            reference_pitch_angle_latitude=Mishev_params_to_run_for["Phi"],
            reference_pitch_angle_longitude=Mishev_params_to_run_for["Lambda"],
            Kp_index=Kp_index_at_time,
            date_and_time=date_and_time_to_run_for,
            array_of_lats_and_longs=lats_and_longs_to_run_for,
            altitudes_in_km=altitudes_to_run_for_in_km,

            highestMaxRigValue=highestMaxRigValue,
            maxRigValue=maxRigValue,
            minRigValue=minRigValue,
            nIncrements_high=nIncrements_high,
            nIncrements_low=nIncrements_low,

            **kwargs
            #array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
        )

    else:
        raise Exception("Error: anisotropy mode not among valid options!")
    #exit()
    #indfidfbdbudufduf

    
    
    return flat_nm_Mishev_dose_rate

def find_dose_rate_line(nm_dose_rates, station_data_to_run_for):
    nm_dose_rate_single = nm_dose_rates[(nm_dose_rates["latitude"] == station_data_to_run_for["latitude"].iloc[0]) &
                                               (nm_dose_rates["longitude"] == station_data_to_run_for["longitude"].iloc[0]) &
                                               ((nm_dose_rates["altitude (km)"] * 1000).apply(round) == round(station_data_to_run_for["altitude(m)"].iloc[0]))]
    return nm_dose_rate_single

def reverse_symmetry_dir(inputDF):
    outputDF = inputDF.copy()

    outputDF["Phi"] = inputDF["Phi"] * -1
    outputDF["Lambda"] = inputDF["Lambda"] + 180.0

    return outputDF

def get_percent_increases_for_many_datetimes(station_to_run_for, 
                                             params_to_run_across, 
                                             single_station_run_mode=False, 
                                             anisotropy_mode = "power_law_gaussian",
                                             **kwargs):
    percent_increases_list = []

    for index, Mishev_params_for_datetime in params_to_run_across.iterrows():

        percent_increase_for_datetime = get_percentage_increase_at_NM(Mishev_params_for_datetime["datetime"], 
                                                                      Mishev_params_for_datetime,station_to_run_for, 
                                                                      single_station_run_mode=single_station_run_mode,
                                                                      anisotropy_mode=anisotropy_mode,
                                                                      **kwargs)
        percent_increases_list.append(percent_increase_for_datetime)

    percent_increases_DF = pd.concat(percent_increases_list).reset_index(drop=True)

    percent_increases_DF["STATION"] = station_to_run_for

    return percent_increases_DF

def divide_dose_rates_by_factor(df_of_dose_rates,factor):

    output_df = df_of_dose_rates.copy()

    output_df[list_of_dose_columns + list_of_flat_dose_columns] = output_df[list_of_dose_columns + list_of_flat_dose_columns] / factor

    return output_df

def divide_flat_dose_rates_by_factor(df_of_dose_rates,factor):

    output_df = df_of_dose_rates.copy()

    output_df[list_of_flat_dose_columns] = output_df[list_of_flat_dose_columns] / factor

    return output_df

def plot_nm_dose_rate_comp(percent_increase_DF_raw,divide_MAIRE_data_by = 1,metric="NM64_cr_unnorm",plot_log=False):

    percent_increase_DF = divide_dose_rates_by_factor(percent_increase_DF_raw,divide_MAIRE_data_by)

    station_name = percent_increase_DF["STATION"].iloc[0]
    longitude = percent_increase_DF["longitude"].iloc[0]
    latitude = percent_increase_DF["latitude"].iloc[0]
    altitude_km = percent_increase_DF["altitude (km)"].iloc[0]

    title_string = f"{station_name}\n" + \
              f"longitude = {longitude}°E, latitude = {latitude}°N, altitude = {altitude_km} km"

    percent_increase_DF.plot(x="datetime",
                             y=[f"{metric}",f"flat_{metric}","nm_percent_increase"],
                             marker="o",ls="-")
    plt.grid(True)
    plt.title(title_string)
    plt.ylabel("percentage increase")

    plt.legend(["tn3, anisotropic","tn3, isotropic","actual percentage increase"])

    if plot_log == True:

        percent_increase_DF.plot(x="datetime",y=[f"{metric}",f"flat_{metric}","nm_percent_increase"],marker="o")
        plt.yscale("log")
        plt.grid(True)

        plt.title(title_string)
    