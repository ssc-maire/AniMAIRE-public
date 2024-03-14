import numpy as np
import pandas as pd
import datetime as dt
import sys
import ParticleRigidityCalculationTools as PRCT
from joblib import Memory
import tqdm
tqdm.tqdm.pandas()
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)
from spacepy.coordinates import Coords as spaceCoords
from spacepy.time import Ticktock as spaceTicktock
import numba

from .spectralCalculations.particleDistribution import particleDistribution

memory_asymp_dirs = Memory("./cacheAsymptoticDirectionOutputs", verbose=0)

m0 = 1.67262192e-27 #kg
c = 299792458.0 #m/s
protonCharge = 1.60217663e-19 #C
protonRestEnergy = m0 * (c**2)

global get_apply_method
def get_apply_method(DF_or_Series):
    if hasattr(sys, 'gettrace') and sys.gettrace() is not None:
        print("debug mode being used: setting AniMAIRE to use progress_apply rather than running in parallel!")
        apply_method = DF_or_Series.progress_apply
    else:
        #print("not in debug mode: setting AniMAIRE to use parallel_apply!")
        apply_method = DF_or_Series.parallel_apply
        
    return apply_method

def generate_asymp_dir_DF(dataframeToFillFrom:pd.DataFrame, IMFlatitude:float, IMFlongitude:float, datetime_to_run_across_UTC, cache:bool):

    new_asymp_dir_DF = dataframeToFillFrom.copy()

    #new_asymp_dir_DF["Energy"] = PRCT.convertParticleRigidityToEnergy(dataframeToFillFrom["Rigidity"], particleMassAU = 1, particleChargeAU = 1)

    print("assigning asymptotic coordinates")
    if cache == False:
        asymptoticDirectionList = convertAsymptoticDirectionsToPitchAngle(dataframeToFillFrom, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC)
    else:
        cachedConvertAsympDirFunc = memory_asymp_dirs.cache(convertAsymptoticDirectionsToPitchAngle)
        asymptoticDirectionList = cachedConvertAsympDirFunc(dataframeToFillFrom, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC)

    print("successfully converted asymptotic directions")

    new_asymp_dir_DF["angleBetweenIMFinRadians"] = asymptoticDirectionList

    return new_asymp_dir_DF

def convertAsymptoticDirectionsToPitchAngle(dataframeToFillFrom:pd.DataFrame, IMFlatitude:float, IMFlongitude:float, datetime_to_run_across_UTC:dt.datetime):

    print("acquiring pitch angles...")
    pitch_angle_list = get_apply_method(dataframeToFillFrom)(lambda dataframe_row:get_pitch_angle_for_DF_analytic(IMFlatitude, IMFlongitude, dataframe_row["Lat"], dataframe_row["Long"]),axis=1)

    return pitch_angle_list

@numba.jit(nopython=True)
def get_pitch_angle_for_DF_analytic(IMFlatitude:float, IMFlongitude:float, asymptotic_dir_latitude:float, asymptotic_dir_longitude:float):

    IMFlatitude_rad = IMFlatitude * (np.pi/180.0)
    IMFlongitude_rad = IMFlongitude * (np.pi/180.0)
    asymptotic_dir_latitude_rad = asymptotic_dir_latitude * (np.pi/180.0)
    asymptotic_dir_longitude_rad = asymptotic_dir_longitude * (np.pi/180.0)

    cos_pitch_angle = (np.sin(asymptotic_dir_latitude_rad) * np.sin(IMFlatitude_rad)) + \
                      (np.cos(asymptotic_dir_latitude_rad) * np.cos(IMFlatitude_rad) * np.cos(asymptotic_dir_longitude_rad - IMFlongitude_rad))
    
    pitch_angle = np.arccos(cos_pitch_angle)

    return pitch_angle

def acquireWeightingFactors(asymptotic_direction_DF:pd.DataFrame, particle_dist:particleDistribution):

    momentaDist = particle_dist.momentum_distribution
    new_asymptotic_direction_DF = asymptotic_direction_DF.copy()

    # find weighting factors from the angles and rigidities
    #jacobian_function_to_use = lambda pitch_angle_in_radians:1/np.sin(2.0 * pitch_angle_in_radians)
        
    pitchAngleFunctionToUse = lambda row : momentaDist.getPitchAngleDistribution()(row["angleBetweenIMFinRadians"],row["Rigidity"])
    fullRigidityPitchWeightingFactorFunctionToUse = lambda row : momentaDist(row["angleBetweenIMFinRadians"],row["Rigidity"])

    print("calculating pitch angle weighting factors...")
    new_asymptotic_direction_DF["PitchAngleWeightingFactor"] = get_apply_method(new_asymptotic_direction_DF)(pitchAngleFunctionToUse,axis=1)

    new_asymptotic_direction_DF["Filter"] = (new_asymptotic_direction_DF["Filter"] == 1) * 1

    print("calculating rigidity weighting factors...")
    new_asymptotic_direction_DF["RigidityWeightingFactor"] = get_apply_method(new_asymptotic_direction_DF["Rigidity"])(momentaDist.getRigiditySpectrum())

    print("calculating rigidity + pitch combined weighting factors...")
    all_weighted_asymp_dirs = get_apply_method(new_asymptotic_direction_DF)(fullRigidityPitchWeightingFactorFunctionToUse,axis=1)
    new_asymptotic_direction_DF["fullRigidityPitchWeightingFactor"] = all_weighted_asymp_dirs * (new_asymptotic_direction_DF["Filter"] == 1)
    
    print("calculating energy + pitch combined weighting factors...")
    new_asymptotic_direction_DF["Energy"] = PRCT.convertParticleRigidityToEnergy(new_asymptotic_direction_DF["Rigidity"], 
                                                                      particleMassAU = particle_dist.particle_species.atomicMass, 
                                                                      particleChargeAU = particle_dist.particle_species.atomicNumber)
    energySpectrum = PRCT.convertParticleRigiditySpecToEnergySpec(new_asymptotic_direction_DF["Rigidity"],
                                                                  new_asymptotic_direction_DF["fullRigidityPitchWeightingFactor"],
                                                                  particleMassAU=particle_dist.particle_species.atomicMass, 
                                                                  particleChargeAU=particle_dist.particle_species.atomicNumber)

    new_asymptotic_direction_DF["fullEnergyPitchWeightingFactor"] = energySpectrum["Energy distribution values"]

    return new_asymptotic_direction_DF

def calculatePitchAngle_from_IMF_dir(interplanetary_mag_field:spaceCoords, asymptotic_direction:spaceCoords, datetime_to_run_across_UTC):

    cartesianAsympDir = asymptotic_direction.convert("GEO","car").data[0]
    GEOIMF = interplanetary_mag_field
    GEOIMF.ticks = spaceTicktock(datetime_to_run_across_UTC,"UTC")
    cartesianIMF = GEOIMF.convert("GEO","car").data[0]
    return calculateAngleBetweenTheSpaceVectors(cartesianAsympDir, cartesianIMF)

def calculatePitchAngle(momentaDist, dfRow, datetime_to_run_across_UTC):

    cartesianAsympDir = dfRow["Asymptotic Direction"].convert("GEO","car").data[0]
    GEOIMF = momentaDist.getPitchAngleDistribution().interplanetary_mag_field
    GEOIMF.ticks = spaceTicktock(datetime_to_run_across_UTC,"UTC")
    cartesianIMF = GEOIMF.convert("GEO","car").data[0]
    return calculateAngleBetweenTheSpaceVectors(cartesianAsympDir, cartesianIMF)

def calculateAngleBetweenTheSpaceVectors(cartesianAsympDir, cartesianIMF):
    dotProduct = np.dot(cartesianAsympDir, cartesianIMF)
    amplitude = np.linalg.norm(cartesianAsympDir)*np.linalg.norm(cartesianIMF)
    angleBetweenIMF = np.arccos(dotProduct/amplitude)
    return angleBetweenIMF
