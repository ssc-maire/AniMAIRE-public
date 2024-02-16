import numpy as np
import pandas as pd
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

memory_asymp_dirs = Memory("./cacheAsymptoticDirectionOutputs", verbose=0)

m0 = 1.67262192e-27 #kg
c = 299792458.0 #m/s
protonCharge = 1.60217663e-19 #C
protonRestEnergy = m0 * (c**2)

gettrace = getattr(sys, 'gettrace', None)
global get_apply_method
if gettrace is None:
    #print('No sys.gettrace')
    get_apply_method = lambda DF_or_Series:DF_or_Series.parallel_apply
else:
    get_apply_method = lambda DF_or_Series:DF_or_Series.progress_apply
    

def generate_asymp_dir_DF(dataframeToFillFrom:pd.DataFrame, IMFlatitude:float, IMFlongitude:float, datetime_to_run_across_UTC, cache:bool):

    new_asymp_dir_DF = dataframeToFillFrom.copy()

    new_asymp_dir_DF["Energy"] = PRCT.convertParticleRigidityToEnergy(dataframeToFillFrom["Rigidity"], particleMassAU = 1, particleChargeAU = 1)

    print("assigning asymptotic coordinates")
    if cache == False:
        asymptoticDirectionList = convertAsymptoticDirectionsToPitchAngle(dataframeToFillFrom, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC)
    else:
        cachedConvertAsympDirFunc = memory_asymp_dirs.cache(convertAsymptoticDirectionsToPitchAngle)
        asymptoticDirectionList = cachedConvertAsympDirFunc(dataframeToFillFrom, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC)

    print("successfully converted asymptotic directions")

    new_asymp_dir_DF["angleBetweenIMFinRadians"] = asymptoticDirectionList

    return new_asymp_dir_DF

#@memory_asymp_dirs.cache()
def convertInitialCoordSetToSpacePy(dataframeRow:pd.Series, radialPosition:float, firstCoordinateLabel:str, secondCoordinateLabel:str, datetime_to_run_across_UTC):

    rowLatitude = dataframeRow[firstCoordinateLabel]
    rowLongitude = dataframeRow[secondCoordinateLabel]

    rowInSpacepyCoords = spaceCoords([radialPosition,rowLatitude,rowLongitude],
                                     "GEO","sph",
                                     units=['Re', 'deg', 'deg'],
                                     use_irbem=False)
    rowInSpacepyCoords.ticks = spaceTicktock(datetime_to_run_across_UTC,"UTC")

    return rowInSpacepyCoords

def convertAsymptoticDirectionsToPitchAngle(dataframeToFillFrom, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC):
    
    # pitch_angle_list = []
    # for index, dataframeRow in tqdm.tqdm(dataframeToFillFrom.iterrows(),total=len(dataframeToFillFrom)):
    #     pitch_angle_for_row = get_pitch_angle_for_DF_row(IMFlatitude, IMFlongitude, dataframeRow)
    #     pitch_angle_list.append(pitch_angle_for_row)

    print("acquiring pitch angles...")
    #print(dataframeToFillFrom)
    #pitch_angle_list = dataframeToFillFrom.parallel_apply(lambda dataframe_row:get_pitch_angle_for_DF_row(IMFlatitude, IMFlongitude, dataframe_row, datetime_to_run_across_UTC),axis=1)
    
    #pitch_angle_list = dataframeToFillFrom.parallel_apply(lambda dataframe_row:get_pitch_angle_for_DF_analytic(IMFlatitude, IMFlongitude, dataframe_row["Lat"], dataframe_row["Long"]),axis=1)
    pitch_angle_list = dataframeToFillFrom.progress_apply(lambda dataframe_row:get_pitch_angle_for_DF_analytic(IMFlatitude, IMFlongitude, dataframe_row["Lat"], dataframe_row["Long"]),axis=1)

    #pitch_angle_list = dataframeToFillFrom.parallel_apply(lambda dataframe_row:get_pitch_angle_for_DF_Mishev(IMFlatitude, IMFlongitude, dataframe_row["Lat"], dataframe_row["Long"]),axis=1)

    #print(pitch_angle_list)

    #asymptoticDirectionList = list(dataframeToFillFrom.progress_apply(lambda dataframeRow:convertInitialCoordSetToSpacePy(dataframeRow[1], 100.0,"Lat","Long")))

    # with pmp.Pool(processes = 4) as pmpPool:
    #     asymptoticDirectionList = pmpPool.map(lambda dataframeRow:convertInitialCoordSetToSpacePy(dataframeRow[1], 100.0,"Lat","Long"),
    #                                         tqdm.tqdm(list(dataframeToFillFrom.iterrows()),total=len(dataframeToFillFrom)))

    return pitch_angle_list

@numba.jit(nopython=True)
def get_pitch_angle_for_DF_analytic(IMFlatitude, IMFlongitude, asymptotic_dir_latitude, asymptotic_dir_longitude):

    IMFlatitude_rad = IMFlatitude * (np.pi/180.0)
    IMFlongitude_rad = IMFlongitude * (np.pi/180.0)
    asymptotic_dir_latitude_rad = asymptotic_dir_latitude * (np.pi/180.0)
    asymptotic_dir_longitude_rad = asymptotic_dir_longitude * (np.pi/180.0)

    cos_pitch_angle = (np.sin(asymptotic_dir_latitude_rad) * np.sin(IMFlatitude_rad)) + \
                      (np.cos(asymptotic_dir_latitude_rad) * np.cos(IMFlatitude_rad) * np.cos(asymptotic_dir_longitude_rad - IMFlongitude_rad))
    
    pitch_angle = np.arccos(cos_pitch_angle)

    return pitch_angle

def get_pitch_angle_for_DF_Mishev(IMFlatitude, IMFlongitude, asymptotic_dir_latitude, asymptotic_dir_longitude):

    IMFlatitude_rad = IMFlatitude * (np.pi/180.0)
    IMFlongitude_rad = IMFlongitude * (np.pi/180.0)
    asymptotic_dir_latitude_rad = asymptotic_dir_latitude * (np.pi/180.0)
    asymptotic_dir_longitude_rad = asymptotic_dir_longitude * (np.pi/180.0)

    cos_pitch_angle = (np.sin(asymptotic_dir_latitude_rad) * np.sin(IMFlatitude_rad)) + \
                      (np.cos(asymptotic_dir_latitude_rad) * np.cos(IMFlatitude_rad) * np.cos(IMFlongitude_rad) * (np.cos(asymptotic_dir_longitude_rad) + np.sin(asymptotic_dir_longitude_rad)))
    
    pitch_angle = np.arccos(cos_pitch_angle)

    return pitch_angle

def get_pitch_angle_for_DF_row(IMFlatitude, IMFlongitude, dataframeRow, datetime_to_run_across_UTC):
    #interplanetary_mag_field = spaceCoords([100.0,IMFlatitude, IMFlongitude],"GSM","sph",units=['Re', 'deg', 'deg'],use_irbem=False)
    interplanetary_mag_field = spaceCoords([100.0,IMFlatitude, IMFlongitude],"GEO","sph",units=['Re', 'deg', 'deg'],use_irbem=False)
    rowInSpacepyCoords = convertInitialCoordSetToSpacePy(dataframeRow, 100.0,"Lat","Long", datetime_to_run_across_UTC)
    pitch_angle_for_row = calculatePitchAngle_from_IMF_dir(interplanetary_mag_field, rowInSpacepyCoords, datetime_to_run_across_UTC)
    return pitch_angle_for_row

def acquireWeightingFactors(asymptotic_direction_DF:pd.DataFrame, momentaDist):

    new_asymptotic_direction_DF = asymptotic_direction_DF.copy()

    # # find angle between asymptotic directions and IMF using dot product
    # print("determining angles between asympotic directions and IMF...")
    # angleColumnTitle = "angleBetweenIMFinRadians"
    # new_asymptotic_direction_DF = assignPitchAngles(new_asymptotic_direction_DF, momentaDist, angleColumnTitle)

    # find weighting factors from the angles and rigidities
    #jacobian_function_to_use = lambda pitch_angle_in_radians:1/np.sin(2.0 * pitch_angle_in_radians)
        
    pitchAngleFunctionToUse = lambda row : momentaDist.getPitchAngleDistribution()(row["angleBetweenIMFinRadians"],row["Rigidity"])
    fullRigidityPitchWeightingFactorFunctionToUse = lambda row : momentaDist(row["angleBetweenIMFinRadians"],row["Rigidity"])

    #print(new_asymptotic_direction_DF)
    print("calculating pitch angle weighting factors...")
    #new_asymptotic_direction_DF["PitchAngleWeightingFactor"] = new_asymptotic_direction_DF.progress_apply(pitchAngleFunctionToUse,axis=1)
    #new_asymptotic_direction_DF["PitchAngleWeightingFactor"] = new_asymptotic_direction_DF.parallel_apply(pitchAngleFunctionToUse,axis=1)
    new_asymptotic_direction_DF["PitchAngleWeightingFactor"] = get_apply_method(new_asymptotic_direction_DF)(pitchAngleFunctionToUse,axis=1)
    #print("pitch angle weights:")
    #print(new_asymptotic_direction_DF[["angleBetweenIMFinRadians","PitchAngleWeightingFactor"]])
    #new_asymptotic_direction_DF["PitchAngleWeightingFactor"] = 1
    #print(new_asymptotic_direction_DF)

    #new_asymptotic_direction_DF["Filter"] = new_asymptotic_direction_DF["Filter"].loc[::-1].cummax().loc[::-1]
    #new_asymptotic_direction_DF = new_asymptotic_direction_DF[new_asymptotic_direction_DF["Filter"] == 1]
    new_asymptotic_direction_DF["Filter"] = (new_asymptotic_direction_DF["Filter"] == 1) * 1

    print("calculating rigidity weighting factors...")
    new_asymptotic_direction_DF["RigidityWeightingFactor"] = get_apply_method(new_asymptotic_direction_DF["Rigidity"])(momentaDist.getRigiditySpectrum())
    #print(new_asymptotic_direction_DF)

    print("calculating rigidity + pitch combined weighting factors...")
    new_asymptotic_direction_DF["fullRigidityPitchWeightingFactor"] = get_apply_method(new_asymptotic_direction_DF)(fullRigidityPitchWeightingFactorFunctionToUse,axis=1) * (new_asymptotic_direction_DF["Filter"] == 1)
    #print(new_asymptotic_direction_DF)

    print("calculating energy + pitch combined weighting factors...")
    energySpectrum = PRCT.convertParticleRigiditySpecToEnergySpec(new_asymptotic_direction_DF["Rigidity"],
                                                                  new_asymptotic_direction_DF["fullRigidityPitchWeightingFactor"],
                                                                  particleMassAU=1, particleChargeAU=1)

    new_asymptotic_direction_DF["fullEnergyPitchWeightingFactor"] = energySpectrum["Energy distribution values"]
    #print(new_asymptotic_direction_DF)

    return new_asymptotic_direction_DF

def assignPitchAngles(asymp_dir_DF:pd.DataFrame, momentaDist, angleColumnTitle):

    new_asymp_dir_DF = asymp_dir_DF.copy()
    
    # with pmp.Pool(processes = 8) as pmpPool:
    #     listOfOutputDoseDFs = pmpPool.map(lambda dfRow:self.calculatePitchAngle(momentaDist, dfRow),
    #                                         tqdm.tqdm(self.iterrows(),total=len(self)))

    angleArray = []
    print("WARNING: this loop could become slow for large dataframes. Optimising this dataframe could be done through pandas array manipulation or through numba or cpython")
    for index, dfRow in tqdm.tqdm(new_asymp_dir_DF.iterrows(),total=len(new_asymp_dir_DF)):
        angleBetweenIMF = calculatePitchAngle(momentaDist, dfRow)
        angleArray.append(angleBetweenIMF)
    new_asymp_dir_DF[angleColumnTitle] = angleArray

    #self[angleColumnTitle] = self.to_DataFrame()["Asymptotic Direction"].apply(lambda x:self.calculatePitchAngleApply(momentaDist,x))
    
    #pitchAngleDF = self.to_DataFrame()["Asymptotic Direction"].progress_apply(lambda x:self.calculatePitchAngleApply(momentaDist,x))
    #self[angleColumnTitle] = pitchAngleDF #for some reason I believe a dummy variable is REQUIRED - otherwise the series does not get set correctly. However for some reason this function is still not working.
    
    #self[angleColumnTitle] = self.to_DataFrame()["Asymptotic Direction"].parallel_apply(lambda x:self.calculatePitchAngleApply(momentaDist,x))

    return new_asymp_dir_DF

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

#@staticmethod
#@jit(nopython=True)
def calculateAngleBetweenTheSpaceVectors(cartesianAsympDir, cartesianIMF):
    dotProduct = np.dot(cartesianAsympDir, cartesianIMF)
    amplitude = np.linalg.norm(cartesianAsympDir)*np.linalg.norm(cartesianIMF)
    angleBetweenIMF = np.arccos(dotProduct/amplitude)
    return angleBetweenIMF

def splitDataframeIntoUniqueCoordinates(self):

    columnsToSplitBy = ["initialLatitude","initialLongitude","initialZenith","initialAzimuth"]

    vanillaSelfDataFrame = self.to_DataFrame()

    uniqueRelevantRowsDF = vanillaSelfDataFrame[columnsToSplitBy].drop_duplicates()
    
    outputListOfSplitDFs = []
    for index, row in tqdm.tqdm(uniqueRelevantRowsDF.iterrows(),total=len(uniqueRelevantRowsDF)):
        useRowBoolean = True
        for column in columnsToSplitBy:
            useRowBoolean = (vanillaSelfDataFrame[column] == row[column]) & useRowBoolean
        
        outputListOfSplitDFs.append(vanillaSelfDataFrame[useRowBoolean])
    
    return outputListOfSplitDFs
