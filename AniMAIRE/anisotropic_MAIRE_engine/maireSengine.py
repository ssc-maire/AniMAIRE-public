#!/bin/python3
import numpy as np
import pandas as pd
from tqdm import tqdm

from .singleParticleEngineInstance import singleParticleEngineInstance

tqdm.pandas()


from AsympDirsCalculator import AsympDirsTools

from .utils.AsymptoticDirectionDataframe import generate_asymp_dir_DF

import datetime as dt

from joblib import Memory
MAGCOScachedir = 'cachedMagnetocosmicsRunData'
MAGCOSmemory = Memory(MAGCOScachedir, verbose=0)
default_array_of_lats_and_longs = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

class generalEngineInstance():

    def __init__(self, 
                 list_of_particle_distributions:list,
                 list_of_altitudes_km:list, 
                 Kp_index:int, 
                 date_and_time:dt.datetime,
                 reference_latitude=0.0, #None,
                 reference_longitude=45.0, #None,
                 array_of_lats_and_longs=default_array_of_lats_and_longs,
                 cache_magnetocosmics_runs=True,
                 generate_NM_count_rates=False):
        
        self.rigiditySpectrumParamDict = {}
        self.pitchAngleDistributionParamDict = {}

        self.list_of_particle_distributions = list_of_particle_distributions
        self.list_of_altitudes_km = list_of_altitudes_km
        self.Kp_index = Kp_index
        self.date_and_time = date_and_time
        self.reference_latitude = reference_latitude
        self.reference_longitude = reference_longitude
        self.array_of_lats_and_longs = array_of_lats_and_longs

        self.cache_magnetocosmics_runs = cache_magnetocosmics_runs
        self.generate_NM_count_rates = generate_NM_count_rates

    def getAsymptoticDirsAndRun(self,**mag_cos_kwargs)->pd.DataFrame:
        
        self.acquireDFofAllAsymptoticDirections(**mag_cos_kwargs)

        fullDoseRateList = []

        for incoming_particle_distribution in self.list_of_particle_distributions:
            singleParticleEngine = singleParticleEngineInstance(incoming_particle_distribution, 
                                                                self.df_of_asymptotic_directions,
                                                                self.list_of_altitudes_km,
                                                                self.generate_NM_count_rates)
            
            doseRateDFforParticleSpecies = singleParticleEngine.runOverSpecifiedAltitudes()

            fullDoseRateList.append(doseRateDFforParticleSpecies)

        summedDoseRateDF = fullDoseRateList[0]
        for doseRateDF in fullDoseRateList[1:]:
            for doseRateName in ["adose","edose","dosee","SEU","SEL"]:
                summedDoseRateDF[doseRateName] += doseRateDF[doseRateName]

        return summedDoseRateDF
    
    def acquireDFofAllAsymptoticDirections(self, **mag_cos_kwargs):
        # run magnetocosmics to get asymptotic directions
        raw_asymp_dir_DF = AsympDirsTools.get_magcos_asymp_dirs(
                                                array_of_lats_and_longs=self.array_of_lats_and_longs,
                                                KpIndex=self.Kp_index,
                                                dateAndTime=self.date_and_time,
                                                cache=self.cache_magnetocosmics_runs,
                                                full_output=True,
                                                **mag_cos_kwargs,
                                                )
        raw_asymp_dir_DF.to_pickle("raw_asymp_dir_DF.pkl")
        self.df_of_asymptotic_directions = generate_asymp_dir_DF(raw_asymp_dir_DF, 
                                                                 self.reference_latitude, 
                                                                 self.reference_longitude, 
                                                                 self.date_and_time,
                                                                 cache = False)
        self.df_of_asymptotic_directions.to_pickle("self_df_of_asymptotic_directions.pkl")





