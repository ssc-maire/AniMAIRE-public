#!/bin/python3
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Memory
import datetime as dt
from typing import Optional

from .singleParticleEngineInstance import singleParticleEngineInstance
from AsympDirsCalculator import AsympDirsTools
from .AsymptoticDirectionProcessing import generate_asymp_dir_DF
import os

# Initialize tqdm for progress bars
tqdm.pandas()

# Set up caching for Magnetocosmics run data
MAGCOScachedir = 'cachedMagnetocosmicsRunData'
MAGCOSmemory = Memory(MAGCOScachedir, verbose=0)

# Default array of latitudes and longitudes
default_array_of_lats_and_longs = np.array(np.meshgrid(np.linspace(-90.0, 90.0, 37), np.linspace(0.0, 355.0, 72))).T.reshape(-1, 2)

class generalEngineInstance:
    """
    General engine instance for running dose rate calculations.
    """

    def __init__(self, 
                 list_of_particle_distributions: list,
                 list_of_altitudes_km: list[float], 
                 Kp_index: int, 
                 date_and_time: dt.datetime,
                 reference_latitude: float = 0.0,
                 reference_longitude: float = 45.0,
                 array_of_lats_and_longs: np.ndarray = default_array_of_lats_and_longs,
                 cache_magnetocosmics_runs: bool = True,
                 generate_NM_count_rates: bool = False,
                 asymp_dir_file: Optional[str] = None):
        """
        Initialize the general engine instance with necessary parameters.

        Parameters:
        - list_of_particle_distributions: list
            List of particle distributions.
        - list_of_altitudes_km: list[float]
            List of altitudes in kilometers.
        - Kp_index: int
            Kp index for the calculations.
        - date_and_time: dt.datetime
            Date and time for the calculations.
        - reference_latitude: float, optional
            Reference latitude for pitch angle distribution.
        - reference_longitude: float, optional
            Reference longitude for pitch angle distribution.
        - array_of_lats_and_longs: np.ndarray, optional
            Array of latitudes and longitudes.
        - cache_magnetocosmics_runs: bool, optional
            Whether to cache Magnetocosmics runs.
        - generate_NM_count_rates: bool, optional
            Whether to generate neutron monitor count rates.
        - asymp_dir_file: str, optional
            File path to read asymptotic directions from.
        """
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
        self.asymp_dir_file = asymp_dir_file

    def getAsymptoticDirsAndRun(self, use_default_9_zeniths_azimuths: bool, **mag_cos_kwargs) -> pd.DataFrame:
        """
        Acquire asymptotic directions and run calculations.

        Parameters:
        - use_default_9_zeniths_azimuths: bool
            Whether to use default 9 zeniths and azimuths.
        - **mag_cos_kwargs: additional keyword arguments
            Additional arguments for Magnetocosmics.

        Returns:
        - pd.DataFrame
            DataFrame containing the dose rate calculations.
        """
        self.acquireDFofAllAsymptoticDirections(use_default_9_zeniths_azimuths, **mag_cos_kwargs)

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
            for doseRateName in ["adose", "edose", "dosee", "SEU", "SEL"]:
                summedDoseRateDF[doseRateName] += doseRateDF[doseRateName]

        return summedDoseRateDF
    
    def acquireDFofAllAsymptoticDirections(self, use_default_9_zeniths_and_azimuths: bool, **magneto_kwargs):
        """
        Acquire DataFrame of all asymptotic directions.

        Parameters:
        - use_default_9_zeniths_and_azimuths: bool
            Whether to use the default 9 zeniths and azimuths.
        - **magneto_kwargs: additional keyword arguments for Magnetocosmics.
        """
        if self.asymp_dir_file:
            raw_asymp_df = self.get_raw_asymp_DF_from_file(self.asymp_dir_file)
        else:
            if use_default_9_zeniths_and_azimuths and "array_of_zeniths_and_azimuths" in magneto_kwargs:
                raise Exception("Error: Both use_default_9_zeniths_and_azimuths is True and 'array_of_zeniths_and_azimuths' is specified.")
            
            if use_default_9_zeniths_and_azimuths:
                default_zeniths_azimuths = [
                    [0.0, 0.0],
                    [16.0, 0.0],
                    [16.0, 90.0],
                    [16.0, 180.0],
                    [16.0, 270.0],
                    [32.0, 0.0],
                    [32.0, 90.0],
                    [32.0, 180.0],
                    [32.0, 270.0],
                ]
                raw_asymp_df = AsympDirsTools.get_magcos_asymp_dirs(
                    array_of_lats_and_longs=self.array_of_lats_and_longs,
                    KpIndex=self.Kp_index,
                    dateAndTime=self.date_and_time,
                    cache=self.cache_magnetocosmics_runs,
                    full_output=True,
                    array_of_zeniths_and_azimuths=default_zeniths_azimuths,
                    **magneto_kwargs,
                )
            else:
                raw_asymp_df = AsympDirsTools.get_magcos_asymp_dirs(
                    array_of_lats_and_longs=self.array_of_lats_and_longs,
                    KpIndex=self.Kp_index,
                    dateAndTime=self.date_and_time,
                    cache=self.cache_magnetocosmics_runs,
                    full_output=True,
                    **magneto_kwargs,
                )
                
        raw_asymp_df.to_pickle("raw_asymp_dir_DF.pkl")
        processed_df = generate_asymp_dir_DF(
            raw_asymp_df,
            self.reference_latitude,
            self.reference_longitude,
            self.date_and_time,
            cache=False
        )
        processed_df.to_csv("self_df_of_asymptotic_directions.csv", index=False)
        self.df_of_asymptotic_directions = processed_df

    def get_raw_asymp_DF_from_file(self,file_path):
        if isinstance(file_path, list):
            raw_dfs = []
            for file_path in file_path:
                init_lat, init_lon = self._parse_initial_coordinates(file_path)
                df = pd.read_csv(file_path, skipfooter=1, engine='python')
                df["initialLatitude"] = init_lat
                df["initialLongitude"] = init_lon
                raw_dfs.append(self.validate_asymp_dir_df(df))
            raw_asymp_df = pd.concat(raw_dfs, ignore_index=True)
        else:
            init_lat, init_lon = self._parse_initial_coordinates(file_path)
            raw_df_from_file = pd.read_csv(file_path, skipfooter=1, engine='python')
            raw_df_from_file["initialLatitude"] = init_lat
            raw_df_from_file["initialLongitude"] = init_lon
            raw_asymp_df = self.validate_asymp_dir_df(raw_df_from_file)

        return raw_asymp_df

    def _parse_initial_coordinates(self, file_path: str) -> tuple:
        """
        Parse initial latitude and longitude from the filename.

        Expects filename format: "latitude_longitude*.csv" where latitude and longitude are numeric.

        Parameters:
        - file_path: str
            Path to the asymptotic directions file.

        Returns:
        - tuple: (initial latitude, initial longitude)
        """
        file_name = os.path.basename(file_path)
        name_without_ext, _ = os.path.splitext(file_name)
        try:
            latitude_str, longitude_str = name_without_ext.split('_')
            initial_lat = float(latitude_str)
            initial_lon = float(longitude_str)
        except Exception as error:
            raise ValueError("Filename must be in format 'latitude_longitude*.csv' with numeric latitude and longitude.") from error
        return initial_lat, initial_lon

    def validate_asymp_dir_df(self, df: pd.DataFrame):
        """
        Validate the structure of the asymptotic directions DataFrame.

        Parameters:
        - df: pd.DataFrame
            DataFrame to validate.
        """
        required_columns_set1 = [
            "initialLatitude", "initialLongitude", "Rigidity", "Lat", "Long", "Filter"
        ]
        required_columns_set2 = [
            "Rigidity(GV)", "Filter", "Latitude", "Longitude", "X", "Y", "Z"
        ]
        
        if all(col in df.columns for col in required_columns_set1):
            print("Asymptotic directions DataFrame validated successfully with the first format.")
            return df
        elif all(col in df.columns for col in required_columns_set2):
            print("Asymptotic directions DataFrame validated successfully with the second format.")
            transformed_df = pd.DataFrame({
                "initialLatitude": df["initialLatitude"],
                "initialLongitude": df["initialLongitude"],
                "Rigidity": df["Rigidity(GV)"],
                "Lat": df["Latitude"],
                "Long": df["Longitude"],
                "Filter": df["Filter"]
            })
            print("Asymptotic directions DataFrame converted from second format to first format with X and Y converted to new latitudes and longitudes.")
            return transformed_df
        else:
            raise ValueError("Missing required columns for both formats.")





