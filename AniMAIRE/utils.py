from .anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution
from .anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import isotropicPitchAngleDistribution, pitchAngleDistribution
import numpy as np
import datetime as dt
import spaceweather as sw
import pandas as pd
from typing import List, Optional

def get_correctly_formatted_particle_dist_list(proton_rigidity_spectrum, 
                                               alpha_rigidity_spectrum, 
                                               reference_pitch_angle_latitude, 
                                               reference_pitch_angle_longitude, 
                                               proton_pitch_angle_distribution, 
                                               alpha_pitch_angle_distribution):
    list_of_particle_distributions = []
    for particle_name, particle_rigidity_spectrum, particle_pitch_angle_distribution in [["proton",proton_rigidity_spectrum,proton_pitch_angle_distribution], 
                                                                                         ["alpha",alpha_rigidity_spectrum,alpha_pitch_angle_distribution]]:
        if not particle_rigidity_spectrum is None:
            if not type(particle_pitch_angle_distribution).__name__ == 'isotropicPitchAngleDistribution':
                if (reference_pitch_angle_latitude is None) or (reference_pitch_angle_longitude is None):
                    raise Exception("Error: a non-isotropic pitch angle distribution was used but reference longitudes and latitudes were not set!")

                particle_pitch_angle_distribution = pitchAngleDistribution(particle_pitch_angle_distribution,
                                                                  reference_pitch_angle_latitude,
                                                                  reference_pitch_angle_longitude)
                
            particle_distribution = particleDistribution(
                                                particle_name,
                                                particle_rigidity_spectrum,
                                                particle_pitch_angle_distribution)
            
            list_of_particle_distributions.append(particle_distribution)
    return list_of_particle_distributions

def get_kp_index(date_and_time: dt.datetime) -> int:
    return sw.ap_kp_3h()[pd.to_datetime(sw.ap_kp_3h().index) < date_and_time].iloc[-1]["Kp"]

default_altitudes_in_kft = [0,10,20] + [i for i in range(25, 61 + 1, 3)]

def validate_altitudes(altitudes_in_km: Optional[List[float]], altitudes_in_kft: List[float]) -> List[float]:
    if (altitudes_in_km is not None) and (altitudes_in_kft is not None):
        raise Exception("Error: only one of altitudes_in_km and altitudes_in_kft should be supplied!")
    elif (altitudes_in_km is None) and (altitudes_in_kft is None):
        return np.array(default_altitudes_in_kft) * 0.3048
    elif altitudes_in_km is not None:
        return altitudes_in_km
    elif altitudes_in_kft is not None:
        return np.array(altitudes_in_kft) * 0.3048 if altitudes_in_km is None else altitudes_in_km

