import pytest
import datetime as dt
import numpy as np
from AniMAIRE.anisotropic_MAIRE_engine.generalEngineInstance import generalEngineInstance
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import rigiditySpectrum
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import isotropicPitchAngleDistribution

@pytest.fixture
def sample_particle_distribution():
    spectrum = rigiditySpectrum()
    pitch_angle_dist = isotropicPitchAngleDistribution()
    return particleDistribution("proton", spectrum, pitch_angle_dist)

def test_generalEngineInstance(sample_particle_distribution):
    engine = generalEngineInstance(
        list_of_particle_distributions=[sample_particle_distribution],
        list_of_altitudes_km=[0.0, 10.0, 20.0],
        Kp_index=3,
        date_and_time=dt.datetime(2023, 1, 1),
        reference_latitude=0.0,
        reference_longitude=45.0,
        array_of_lats_and_longs=np.array([[0.0, 0.0], [10.0, 10.0]]),
        cache_magnetocosmics_runs=False,
        generate_NM_count_rates=False
    )
    
    assert engine.Kp_index == 3
    assert engine.date_and_time == dt.datetime(2023, 1, 1)
    assert engine.reference_latitude == 0.0
    assert engine.reference_longitude == 45.0
    assert engine.cache_magnetocosmics_runs == False
    assert engine.generate_NM_count_rates == False
