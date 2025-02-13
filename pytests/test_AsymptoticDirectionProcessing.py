import pytest
import pandas as pd
import numpy as np
import datetime as dt
from AniMAIRE.anisotropic_MAIRE_engine.AsymptoticDirectionProcessing import (
    generate_asymp_dir_DF,
    acquireWeightingFactors,
    get_pitch_angle_for_DF_analytic
)
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution

@pytest.fixture
def sample_dataframe():
    data = {
        "Rigidity": [1.0, 2.0, 3.0],
        "Lat": [10.0, 20.0, 30.0],
        "Long": [100.0, 110.0, 120.0],
        "Filter": [1, 1, 1]
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_particle_distribution():
    # Create a mock particle distribution object with correct parameters
    class MockRigiditySpectrum:
        def __call__(self, rigidity):
            return rigidity

    class MockPitchAngleDistribution:
        def getPitchAngleDistribution(self, angle, rigidity):
            return angle * rigidity

        def __call__(self, angle, rigidity):
            return angle * rigidity

    return particleDistribution(
        particle_name="proton",
        rigidity_spectrum=MockRigiditySpectrum(),
        pitch_angle_distribution=MockPitchAngleDistribution()
    )

def test_generate_asymp_dir_DF(sample_dataframe):
    IMFlatitude = 0.0
    IMFlongitude = 0.0
    datetime_to_run_across_UTC = dt.datetime.utcnow()
    cache = False

    result = generate_asymp_dir_DF(sample_dataframe, IMFlatitude, IMFlongitude, datetime_to_run_across_UTC, cache)
    assert "angleBetweenIMFinRadians" in result.columns

def test_acquireWeightingFactors(sample_dataframe, sample_particle_distribution):
    result = acquireWeightingFactors(sample_dataframe, sample_particle_distribution)
    assert "PitchAngleWeightingFactor" in result.columns
    assert "RigidityWeightingFactor" in result.columns
    assert "fullRigidityPitchWeightingFactor" in result.columns
    assert "fullEnergyPitchWeightingFactor" in result.columns

def test_get_pitch_angle_for_DF_analytic():
    IMFlatitude = 0.0
    IMFlongitude = 0.0
    asymptotic_dir_latitude = 10.0
    asymptotic_dir_longitude = 100.0

    result = get_pitch_angle_for_DF_analytic(IMFlatitude, IMFlongitude, asymptotic_dir_latitude, asymptotic_dir_longitude)
    assert isinstance(result, float)
