import pytest
import pandas as pd
import numpy as np
from AniMAIRE.anisotropic_MAIRE_engine.singleParticleEngineInstance import singleParticleEngineInstance
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import rigiditySpectrum
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import isotropicPitchAngleDistribution

@pytest.fixture
def sample_particle_distribution():
    spectrum = rigiditySpectrum()
    pitch_angle_dist = isotropicPitchAngleDistribution()
    return particleDistribution("proton", spectrum, pitch_angle_dist)

@pytest.fixture
def sample_asymptotic_directions():
    data = {
        "initialLatitude": [0.0, 10.0],
        "initialLongitude": [0.0, 10.0],
        "Rigidity": [1.0, 2.0],
        "fullRigidityPitchWeightingFactor": [0.5, 0.5]
    }
    return pd.DataFrame(data)

def test_singleParticleEngineInstance(sample_particle_distribution, sample_asymptotic_directions):
    engine = singleParticleEngineInstance(
        particle_distribution=sample_particle_distribution,
        dfofAsymptoticDirections=sample_asymptotic_directions,
        list_of_altitudes_in_km=[0.0, 10.0, 20.0],
        generate_NM_count_rates=False
    )
    
    assert engine.particle_distribution == sample_particle_distribution
    assert engine.list_of_altitudes_in_km == [0.0, 10.0, 20.0]
    assert engine.generate_NM_count_rates == False
