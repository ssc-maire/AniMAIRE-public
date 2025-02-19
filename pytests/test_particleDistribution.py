import pytest
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import rigiditySpectrum
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import isotropicPitchAngleDistribution

def test_particleDistribution():
    spectrum = rigiditySpectrum()
    pitch_angle_dist = isotropicPitchAngleDistribution()
    dist = particleDistribution("proton", spectrum, pitch_angle_dist)
    
    assert dist.particle_species.particleName == "proton"
    assert dist.momentum_distribution.rigidity_spectrum == spectrum
    assert dist.momentum_distribution.pitch_angle_distribution == pitch_angle_dist
