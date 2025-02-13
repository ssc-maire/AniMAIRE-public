import pytest
import numpy as np
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import (
    pitchAngleDistribution,
    cosinePitchAngleDistribution,
    isotropicPitchAngleDistribution,
    gaussianPitchAngleDistribution,
    gaussianBeeckPitchAngleDistribution
)

def test_pitchAngleDistribution():
    dist = pitchAngleDistribution(lambda x, y: x * y, 0.0, 0.0)
    assert np.isclose(dist(1.0, 2.0), 2.0)

def test_cosinePitchAngleDistribution():
    dist = cosinePitchAngleDistribution()
    assert np.isclose(dist(0.0, 1.0), 0.0)
    assert np.isclose(dist(np.pi / 2, 1.0), 0.5)

def test_isotropicPitchAngleDistribution():
    dist = isotropicPitchAngleDistribution()
    assert np.isclose(dist(0.0, 1.0), 1.0)
    assert np.isclose(dist(np.pi / 2, 1.0), 1.0)

def test_gaussianPitchAngleDistribution():
    dist = gaussianPitchAngleDistribution(normFactor=1.0, sigma=1.0)
    assert np.isclose(dist(0.0, 1.0), 1.0)
    assert np.isclose(dist(1.0, 1.0), np.exp(-0.5))

def test_gaussianBeeckPitchAngleDistribution():
    dist = gaussianBeeckPitchAngleDistribution(normFactor=1.0, A=1.0, B=1.0)
    assert np.isclose(dist(0.0, 1.0), 1.0)
    assert np.isclose(dist(np.pi / 2, 1.0), np.exp(-0.5))
