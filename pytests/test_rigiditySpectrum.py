import pytest
import numpy as np
from AniMAIRE.anisotropic_MAIRE_engine.spectralCalculations.rigiditySpectrum import (
    rigiditySpectrum,
    powerLawSpectrum,
    interpolatedInputFileSpectrum,
    DLRmodelSpectrum,
    CommonModifiedPowerLawSpectrum,
    CommonModifiedPowerLawSpectrumSplit
)

def test_powerLawSpectrum():
    spectrum = powerLawSpectrum(normalisationFactor=1.0, spectralIndex=-2.7)
    assert np.isclose(spectrum(1.0), 1.0)
    assert np.isclose(spectrum(2.0), 2.0**(-2.7))

def test_interpolatedInputFileSpectrum(tmp_path):
    data = "1.0,10.0\n2.0,20.0\n3.0,30.0\n"
    file_path = tmp_path / "spectrum.csv"
    file_path.write_text(data)
    
    spectrum = interpolatedInputFileSpectrum(str(file_path))
    assert np.isclose(spectrum(1.0), 10.0 / (100**2))
    assert np.isclose(spectrum(2.0), 20.0 / (100**2))

def test_DLRmodelSpectrum():
    spectrum = DLRmodelSpectrum(atomicNumber=1, OULUcountRateInSeconds=100.0)
    assert spectrum(1.0) is not None

def test_CommonModifiedPowerLawSpectrum():
    spectrum = CommonModifiedPowerLawSpectrum(J0=1.0, gamma=2.7, deltaGamma=0.1)
    assert np.isclose(spectrum(1.0), 1.0 / (100**2))
    assert np.isclose(spectrum(2.0), 2.0**(-2.8) / 100**2)

def test_CommonModifiedPowerLawSpectrumSplit():
    spectrum = CommonModifiedPowerLawSpectrumSplit(J0=1.0, gamma=2.7, deltaGamma=0.1)
    assert np.isclose(spectrum(1.0), 1.0 / (100**2))
    assert np.isclose(spectrum(2.0), 2.0**(-2.8) / 100**2)
