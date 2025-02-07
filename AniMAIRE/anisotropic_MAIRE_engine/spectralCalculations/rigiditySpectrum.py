from scipy.interpolate import interp1d
import pandas as pd
import scipy
import pickle as pkl
import numpy as np
import datetime

from CosRayModifiedISO import CosRayModifiedISO

print("WARNING: currently unknown whether the reported spectral weighting factor should be in terms of energy or rigidity")

class rigiditySpectrum():
    """
    Base class for rigidity spectra.
    """

    def __init__(self):
        """
        Initialize the rigidity spectrum.
        """
        pass

    def __call__(self, x: float) -> float:
        """
        Evaluate the rigidity spectrum at a given rigidity.

        Parameters:
        - x: float
            The rigidity.

        Returns:
        - float
            The value of the rigidity spectrum.
        """
        return self.rigiditySpec(x)
    
    def __add__(self, right: 'rigiditySpectrum') -> 'rigiditySpectrum':
        """
        Add two rigidity spectra.

        Parameters:
        - right: rigiditySpectrum
            The other rigidity spectrum.

        Returns:
        - rigiditySpectrum
            The sum of the two rigidity spectra.
        """
        summed_spectrum = rigiditySpectrum()
        summed_spectrum.rigiditySpec = lambda x: self.rigiditySpec(x) + right.rigiditySpec(x)
        return summed_spectrum

class powerLawSpectrum(rigiditySpectrum):
    """
    Power law rigidity spectrum.
    """

    rigiditySpec = lambda self, x: self.normalisationFactor * (x ** self.spectralIndex)

    def __init__(self, normalisationFactor: float, spectralIndex: float):
        """
        Initialize the power law spectrum.

        Parameters:
        - normalisationFactor: float
            The normalization factor.
        - spectralIndex: float
            The spectral index.
        """
        self.normalisationFactor = normalisationFactor
        self.spectralIndex = spectralIndex

class interpolatedInputFileSpectrum(rigiditySpectrum):
    """
    Interpolated rigidity spectrum from an input file.
    """

    def __init__(self, inputFileName: str):
        """
        Initialize the interpolated spectrum.

        Parameters:
        - inputFileName: str
            The path to the input file.
        """
        self.inputFilename = inputFileName
        self.rigiditySpec = self.readSpecFromCSV(self.inputFilename)

    def readSpecFromCSV(self, inputFileName: str) -> callable:
        """
        Read the spectrum from a CSV file.

        Parameters:
        - inputFileName: str
            The name of the input file.

        Returns:
        - callable
            The interpolated spectrum.
        """
        inputDF = pd.read_csv(inputFileName, header=None)
        rigidityList = inputDF[0]  # GV
        fluxList = inputDF[1]  # p/m2/sr/s/GV
        fluxListcm2 = fluxList / (100 ** 2)

        rigiditySpec = scipy.interpolate.interp1d(rigidityList, fluxListcm2, kind="linear",
                                                  fill_value=0.0, bounds_error=False)
        
        return rigiditySpec

class DLRmodelSpectrum(rigiditySpectrum):
    """
    DLR model rigidity spectrum.
    """

    def __init__(self, atomicNumber: int, date_and_time: 'datetime' = None, OULUcountRateInSeconds: float = None, W_parameter: float = None):
        """
        Initialize the DLR model spectrum.

        Parameters:
        - atomicNumber: int
            The atomic number of the particle.
        - date_and_time: datetime, optional
            The date and time for the spectrum.
        - OULUcountRateInSeconds: float, optional
            The OULU count rate in seconds.
        - W_parameter: float, optional
            The W parameter for the DLR model.
        """
        if not sum([(date_and_time is not None), (OULUcountRateInSeconds is not None), (W_parameter is not None)]) == 1:
            print("Error: exactly one supplied input out of the date and time, OULU count rate per second or the W parameter, to the DLR model spectrum must be given!")
            raise Exception

        if date_and_time is not None:
            self._generatedSpectrumDF = CosRayModifiedISO.getSpectrumUsingTimestamp(timestamp=date_and_time, atomicNumber=atomicNumber)

        if W_parameter is not None:
            self._generatedSpectrumDF = CosRayModifiedISO.getSpectrumUsingSolarModulation(solarModulationWparameter=W_parameter, atomicNumber=atomicNumber)

        if OULUcountRateInSeconds is not None:
            self._generatedSpectrumDF =CosRayModifiedISO.getSpectrumUsingOULUcountRate(OULUcountRatePerSecond=OULUcountRateInSeconds, atomicNumber=atomicNumber)

        self.rigiditySpec = interp1d(x=self._generatedSpectrumDF["Rigidity (GV/n)"],
                                     y=self._generatedSpectrumDF["d_Flux / d_R (cm-2 s-1 sr-1 (GV/n)-1)"],
                                     kind="linear",
                                     bounds_error=False,
                                     fill_value=(0.0, 0.0))

        with open("DLRsavedSpectrum.pkl", "wb") as DLRspecFile:
            pkl.dump(self, DLRspecFile)

class CommonModifiedPowerLawSpectrum(rigiditySpectrum):
    """
    Common modified power law rigidity spectrum.
    """

    specIndexModification = lambda self, P: self.deltaGamma * (P - 1)
    rigiditySpec = lambda self, P: self.J0 * self.step_function(P, self.lowerLimit, self.upperLimit) * (P ** (-(self.gamma + self.specIndexModification(P)))) / (100 ** 2)  # cm-2 s-1 sr-1 GV-1 : converted from m-2 to cm-2

    def __init__(self, J0: float, gamma: float, deltaGamma: float, lowerLimit: float = -np.inf, upperLimit: float = np.inf):
        """
        Initialize the common modified power law spectrum.

        Parameters:
        - J0: float
            The normalization factor.
        - gamma: float
            The spectral index.
        - deltaGamma: float
            The modification factor for the spectral index.
        - lowerLimit: float, optional
            The lower limit for the rigidity.
        - upperLimit: float, optional
            The upper limit for the rigidity.
        """
        self.lowerLimit = lowerLimit
        self.upperLimit = upperLimit

        self.J0 = J0  # m-2 s-1 sr-1 GV-1
        self.gamma = gamma
        self.deltaGamma = deltaGamma

        with open("CommonsavedSpectrum.pkl", "wb") as CommonspecFile:
            pkl.dump(self, CommonspecFile)

    def step_function(self, rigidity: float, lowerLimit: float, upperLimit: float) -> float:
        """
        Step function for the rigidity spectrum.

        Parameters:
        - rigidity: float
            The rigidity.
        - lowerLimit: float
            The lower limit for the rigidity.
        - upperLimit: float
            The upper limit for the rigidity.

        Returns:
        - float
            The value of the step function.
        """
        if (rigidity >= lowerLimit) and (rigidity <= upperLimit):
            return 1.0
        else:
            return 0.0

class CommonModifiedPowerLawSpectrumSplit(rigiditySpectrum):
    """
    Common modified power law rigidity spectrum with split spectral index modification.
    """

    specIndexModification_high = lambda self, P: self.deltaGamma * (P - 1)
    specIndexModification_low = lambda self, P: self.deltaGamma * (P)
    specIndexModification = lambda self, P: self.specIndexModification_high(P) if P > 1.0 else self.specIndexModification_low(P)
    rigiditySpec = lambda self, P: self.J0 * (P ** (-(self.gamma + self.specIndexModification(P)))) / (100 ** 2)  # cm-2 s-1 sr-1 GV-1 : converted from m-2 to cm-2

    def __init__(self, J0: float, gamma: float, deltaGamma: float):
        """
        Initialize the common modified power law spectrum with split spectral index modification.

        Parameters:
        - J0: float
            The normalization factor.
        - gamma: float
            The spectral index.
        - deltaGamma: float
            The modification factor for the spectral index.
        """
        self.J0 = J0  # m-2 s-1 sr-1 GV-1
        self.gamma = gamma
        self.deltaGamma = deltaGamma

        with open("CommonsavedSpectrum.pkl", "wb") as CommonspecFile:
            pkl.dump(self, CommonspecFile)
