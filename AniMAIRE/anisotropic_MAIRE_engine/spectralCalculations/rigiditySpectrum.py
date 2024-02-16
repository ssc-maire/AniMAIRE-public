
from scipy.interpolate import interp1d
import pandas as pd
import scipy
import copy

import pickle as pkl
import numpy as np

from CosRayModifiedISO import CosRayModifiedISO

print("WARNING: currently unknown whether the reported spectral weighting factor should be in terms of energy or rigidityd")

class rigiditySpectrum():

    def __init__(self):
        pass

    def __call__(self, x):
        return self.rigiditySpec(x)
    
    def __add__(self, right):
        summed_spectrum = rigiditySpectrum()
        summed_spectrum.rigiditySpec = lambda x:self.rigiditySpectrum(x) + right.rigiditySpectrum(x)
        return summed_spectrum

class powerLawSpectrum(rigiditySpectrum):

    rigiditySpec = lambda self,x:self.normalisationFactor * (x**self.spectralIndex)

    def __init__(self, normalisationFactor, spectralIndex):

        self.normalisationFactor = normalisationFactor
        self.spectralIndex = spectralIndex

class interpolatedInputFileSpectrum(rigiditySpectrum):

    def __init__(self, inputFileName):

        self.inputFilename = inputFileName
        self.rigiditySpec = self.readSpecFromCSV(self.inputFilename)

    def readSpecFromCSV(self, inputFileName):

        inputDF = pd.read_csv(inputFileName,header=None)
        rigidityList = inputDF[0] # GV
        fluxList = inputDF[1] # p/m2/sr/s/GV
        fluxListcm2 = fluxList / (100**2)

        #rigiditySpec = scipy.interpolate.interp1d(rigidityList,fluxListcm2,kind="linear",fill_value="extrapolate")
        rigiditySpec = scipy.interpolate.interp1d(rigidityList,fluxListcm2,kind="linear",
                                                  fill_value=0.0, bounds_error=False)
        
        return rigiditySpec

class DLRmodelSpectrum(rigiditySpectrum):

    def __init__(self, atomicNumber, date_and_time=None, OULUcountRateInSeconds=None, W_parameter=None):

        if not sum([(date_and_time is not None),(OULUcountRateInSeconds is not None),(W_parameter is not None)]) == 1:
            print("Error: exactly one supplied input out of the date and time, OULU count rate per second or the W paramer, to the DLR model spectrum must be given!")
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
                                                        fill_value = (0.0,0.0))

        with open("DLRsavedSpectrum.pkl","wb") as DLRspecFile:
            pkl.dump(self,DLRspecFile)

class MishevModifiedPowerLawSpectrum(rigiditySpectrum):

    specIndexModification = lambda self,P:self.deltaGamma * (P-1)
    rigiditySpec = lambda self,P:self.J0 * self.step_function(P, self.lowerLimit, self.upperLimit) * (P**(-(self.gamma + self.specIndexModification(P)))) / (100**2) #cm-2 s-1 sr-1 GV-1 : converted from m-2 to cm-2

    def __init__(self,J0,gamma,deltaGamma,lowerLimit=-1 * np.inf,upperLimit=np.inf):

        self.lowerLimit = lowerLimit
        self.upperLimit = upperLimit

        self.J0 = J0 #m-2 s-1 sr-1 GV-1
        self.gamma = gamma
        self.deltaGamma = deltaGamma

        with open("MishevsavedSpectrum.pkl","wb") as MishevspecFile:
            pkl.dump(self,MishevspecFile)

    def step_function(self, rigidity, lowerLimit, upperLimit):
        
        if (rigidity >= lowerLimit) and (rigidity <= upperLimit):
            return 1.0
        else:
            return 0.0

class MishevModifiedPowerLawSpectrumSplit(rigiditySpectrum):

    specIndexModification_high = lambda self,P:self.deltaGamma * (P-1)
    specIndexModification_low = lambda self,P:self.deltaGamma * (P)
    specIndexModification = lambda self,P:self.specIndexModification_high(P) if P > 1.0 else self.specIndexModification_low(P)
    rigiditySpec = lambda self,P:self.J0 * (P**(-(self.gamma + self.specIndexModification(P)))) / (100**2) #cm-2 s-1 sr-1 GV-1 : converted from m-2 to cm-2

    def __init__(self,J0,gamma,deltaGamma):

        self.J0 = J0 #m-2 s-1 sr-1 GV-1
        self.gamma = gamma
        self.deltaGamma = deltaGamma

        with open("MishevsavedSpectrum.pkl","wb") as MishevspecFile:
            pkl.dump(self,MishevspecFile)
