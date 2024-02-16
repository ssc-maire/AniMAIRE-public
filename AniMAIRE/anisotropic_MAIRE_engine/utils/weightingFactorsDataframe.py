import pandas as pd
import tqdm

import ParticleRigidityCalculationTools as PRCT

from pythonScripts.spectralCalculations.momentaDistribution import momentaDistribution

class WeightingFactorsSeries(pd.Series):
    @property
    def _constructor(self):
        return WeightingFactorsSeries

    @property
    def _constructor_expanddim(self):
        return WeightingFactorsDataFrame

class WeightingFactorsDataFrame(pd.DataFrame):

    @property
    def _constructor(self):
        return WeightingFactorsDataFrame

    @property
    def _constructor_sliced(self):
        return WeightingFactorsSeries

    def __init__(self,  AsymptoticDirectionDataframe, momentaDist:momentaDistribution): #formerly acquireWeightingFactors

        # initialise as empty dataframe
        super().__init__()

        # fill with values from original dataframe
        self.fillValuesFromDataFrame(AsymptoticDirectionDataframe)

        # find angle between asymptotic directions and IMF using dot product
        angleColumnTitle = "angleBetweenIMFinRadians"
        self.assignPitchAngles(momentaDist, angleColumnTitle)

        # find weighting factors from the angles and rigidities
        pitchAngleFunctionToUse = lambda row : momentaDist.getPitchAngleDistribution()(row[angleColumnTitle],row["Rigidity"])
        fullRigidityPitchWeightingFactorFunctionToUse = lambda row : momentaDist(row[angleColumnTitle],row["Rigidity"])

        self["PitchAngleWeightingFactor"] = self.apply(pitchAngleFunctionToUse,axis=1)
        self["RigidityWeightingFactor"] = self["Rigidity"].apply(momentaDist.getRigiditySpectrum()) * (self["Filter"] == 1)
        self["fullRigidityPitchWeightingFactor"] = self.apply(fullRigidityPitchWeightingFactorFunctionToUse,axis=1) * (self["Filter"] == 1)
        self["fullEnergyPitchWeightingFactor"] = PRCT.convertParticleRigiditySpecToEnergySpec(self["Rigidity"],self["fullRigidityPitchWeightingFactor"],
                                                                                       particleMassAU=1, particleChargeAU=1)

    def fillValuesFromDataFrame(self,dataframeToFillFrom:pd.DataFrame):

        self["Rigidity"] = dataframeToFillFrom["Rigidity"]
        self["Energy"] = PRCT.convertParticleRigidityToEnergy(dataframeToFillFrom["Rigidity"], particleMassAU = 1, particleChargeAU = 1)
        self["Filter"] = dataframeToFillFrom["Filter"]
        #print(spacepy.coordinates.Coords(dataframeToFillFrom.iloc[0]["Lat"],))

        print("assigning asymptotic coordinates")
        asymptoticDirectionList = []
        initialPositionList = []
        initialMomDirList = []
        for dataframeRow in tqdm.tqdm(dataframeToFillFrom.iterrows(),total=len(dataframeToFillFrom)):
          rowInSpacepyCoords = self.convertInitialCoordSetToSpacePy(dataframeRow[1], 100.0,"Lat","Long")
          asymptoticDirectionList.append(rowInSpacepyCoords)

        print("successfully converted asymptotic directions")

        self["initialLatitude"] = dataframeToFillFrom["initialLatitude"]
        self["initialLongitude"] = dataframeToFillFrom["initialLongitude"]
        self["initialZenith"] = dataframeToFillFrom["initialZenith"]
        self["initialAzimuth"] = dataframeToFillFrom["initialAzimuth"]

        self["Asymptotic Direction"] = asymptoticDirectionList