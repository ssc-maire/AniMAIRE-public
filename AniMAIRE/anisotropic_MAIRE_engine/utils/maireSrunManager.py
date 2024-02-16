import os
import numpy as np

onekftinkm = 0.3048

class maireSrunManager():

    MAIREHOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    WORKDIR =MAIREHOME + "/developmentWorkDir/"
    #datetimeToRunAcross = dt.datetime(year=2010, month=5,day=24, hour=3, minute=20,second=34)
    _listOfDoseRateTypes = ["adose","edose","dosee","tn1","tn2","tn3","SEU","SEL"]

    generalIterativeInteger = 0

    def __init__(self):
        print("MAIREHOME is",self.MAIREHOME)

    def getMAIREHOME(self):
        return self.MAIREHOME

    def getWORKDIR(self):
        return self.WORKDIR

    def setRunMode(self, runMode):
        self._runMode = runMode

    def getKpIndex(self):
        return self._runMode

    def setDATETIMETORUNACROSS(self, datetimeToRunAcross):
        self._datetimeToRunAcross = datetimeToRunAcross

    def getDATETIMETORUNACROSS(self):
        return self._datetimeToRunAcross

    def setKpIndex(self, KpIndex):
        self._KpIndex = KpIndex

    def getKpIndex(self):
        return self._KpIndex

    def getGeneralIterativeInteger(self):

        intToReturn = self.generalIterativeInteger
        self.generalIterativeInteger = self.generalIterativeInteger + 1

        return intToReturn

    def setRigidityValues(self, rigidityValues:np.array):
        self._rigidityValues = rigidityValues

_maireSrunManager = maireSrunManager()