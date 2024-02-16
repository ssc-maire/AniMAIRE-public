import numpy as np
from spacepy.coordinates import Coords as spaceCoords
import pickle as pkl
import copy

class pitchAngleDistribution():

    def __init__(self, 
                 pitch_angle_distribution:callable, 
                 reference_latitude_in_GSM, 
                 reference_longitude_in_GSM):

        self.pitchAngleDistFunction = pitch_angle_distribution

        self.interplanetary_mag_field = spaceCoords([100.0,
                                                     reference_latitude_in_GSM, 
                                                     reference_longitude_in_GSM],
                                                     "GSM","sph")

    def __call__(self,pitchAngle, rigidity):
        #print(f"pitch weighting factor: {self.pitchAngleDistFunction(pitchAngle, rigidity)}")
        return self.pitchAngleDistFunction(pitchAngle, rigidity)
    
    def __add__(self, right):
        summed_dist = copy.deepcopy(self)
        summed_dist.pitchAngleDistFunction = lambda x,y:self.pitchAngleDistFunction(x,y) + right.pitchAngleDistFunction(x,y)
        return summed_dist
    
    def __mul__(self, right:float):
        multiplied_dist = copy.deepcopy(self)
        multiplied_dist.pitchAngleDistFunction = lambda x,y:right * self.pitchAngleDistFunction(x,y)
        return multiplied_dist

    __rmul__ = __mul__

class cosinePitchAngleDistribution(pitchAngleDistribution):

    pitchAngleDistFunction = lambda self,pitchAngle, rigidity:np.abs(0.5*np.sin(2*pitchAngle))

    def __init__(self):
        pass

class isotropicPitchAngleDistribution(pitchAngleDistribution):

    #pitchAngleDistFunction = lambda self,pitchAngle,theta:1/np.pi
    pitchAngleDistFunction = lambda self,pitchAngle,rigidity:1
    #pitchAngleDistFunction = lambda self,pitchAngle,rigidity:np.sin(2.0 * pitchAngle)

    def __init__(self):
        # with open("DLRsavedPAD.pkl","wb") as DLRPADFile:
        #     pkl.dump(self,DLRPADFile)
        pass

class gaussianPitchAngleDistribution(pitchAngleDistribution):

    pitchAngleDistFunction = lambda self,pitchAngle,rigidity:self.normFactor * np.exp(-(pitchAngle - self.alpha)**2/(self.sigma**2))

    def __init__(self,normFactor,sigma, alpha=0.0):
        self.normFactor = normFactor
        self.sigma = sigma
        self.alpha = alpha
        with open("MishevsavedPAD.pkl","wb") as MishevPADFile:
            pkl.dump(self,MishevPADFile)

class gaussianBeeckPitchAngleDistribution(pitchAngleDistribution):

    pitchAngleDistFunction = lambda self,pitchAngle_radians,rigidity:self.normFactor * np.exp((-0.5 * (pitchAngle_radians - (np.sin(pitchAngle_radians) * np.cos(pitchAngle_radians)))) / \
                                                                (self.A - (0.5 * (self.A - self.B) * (1-np.cos(pitchAngle_radians)))))

    def __init__(self,normFactor,A, B):
        self.normFactor = normFactor
        self.A = A
        self.B = B