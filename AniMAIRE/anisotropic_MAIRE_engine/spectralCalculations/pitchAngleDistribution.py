import numpy as np
from spacepy.coordinates import Coords as spaceCoords
import pickle as pkl
import copy

class pitchAngleDistribution():
    """
    Base class for pitch angle distributions.
    """

    def __init__(self, 
                 pitch_angle_distribution: callable, 
                 reference_latitude_in_GSM: float, 
                 reference_longitude_in_GSM: float):
        """
        Initialize the pitch angle distribution.

        Parameters:
        - pitch_angle_distribution: callable
            Function describing the pitch angle distribution.
        - reference_latitude_in_GSM: float
            Reference latitude in GSM coordinates.
        - reference_longitude_in_GSM: float
            Reference longitude in GSM coordinates.
        """
        self.pitchAngleDistFunction = pitch_angle_distribution

        self.interplanetary_mag_field = spaceCoords([100.0,
                                                     reference_latitude_in_GSM, 
                                                     reference_longitude_in_GSM],
                                                     "GSM","sph")

    def __call__(self, pitchAngle: float, rigidity: float) -> float:
        """
        Evaluate the pitch angle distribution at a given pitch angle and rigidity.

        Parameters:
        - pitchAngle: float
            The pitch angle.
        - rigidity: float
            The rigidity.

        Returns:
        - float
            The value of the pitch angle distribution.
        """
        return self.pitchAngleDistFunction(pitchAngle, rigidity)
    
    def __add__(self, right: 'pitchAngleDistribution') -> 'pitchAngleDistribution':
        """
        Add two pitch angle distributions.

        Parameters:
        - right: pitchAngleDistribution
            The other pitch angle distribution.

        Returns:
        - pitchAngleDistribution
            The sum of the two pitch angle distributions.
        """
        summed_dist = copy.deepcopy(self)
        summed_dist.pitchAngleDistFunction = lambda x, y: self.pitchAngleDistFunction(x, y) + right.pitchAngleDistFunction(x, y)
        return summed_dist
    
    def __mul__(self, right: float) -> 'pitchAngleDistribution':
        """
        Multiply the pitch angle distribution by a scalar.

        Parameters:
        - right: float
            The scalar.

        Returns:
        - pitchAngleDistribution
            The scaled pitch angle distribution.
        """
        multiplied_dist = copy.deepcopy(self)
        multiplied_dist.pitchAngleDistFunction = lambda x, y: right * self.pitchAngleDistFunction(x, y)
        return multiplied_dist

    __rmul__ = __mul__

class cosinePitchAngleDistribution(pitchAngleDistribution):
    """
    Cosine pitch angle distribution.
    """

    pitchAngleDistFunction = lambda self, pitchAngle, rigidity: np.abs(0.5 * np.sin(2 * pitchAngle))

    def __init__(self):
        """
        Initialize the cosine pitch angle distribution.
        """
        pass

class isotropicPitchAngleDistribution(pitchAngleDistribution):
    """
    Isotropic pitch angle distribution.
    """

    pitchAngleDistFunction = lambda self, pitchAngle, rigidity: 1

    def __init__(self):
        """
        Initialize the isotropic pitch angle distribution.
        """
        pass

class gaussianPitchAngleDistribution(pitchAngleDistribution):
    """
    Gaussian pitch angle distribution.
    """

    pitchAngleDistFunction = lambda self, pitchAngle, rigidity: self.normFactor * np.exp(-(pitchAngle - self.alpha)**2 / (self.sigma**2))

    def __init__(self, normFactor: float, sigma: float, alpha: float = 0.0):
        """
        Initialize the Gaussian pitch angle distribution.

        Parameters:
        - normFactor: float
            The normalization factor.
        - sigma: float
            The standard deviation of the Gaussian distribution.
        - alpha: float, optional
            The mean of the Gaussian distribution.
        """
        self.normFactor = normFactor
        self.sigma = sigma
        self.alpha = alpha
        with open("CommonsavedPAD.pkl", "wb") as CommonPADFile:
            pkl.dump(self, CommonPADFile)

class gaussianBeeckPitchAngleDistribution(pitchAngleDistribution):
    """
    Gaussian Beeck pitch angle distribution.
    """

    pitchAngleDistFunction = lambda self, pitchAngle_radians, rigidity: self.normFactor * np.exp((-0.5 * (pitchAngle_radians - (np.sin(pitchAngle_radians) * np.cos(pitchAngle_radians)))) / \
                                                                (self.A - (0.5 * (self.A - self.B) * (1 - np.cos(pitchAngle_radians)))))

    def __init__(self, normFactor: float, A: float, B: float):
        """
        Initialize the Gaussian Beeck pitch angle distribution.

        Parameters:
        - normFactor: float
            The normalization factor.
        - A: float
            Parameter A for the distribution.
        - B: float
            Parameter B for the distribution.
        """
        self.normFactor = normFactor
        self.A = A
        self.B = B