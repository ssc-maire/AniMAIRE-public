import numpy as np

from .pitchAngleDistribution import pitchAngleDistribution
from .rigiditySpectrum import rigiditySpectrum

def jacobian_function_to_use(pitch_angle_in_radians: float) -> float:
    try:
        output_Jacobian_factor = 1 / np.sin(2.0 * pitch_angle_in_radians)
    except ZeroDivisionError:
        output_Jacobian_factor = 0.0
    return output_Jacobian_factor

class momentaDistribution():
    """
    Class representing a momenta distribution.
    """

    def __init__(self, rigidity_spectrum: rigiditySpectrum, pitch_angle_distribution: pitchAngleDistribution):
        """
        Initialize the momenta distribution.

        Parameters:
        - rigidity_spectrum: rigiditySpectrum
            The rigidity spectrum.
        - pitch_angle_distribution: pitchAngleDistribution
            The pitch angle distribution.
        """
        self.setRigiditySpectrum(rigidity_spectrum)
        self.setPitchAngleDistribution(pitch_angle_distribution)

    # setter methods

    def setPitchAngleDistribution(self, pitch_angle_distribution: pitchAngleDistribution):
        """
        Set the pitch angle distribution.

        Parameters:
        - pitch_angle_distribution: pitchAngleDistribution
            The pitch angle distribution.
        """
        self.pitch_angle_distribution = pitch_angle_distribution

    def setRigiditySpectrum(self, rigidity_spectrum: rigiditySpectrum):
        """
        Set the rigidity spectrum.

        Parameters:
        - rigidity_spectrum: rigiditySpectrum
            The rigidity spectrum.
        """
        self.rigidity_spectrum = rigidity_spectrum

    # getter methods

    def getPitchAngleDistribution(self) -> pitchAngleDistribution:
        """
        Get the pitch angle distribution.

        Returns:
        - pitchAngleDistribution
            The pitch angle distribution.
        """
        return self.pitch_angle_distribution

    def getRigiditySpectrum(self) -> rigiditySpectrum:
        """
        Get the rigidity spectrum.

        Returns:
        - rigiditySpectrum
            The rigidity spectrum.
        """
        return self.rigidity_spectrum

    # other Methods

    def __call__(self, pitchAngle: float, rigidity: float) -> float:
        """
        Evaluate the momenta distribution.

        Parameters:
        - pitchAngle: float
            The pitch angle.
        - rigidity: float
            The rigidity.

        Returns:
        - float
            The value of the momenta distribution.
        """
        pitch_angle_weighting_factor = self.pitch_angle_distribution(pitchAngle, rigidity) #* jacobian_function_to_use(pitchAngle)
        rigidity_weighting_factor = self.rigidity_spectrum(rigidity)

        full_rigidity_pitch_weighting_factor = pitch_angle_weighting_factor * rigidity_weighting_factor

        return full_rigidity_pitch_weighting_factor