from .particleSpecies import particleSpecies
from .momentaDistribution import momentaDistribution
from .rigiditySpectrum import rigiditySpectrum
from .pitchAngleDistribution import pitchAngleDistribution

class particleDistribution():
    """
    Class representing a particle distribution.
    """

    def __init__(self, 
                 particle_name: str, 
                 rigidity_spectrum: rigiditySpectrum, 
                 pitch_angle_distribution: pitchAngleDistribution):
        """
        Initialize the particle distribution.

        Parameters:
        - particle_name: str
            The name of the particle.
        - rigidity_spectrum: rigiditySpectrum
            The rigidity spectrum of the particle.
        - pitch_angle_distribution: pitchAngleDistribution
            The pitch angle distribution of the particle.
        """
        self.particle_species = particleSpecies(particle_name)
        self.momentum_distribution = momentaDistribution(rigidity_spectrum,
                                                         pitch_angle_distribution)