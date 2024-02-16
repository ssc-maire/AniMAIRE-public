from .particleSpecies import particleSpecies
from .momentaDistribution import momentaDistribution
from .rigiditySpectrum import rigiditySpectrum
from .pitchAngleDistribution import pitchAngleDistribution

class particleDistribution():

    def __init__(self, 
                 particle_name:str, 
                 rigidity_spectrum:rigiditySpectrum, 
                 pitch_angle_distribution:pitchAngleDistribution):

        self.particle_species = particleSpecies(particle_name)
        self.momentum_distribution = momentaDistribution(rigidity_spectrum,
                                                         pitch_angle_distribution)