

from .anisotropic_MAIRE_engine.spectralCalculations.particleDistribution import particleDistribution
from .anisotropic_MAIRE_engine.spectralCalculations.pitchAngleDistribution import isotropicPitchAngleDistribution, pitchAngleDistribution


def get_correctly_formatted_particle_dist_list(proton_rigidity_spectrum, 
                                               alpha_rigidity_spectrum, 
                                               reference_pitch_angle_latitude, 
                                               reference_pitch_angle_longitude, 
                                               proton_pitch_angle_distribution, 
                                               alpha_pitch_angle_distribution):
    list_of_particle_distributions = []
    for particle_name, particle_rigidity_spectrum, particle_pitch_angle_distribution in [["proton",proton_rigidity_spectrum,proton_pitch_angle_distribution], 
                                                                                         ["alpha",alpha_rigidity_spectrum,alpha_pitch_angle_distribution]]:
        if not particle_rigidity_spectrum is None:
            if not type(particle_pitch_angle_distribution).__name__ == 'isotropicPitchAngleDistribution':
                if (reference_pitch_angle_latitude is None) or (reference_pitch_angle_longitude is None):
                    raise Exception("Error: a non-isotropic pitch angle distribution was used but reference longitudes and latitudes were not set!")

                particle_pitch_angle_distribution = pitchAngleDistribution(particle_pitch_angle_distribution,
                                                                  reference_pitch_angle_latitude,
                                                                  reference_pitch_angle_longitude)
                
            particle_distribution = particleDistribution(
                                                particle_name,
                                                particle_rigidity_spectrum,
                                                particle_pitch_angle_distribution)
            
            list_of_particle_distributions.append(particle_distribution)
    return list_of_particle_distributions

