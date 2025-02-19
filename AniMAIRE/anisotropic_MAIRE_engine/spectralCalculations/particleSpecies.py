class particleSpecies():
    """
    Class representing a particle species.
    """

    particle_atomic_number_dict = {"proton": 1, "alpha": 2}
    particle_atomic_mass_dict = {"proton": 1, "alpha": 4}

    def __init__(self, particleName: str, atomicNumber: int = None):
        """
        Initialize the particle species.

        Parameters:
        - particleName: str
            The name of the particle.
        - atomicNumber: int, optional
            The atomic number of the particle. If not provided, it will be looked up from the dictionary.
        """
        self.particleName = particleName
        if atomicNumber is None:
            atomicNumber = self.particle_atomic_number_dict[particleName]
        
        self.atomicNumber = atomicNumber
        self.atomicMass = self.particle_atomic_mass_dict[particleName]