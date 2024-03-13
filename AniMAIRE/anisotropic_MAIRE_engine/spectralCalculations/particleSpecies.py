

class particleSpecies():

    particle_atomic_number_dict = {"proton":1,"alpha":2}
    particle_atomic_mass_dict = {"proton":1,"alpha":4}

    def __init__(self, particleName, atomicNumber=None):
        self.particleName = particleName
        if atomicNumber == None:
            atomicNumber = self.particle_atomic_number_dict[particleName]
        
        self.atomicNumber = atomicNumber
        self.atomicMass = self.particle_atomic_mass_dict[particleName]