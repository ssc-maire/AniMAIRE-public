

class particleSpecies():

    particle_atomic_number_dict = {"proton":1,"alpha":2}

    def __init__(self, particleName, atomicNumber=None):
        self.particleName = particleName
        if atomicNumber == None:
            atomicNumber = self.particle_atomic_number_dict[particleName]
        
        self.atomicNumber = atomicNumber