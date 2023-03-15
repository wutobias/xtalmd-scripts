import openmm
from openmm import unit


class Logger(object):
    
    def __init__(self, system, filehandle):
        
        self.dof  = system.getNumParticles() * 3.
        self.dof -= system.getNumConstraints()
        if any(isinstance(force, openmm.CMMotionRemover) for force in system.getForces()):
            self.dof -= 3.
        self.total_mass = 0. * unit.dalton
        for i in range(system.getNumParticles()):
            self.total_mass += system.getParticleMass(i)
            
        if isinstance(filehandle, str):
            self.write_str = True
            self.str = ""
        else:
            self.write_str = False
            self.filehandle = filehandle
        self.write_header()

    def write_header(self):
        
        if self.write_str:
            self.str += "# Time,Potential,Kinetic En.,"
            self.str += "Box-XX,Box-XY,Box-XZ,"
            self.str += "Box-YX,Box-YY,Box-YZ,"
            self.str += "Box-ZX,Box-ZY,Box-ZZ,"
            self.str += "Volume,Density,Temperature"
            self.str += "\n"
        else:
            self.filehandle.write("# Time,Potential,Kinetic En.,")
            self.filehandle.write("Box-XX,Box-XY,Box-XZ,")
            self.filehandle.write("Box-YX,Box-YY,Box-YZ,")
            self.filehandle.write("Box-ZX,Box-ZY,Box-ZZ,")
            self.filehandle.write("Volume,Density,Temperature")
            self.filehandle.write("\n")
        
    def write(self, state):
        
        time = state.getTime()
        time = time.in_units_of(unit.picosecond)
        vol = state.getPeriodicBoxVolume()
        vol = vol.in_units_of(unit.nanometer**3)
        boxvectors = state.getPeriodicBoxVectors()
        kin_ene = state.getKineticEnergy()
        pot_ene = state.getPotentialEnergy()
        temp = 2. * kin_ene / (self.dof * unit.MOLAR_GAS_CONSTANT_R)
        dens = self.total_mass / vol
        dens = dens.in_units_of(unit.grams/unit.item/unit.milliliter)

        if self.write_str:
            self.str += f"{time._value},{pot_ene._value},{kin_ene._value},"
            self.str += f"{boxvectors[0][0]._value},{boxvectors[0][1]._value},{boxvectors[0][2]._value},"
            self.str += f"{boxvectors[1][0]._value},{boxvectors[1][1]._value},{boxvectors[1][2]._value},"
            self.str += f"{boxvectors[2][0]._value},{boxvectors[2][1]._value},{boxvectors[2][2]._value},"
            self.str += f"{vol._value},{dens._value},{temp._value}"
            self.str += "\n"
        else:
            self.filehandle.write(f"{time._value},{pot_ene._value},{kin_ene._value},")
            self.filehandle.write(f"{boxvectors[0][0]._value},{boxvectors[0][1]._value},{boxvectors[0][2]._value},")
            self.filehandle.write(f"{boxvectors[1][0]._value},{boxvectors[1][1]._value},{boxvectors[1][2]._value},")
            self.filehandle.write(f"{boxvectors[2][0]._value},{boxvectors[2][1]._value},{boxvectors[2][2]._value},")
            self.filehandle.write(f"{vol._value},{dens._value},{temp._value}")
            self.filehandle.write("\n")

    def write_footer(self, text):

        if self.write_str:
            self.str += str(text)
            self.str += '\n'
        else:
            self.filehandle.write(str(text))
            self.filehandle.write("\n")

