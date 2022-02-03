#!/usr/bin/env python

import openmm
from openmm import unit
import numpy as np
from scipy import optimize

def parse_arguments():

    """
    Parse command line arguments.
    """

    import argparse

    parser = argparse.ArgumentParser(
        description="Python script for minimizing unit cell."
        )

    parser.add_argument(
        '--input', 
        "-i", 
        type=str, 
        help="Input xml file", 
        required=True
        )

    parser.add_argument(
        '--pdb', 
        "-p", 
        type=str, 
        help="Input pdb file", 
        required=True
        )

    parser.add_argument(
        '--output', 
        "-o", 
        type=str, 
        help="Output xml file.", 
        required=True
        )

    return parser.parse_args()


class ContextWrapper(object):

    """
    Wrapping around openmm context and expose handy functionalities
    to facilitate minimization.
    """
    
    def __init__(self, context, number_of_atoms):
        
        self.context = context
        self.number_of_atoms = number_of_atoms
        self.N_crds = self.number_of_atoms * 3
    
    @property
    def pos_box_flat(self):

        """
        Flat vector with concatenated particle and box vectors.
        """
        
        state = self.context.getState(
            getPositions=True,
        )
        pos = state.getPositions(asNumpy=True)
        box = state.getPeriodicBoxVectors(asNumpy=True)
            
        pos_box_flat = self.pos_box_to_flat(pos, box)
        
        return x
    
    @property
    def box_flat(self):

        """
        Flat vector with box vectors.
        """
        
        state = self.context.getState()
        box   = state.getPeriodicBoxVectors(asNumpy=True)

        boxflat = self.box_to_flat(box)
        
        return boxflat
        
    def box_to_flat(self, box):

        """
        Transform 3x3 box vector to flat vector.
        """
        
        boxflat = np.zeros(6, dtype=float) * unit.nanometer
        boxflat[0]   = box[0,0]
        boxflat[1:3] = box[1,:][:2]
        boxflat[3:]  = box[2,:]
        
        return boxflat            
    
    def pos_box_to_flat(self, pos, box):

        """
        Concatenate (N_atoms,3) position vector and box vector to flat
        vector.
        """
        
        boxflat = self.box_to_flat(box)
        pos_box_flat = np.hstack((pos.flatten(), boxflat._value))

        return pos_box_flat                
    
    def flat_to_pos_box(self, flat_pos_box):

        """
        Unravel flat concatenated vector with position and box to
        (N_atoms,3) position vector and 3x3 box vectors
        """
        
        pos_flat = flat_pos_box[:self.N_crds] * unit.nanometer
        pos = pos_flat.reshape((self.number_of_atoms, 3))

        box_flat = flat_pos_box[self.N_crds:] * unit.nanometer
        box = np.zeros((3,3), dtype=float) * unit.nanometer
        box[0,0] = box_flat[0]
        box[1,0] = box_flat[1]
        box[1,1] = box_flat[2]
        box[2,:] = box_flat[3:]
            
        return pos, box

    def flatbox_to_box(self, boxflat):

        """
        Unravel flat box vector to 3x3 box vector
        """
        
        boxflat_cp = np.copy(boxflat[:]) * unit.nanometer
        box = np.zeros((3,3), dtype=float) * unit.nanometer
        box[0,0] = boxflat_cp[0]
        box[1,0] = boxflat_cp[1]
        box[1,1] = boxflat_cp[2]
        box[2,:] = boxflat_cp[3:]
        
        return box 
    
    def ene_box(self, boxflat):

        """
        Compute system energy from flat boxvector. Energy in kJ/mol.
        Note this also alters the positions.
        """
        
        box_new = self.flatbox_to_box(boxflat)
        state = self.context.getState(
            getEnergy=True,
            getPositions=True,
        )

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)

        box_old_inv = np.linalg.inv(box_old)        
        pos_frac    = np.matmul(box_old_inv, pos_old.T).T
        pos_new     = np.matmul(box_new, pos_frac.T).T * unit.nanometer

        self.context.setPositions(pos_new)
        self.context.setPeriodicBoxVectors(*box_new)

        ene = state.getPotentialEnergy()
        return ene._value
    
    def ene(self, pos_box_flat):

        """
        Compute system energy from flat position and boxvector. Energy in kJ/mol
        """
        
        pos, box = self.flat_to_pos_box(pos_box_flat)
        
        self.context.setPositions(pos)
        self.context.setPeriodicBoxVectors(*box)
            
        state = self.context.getState(getEnergy=True)
        ene   = state.getPotentialEnergy()
        
        return ene._value
    
    def grad_box(self, boxflat):

        """
        Compute energy gradient from flat box vector. Gradient in kJ/mol/nm
        """
        
        epsilon = 1e-5 * unit.nanometer
        
        box = self.flatbox_to_box(boxflat)
        state = self.context.getState(
            getPositions=True,
        )
        pos = state.getPositions(asNumpy=True)

        box_inv  = np.linalg.inv(box)        
        pos_frac = np.matmul(box_inv, pos.T).T

        self.context.setPositions(pos)
        self.context.setPeriodicBoxVectors(*box)
        
        boxgrad_flat = np.zeros(6, dtype=float) * unit.kilojoule_per_mole / unit.nanometer
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box[r,c] += epsilon
                pos_new   = np.matmul(box, pos_frac.T).T * unit.nanometer
                x_new     = self.pos_box_to_flat(pos_new, box)
                ene_plus  = self.ene(x_new) * unit.kilojoule_per_mole

                box[r,c] -= 2. * epsilon
                pos_new   = np.matmul(box, pos_frac.T).T * unit.nanometer
                x_new     = self.pos_box_to_flat(pos_new, box)
                ene_minus = self.ene(x_new) * unit.kilojoule_per_mole

                boxgrad_flat[grad_idx]  = (ene_plus - ene_minus)/(epsilon * 2.)

                box[r,c] += epsilon

                grad_idx += 1

        self.context.setPeriodicBoxVectors(*box)
        self.context.setPositions(pos)
        boxgrad_flat = boxgrad_flat._value

        return boxgrad_flat

    def grad(self, pos_box_flat):

        """
        Compute energy gradient from flat position and boxvector. Gradient in kJ/mol/nm
        """        
        
        epsilon = 1e-5 * unit.nanometer
        
        pos, box = self.get_pos_box(pos_box_flat)
        
        box_inv  = np.linalg.inv(box)        
        pos_frac = np.matmul(box_inv, pos.T).T

        self.context.setPositions(pos)
        self.context.setPeriodicBoxVectors(*box)
        state  = self.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True)._value
        
        boxgrad_flat = np.zeros(6, dtype=float) * unit.kilojoule_per_mole / unit.nanometer
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box[r,c] += epsilon
                pos_new   = np.matmul(box, pos_frac.T).T * unit.nanometer
                x_new     = self.pos_box_to_flat(pos_new, box)
                ene_plus  = self.ene(x_new) * unit.kilojoule_per_mole
                
                box[r,c] -= 2. * epsilon
                pos_new   = np.matmul(box, pos_frac.T).T * unit.nanometer
                x_new     = self.pos_box_to_flat(pos_new, box)
                ene_minus = self.ene(x_new) * unit.kilojoule_per_mole
                
                boxgrad_flat[grad_idx]  = (ene_plus - ene_minus)/(epsilon * 2.)

                box[r,c] += epsilon

                grad_idx += 1

        self.context.setPeriodicBoxVectors(*box)
        self.context.setPositions(pos)
        boxgrad_flat = boxgrad_flat._value
        
        forces *= -1
        grad    = np.hstack((forces.flatten(), boxgrad_flat))        

        return grad


def run_xtal_min(
    xml_path, 
    pdb_path,
    platform_name = "CUDA",
    property_dict = {
        "Precision" : "mixed"
    }
    ):

    """
    Run xtal minimizatin. Return final context.
    """
    
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()
    system = openmm.XmlSerializer.deserialize(xml_str)

    integrator = openmm.LangevinIntegrator(
        300. * unit.kelvin,
        1./unit.picoseconds,
        0.002 * unit.picoseconds
    )
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    context = openmm.Context(
        system,
        integrator,
        platform,
    )
    
    context.setPositions(pdbfile.positions)
    context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())

    cw = ContextWrapper(
        context, 
        system.getNumParticles()
        )
    x0 = np.copy(cw.box_flat)
    ene2_t = cw.ene_box(x0)
    for _ in range(500):
        openmm.LocalEnergyMinimizer.minimize(context)
        x0 = np.copy(cw.box_flat)
        ene1 = cw.ene_box(x0)
        result = optimize.minimize(
            fun=cw.ene_box,
            x0=x0, 
            method="BFGS",
            jac=cw.grad_box
        )
        context.setPeriodicBoxVectors(
            *cw.flatbox_to_box(
                result.x
            )
        )
        ene2 = cw.ene_box(result.x)
        delta1 = ene2 - ene1
        delta2 = ene2 - ene2_t
        ene2_t = ene2
        print("Delta1", delta1, "Delta2", delta2)
        

    return context


def main()

    """
    Run the main workflow.
    """

    args = parse_arguments()

    ### Note: CUDA is much faster even for
    ###       the wrapped minimization here (about factor 10)
    context = run_xtal_min(
        args.input,
        args.pdb
        platform_name = "CUDA",
        property_dict = {
            #"Threads"             : '4',
            "DeterministicForces" : "True"
        }
    )

    state = context.getState(getPositions=True)
    with open(args.output, "w") as fopen:
        fopen.write(
            openmm.XmlSerializer.serialize(state)
            )

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()