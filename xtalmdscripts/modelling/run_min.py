#!/usr/bin/env python

import openmm
from openmm import unit
from openmm import app
from openmm.app.internal.unitcell import computePeriodicBoxVectors, computeLengthsAndAngles
import numpy as np
from scipy import optimize
from .utils import Logger

def parse_arguments():

    """
    Parse command line arguments.
    """

    import argparse

    parser = argparse.ArgumentParser(
        description="Python script for minimizing unit cell."
        )

    subparser  = parser.add_subparsers(dest='command')
    subparser.required = True
    yaml_parse = subparser.add_parser("yaml")
    xml_parse  = subparser.add_parser("xml")

    yaml_parse.add_argument(
        '--input', 
        "-i", 
        type=str, 
        help="Input yaml file", 
        required=True
        )

    xml_parse.add_argument(
        '--input', 
        "-i", 
        type=str, 
        help="Input xml file", 
        required=True
        )

    xml_parse.add_argument(
        '--pdb', 
        "-p", 
        type=str, 
        help="Input pdb file", 
        required=True
        )

    xml_parse.add_argument(
        '--prefix', 
        "-pre", 
        type=str, 
        help="Output prefix xml and pdb files.", 
        default="xtal_min",
        required=False
        )

    xml_parse.add_argument(
        '--steps', 
        "-s", 
        type=int, 
        help="Number of iterations.", 
        required=False,
        default=100
        )

    xml_parse.add_argument(
        '--alternating', 
        "-a", 
        action='store_true',
        help="Use alternating minimizations.", 
        required=False,
        default=False,
        )

    xml_parse.add_argument(
        '--use_lengths_and_angles', 
        "-la", 
        action='store_true',
        help="Use cell length and angles for minimizing.", 
        required=False,
        default=False,
        )

    xml_parse.add_argument(
        '--method', 
        "-m", 
        type=str, 
        help="Minimization method for box vectors.", 
        required=False,
        default="Nelder-Mead"
        )

    xml_parse.add_argument(
        '--epsilon', 
        "-e", 
        type=float, 
        help="Epsilon (in nm) used for finite difference gradient calcs.", 
        required=False,
        default=1.e-5
        )

    return parser.parse_args()


class ContextWrapper(object):

    """
    Wrapping around openmm context and expose handy functionalities
    to facilitate minimization.
    """
    
    def __init__(self, context, number_of_atoms, epsilon=1.e-5, use_lengths_and_angles=False):
        
        self.context = context
        self.number_of_atoms = number_of_atoms
        self.N_crds = self.number_of_atoms * 3
        self.use_lengths_and_angles = use_lengths_and_angles

        self._epsilon = epsilon * unit.nanometer
    
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
            
        pos_box_flat = self.pos_box_to_flat(
            pos.value_in_unit(unit.nanometer), 
            box.value_in_unit(unit.nanometer)
            )
        
        return pos_box_flat
    
    @property
    def box_flat(self):

        """
        Flat vector with box vectors.
        """
        
        state = self.context.getState()
        box   = state.getPeriodicBoxVectors(asNumpy=True)

        boxflat = self.box_to_flat(
            box.value_in_unit(unit.nanometer)
            )
        
        return boxflat
        
    def box_to_flat(self, box):

        """
        Transform 3x3 box vector to flat vector.
        """
        if self.use_lengths_and_angles:
            return computeLengthsAndAngles(box)

        boxflat      = np.zeros(6, dtype=float)
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
        pos_box_flat = np.hstack((pos.flatten(), boxflat))

        return pos_box_flat                
    
    def flat_to_pos_box(self, flat_pos_box):

        """
        Unravel flat concatenated vector with position and box to
        (N_atoms,3) position vector and 3x3 box vectors
        """
        
        pos_flat = flat_pos_box[:self.N_crds] * unit.nanometer
        pos = pos_flat.reshape((self.number_of_atoms, 3))
        box = self.flatbox_to_box(flat_pos_box[self.N_crds:])

        return pos, box

    def flatbox_to_box(self, boxflat):

        """
        Unravel flat box vector to 3x3 box vector
        """
        if self.use_lengths_and_angles:
            box = computePeriodicBoxVectors(*boxflat)
        else:
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

        box_new = app.internal.unitcell.reducePeriodicBoxVectors(box_new)
        box_new = np.array(box_new._value) * box_new.unit

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)

        box_old_inv = np.linalg.inv(box_old.value_in_unit(unit.nanometer))        
        pos_frac    = np.matmul(box_old_inv, pos_old.T).T
        pos_new     = np.matmul(box_new, pos_frac.T).T * unit.nanometer

        try:
            self.context.setPositions(pos_new)
            self.context.setPeriodicBoxVectors(*box_new)

            ene = state.getPotentialEnergy()
        except:
            boxold = state.getPeriodicBoxVectors(asNumpy=True)
            raise ValueError(
                f"Cannot set box {boxold} to {box_new}")

        return ene._value
    
    def ene(self, pos_box_flat):

        """
        Compute system energy from flat position and boxvector. Energy in kJ/mol
        """
        
        pos, box = self.flat_to_pos_box(pos_box_flat)
        state = self.context.getState(getEnergy=True)

        box = app.internal.unitcell.reducePeriodicBoxVectors(box)
        box = np.array(box._value) * box.unit

        try:
            self.context.setPositions(pos)
            self.context.setPeriodicBoxVectors(*box)

            ene = state.getPotentialEnergy()
        except:
            boxold = state.getPeriodicBoxVectors(asNumpy=True)
            raise ValueError(
                f"Cannot set box {boxold} to {box}")
        
        return ene._value
    
    def grad_box(self, boxflat):

        """
        Compute energy gradient from flat box vector. Gradient in kJ/mol/nm
        """
        
        epsilon = self._epsilon
        
        box = self.flatbox_to_box(boxflat)
        state = self.context.getState(
            getPositions=True,
        )
        pos = state.getPositions(asNumpy=True)

        box = app.internal.unitcell.reducePeriodicBoxVectors(box)
        box = np.array(box._value) * box.unit

        box_inv  = np.linalg.inv(box)        
        pos_frac = np.matmul(box_inv, pos.T).T

        try:
            self.context.setPositions(pos)
            self.context.setPeriodicBoxVectors(*box)
        except:
            boxold = state.getPeriodicBoxVectors(asNumpy=True)
            raise ValueError(
                f"Cannot set box {boxold} to {box}")
        
        boxgrad_flat = np.zeros(6, dtype=float) * unit.kilojoule_per_mole / unit.nanometer
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box[r,c] += epsilon
                box_plus  = app.internal.unitcell.reducePeriodicBoxVectors(box)
                box_plus  = np.array(box_plus._value) * box_plus.unit
                pos_new   = np.matmul(box_plus, pos_frac.T).T
                x_new     = self.pos_box_to_flat(
                    pos_new, 
                    box_plus.value_in_unit(unit.nanometer)
                    )
                ene_plus  = self.ene(x_new) * unit.kilojoule_per_mole

                box[r,c]  -= 2. * epsilon
                box_minus  = app.internal.unitcell.reducePeriodicBoxVectors(box)
                box_minus  = np.array(box_minus._value) * box_minus.unit
                pos_new    = np.matmul(box_minus, pos_frac.T).T
                x_new      = self.pos_box_to_flat(
                    pos_new, 
                    box_minus.value_in_unit(unit.nanometer)
                    )
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
        
        epsilon = self._epsilon
        
        pos, box = self.flat_to_pos_box(pos_box_flat)

        box = app.internal.unitcell.reducePeriodicBoxVectors(box)
        box = np.array(box._value) * box.unit
        
        box_inv  = np.linalg.inv(box)        
        pos_frac = np.matmul(box_inv, pos.T).T

        try:
            self.context.setPositions(pos)
            self.context.setPeriodicBoxVectors(*box)
        except:
            boxold = state.getPeriodicBoxVectors(asNumpy=True)
            raise ValueError(
                f"Cannot set box {boxold} to {box}")

        state  = self.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True)._value
        
        boxgrad_flat = np.zeros(6, dtype=float) * unit.kilojoule_per_mole / unit.nanometer
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box[r,c] += epsilon
                box_plus  = app.internal.unitcell.reducePeriodicBoxVectors(box)
                box_plus  = np.array(box_plus._value) * box_plus.unit
                pos_new   = np.matmul(box_plus, pos_frac.T).T
                x_new     = self.pos_box_to_flat(
                    pos_new, 
                    box_plus.value_in_unit(unit.nanometer)
                    )
                ene_plus  = self.ene(x_new) * unit.kilojoule_per_mole
                
                box[r,c]  -= 2. * epsilon
                box_minus  = app.internal.unitcell.reducePeriodicBoxVectors(box)
                box_minus  = np.array(box_minus._value) * box_minus.unit
                pos_new    = np.matmul(box_minus, pos_frac.T).T
                x_new      = self.pos_box_to_flat(
                    pos_new, 
                    box_minus.value_in_unit(unit.nanometer)
                    )
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
    steps = 100,
    epsilon = 1.e-5,
    alternating = False,
    use_lengths_and_angles = False,
    method = "Nelder-Mead",
    platform_name = "CUDA",
    property_dict = {
        "Precision" : "mixed"
    },
    prefix="xtal_min"
    ):

    """
    Run xtal minimizatin. Return final context that has minimal energy.
    """

    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np

    time_step = 1. * unit.femtoseconds
    constraint_tolerance = 1.e-6
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()
    system = openmm.XmlSerializer.deserialize(xml_str)

    integrator = openmm.LangevinIntegrator(
        300. * unit.kelvin,
        1./unit.picoseconds,
        time_step
    )
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    context = openmm.Context(
        system,
        integrator,
        platform,
    )
    
    box_reduced = topology.getPeriodicBoxVectors()
    context.setPositions(pdbfile.positions)
    context.setPeriodicBoxVectors(*box_reduced)

    cw = ContextWrapper(
        context, 
        system.getNumParticles(),
        epsilon,
        use_lengths_and_angles
        )
    logfile = Logger(system, "")

    ### Alternate between openmm native minimizer
    ### and xtal cell minimization.
    if alternating:
        x0 = np.copy(cw.box_flat)
        best_ene = 999999999999999999999.
        best_x   = None
        for _ in range(steps):
            openmm.LocalEnergyMinimizer.minimize(context)
            x0 = np.copy(cw.box_flat)
            ene1 = cw.ene_box(x0)
            try:
                result = optimize.minimize(
                    fun=cw.ene_box,
                    x0=x0, 
                    method=method,
                    jac=cw.grad_box,
                    options = {
                        "maxiter" : steps,
                        "disp"    : False,
                        "gtol"    : 1.e-5,
                    }

                )
            except Exception as e:
                return 0, e
            context.setPeriodicBoxVectors(
                *cw.flatbox_to_box(
                    result.x
                )
            )
            ene = cw.ene_box(result.x)
            if ene < best_ene:
                best_ene = ene
                best_x   = cw.pos_box_flat
            else:
                pos, box = cw.flat_to_pos_box(best_x)
                context.setPositions(pos)
                context.setPeriodicBoxVectors(*box)
            state = context.getState(
                getEnergy=True,
                getPositions=True
            )
            logfile.write(state)

        pos, box = cw.flat_to_pos_box(best_x)
        context.setPositions(pos)
        context.setPeriodicBoxVectors(*box)
        state = context.getState(getPositions=True)

    ### Minimize everything together.
    else:
        def callback(xk):
            pos, box = cw.flat_to_pos_box(xk)
            context.setPositions(pos)
            context.setPeriodicBoxVectors(*box)
            state = context.getState(
                getEnergy=True,
                getPositions=True
            )
            logfile.write(state)
            
        x0 = np.copy(cw.pos_box_flat)
        try:
            result = optimize.minimize(
                fun=cw.ene,
                x0=x0, 
                method=method,
                jac=cw.grad,
                callback=callback,
                options = {
                    "maxiter" : steps,
                    "disp"    : True,
                    "gtol"    : 1.e-5,
                }
            )
        except Exception as e:
            return 0, e
        pos, box = cw.flat_to_pos_box(result.x)
        context.setPositions(pos)
        context.setPeriodicBoxVectors(*box)
        state = context.getState(
            getEnergy=True,
            getPositions=True
        )

    return openmm.XmlSerializer.serialize(state), logfile.str


def main():

    """
    Run the main workflow.
    """

    args = parse_arguments()

    HAS_RAY = False
    if args.command == "xml":

        input_dict = {
            "./" : 
                    {
                        "input"       : args.input,
                        "pdb"         : args.pdb,
                        "prefix"      : args.prefix,
                        "steps"       : args.steps,
                        "epsilon"     : args.epsilon,
                        "alternating" : args.alternating,
                        "use_lengths_and_angles" : args.use_lengths_and_angles,
                        "method"      : args.method
                    }
            }

    elif args.command == "yaml":
        import yaml
        with open(args.input, "r") as fopen:
            input_dict = yaml.safe_load(fopen)
        HAS_RAY = True
        try:
            import ray
        except:
            HAS_RAY = False
        if HAS_RAY:
            if "ray_host" in input_dict:
                ray.init(address=input_dict["ray_host"])
            else:
                ray.init()

            ### Wrapper around `run_xtal_min` function for ray
            @ray.remote(num_cpus=input_dict["num_cpus"], num_gpus=1)
            def run_xtal_min_remote(
                xml_path, 
                pdb_path,
                steps = 100,
                epsilon = 1.e-5,
                alternating = False,
                use_lengths_and_angles = False,
                method = "Nelder-Mead",
                platform_name = "CUDA",
                property_dict = {
                    "Precision" : "mixed"
                },
                prefix="xtal_min"):

                return run_xtal_min(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    steps = steps,
                    epsilon = epsilon,
                    alternating = alternating,
                    use_lengths_and_angles  = use_lengths_and_angles,
                    method = method,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix="xtal_min"
                    )

    else:
        raise NotImplementedError(f"Command {args.command} not understood.")

    import os
    import warnings
    worker_id_dict = dict()
    for output_dir in input_dict:
        if output_dir == "num_cpus":
            continue
        if output_dir == "ray_host":
            continue

        if not os.path.exists(input_dict[output_dir]["input"]):
            warnings.warn(f"{input_dict[output_dir]['input']} not found. Skipping.")
            continue
        if not os.path.exists(input_dict[output_dir]["pdb"]):
            warnings.warn(f"{input_dict[output_dir]['pdb']} not found. Skipping.")
            continue

        if HAS_RAY:
            min_func = run_xtal_min_remote.remote
        else:
            min_func = run_xtal_min

        ### Note: CUDA is much faster even for
        ###       the wrapped minimization here (about factor 10)
        worker_id = min_func(
            xml_path = input_dict[output_dir]["input"],
            pdb_path = input_dict[output_dir]["pdb"],
            steps = int(input_dict[output_dir]["steps"]),
            epsilon = float(input_dict[output_dir]["epsilon"]),
            alternating = bool(input_dict[output_dir]["alternating"]),
            use_lengths_and_angles = bool(input_dict[output_dir]["use_lengths_and_angles"]),
            method = input_dict[output_dir]["method"],
            platform_name = "CUDA",
            property_dict = {
                #"Threads"             : '4',
                "DeterministicForces" : "True"
            },
            prefix="xtal_min"
        )
        worker_id_dict[output_dir] = worker_id

    for output_dir in input_dict:
        if output_dir == "num_cpus":
            continue

        if HAS_RAY:
            state, file_str = ray.get(worker_id_dict[output_dir])
        else:
            state, file_str = worker_id_dict[output_dir]
        if state == 0:
            warnings.warn(f"No output for {output_dir}. Minimization failed with {file_str}.")
            continue

        state = openmm.XmlSerializer.deserialize(state)

        os.makedirs(output_dir, exist_ok=True)
        prefix = input_dict[output_dir]["prefix"]
        prefix = f"{output_dir}/{prefix}"

        ### Save state in xml format
        with open(f"{prefix}.xml", "w") as fopen:
            fopen.write(
                openmm.XmlSerializer.serialize(state)
                )
        with open(f"{prefix}.csv", "w") as fopen:
            fopen.write(file_str)
        ### Save State in pdb format
        from openmm import app
        box = state.getPeriodicBoxVectors(asNumpy=True)
        pos = state.getPositions(asNumpy=True)
        pdbfile = app.PDBFile(input_dict[output_dir]["pdb"])
        pdbfile.topology.setPeriodicBoxVectors(box)

        with open(f"{prefix}.pdb", "w") as fopen:
            app.PDBFile.writeHeader(
                pdbfile.topology,
                fopen
                )
            app.PDBFile.writeModel(
                pdbfile.topology,
                pos,
                fopen
                )
            app.PDBFile.writeFooter(
                pdbfile.topology,
                fopen
                )

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()