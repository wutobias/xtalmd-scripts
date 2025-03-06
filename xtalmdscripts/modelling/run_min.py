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

    xml_parse.add_argument(
        '--positions_only',
        '-po',
        type=str,
        help="Only optimize positions, not box vectors",
        required=False,
        default=False
        )

    xml_parse.add_argument(
        '--platform_name',
        '-pl',
        type=str,
        help="Platform name",
        required=False,
        default="CPU",
        choices=["CPU", "CUDA", "OpenCL"],
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

        ### epsilon is nanometer
        self._epsilon = epsilon
    
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
        
    def fix_box(self, box):

        """
        Fix box vector so that diagonal elements are positive.
        """

        box_cp = np.copy(box[:])
        ### Diagonal elements must be positive.
        box_cp[0,0] = np.abs(box_cp[0,0])
        box_cp[1,1] = np.abs(box_cp[1,1])
        box_cp[2,2] = np.abs(box_cp[2,2])

        return box_cp

    def fix_flat_box(self, boxflat):

        """
        Fix flat box vector so that box lengths are positive.
        """

        boxflat_cp = np.copy(boxflat[:])
        if self.use_lengths_and_angles:
            ### Vector lengths must be positive.
            boxflat_cp[:3] = np.abs(boxflat_cp[:3])
        else:
            ### Diagonal elements must be positive.
            boxflat_cp[0] = np.abs(boxflat_cp[0])
            boxflat_cp[2] = np.abs(boxflat_cp[2])
            boxflat_cp[5] = np.abs(boxflat_cp[5])

        return boxflat_cp

    def box_to_flat(self, box):

        """
        Transform 3x3 box vector to flat vector.
        """

        box_fix = self.fix_box(box)
        box_fix = self.reduce_box(box_fix)
        if self.use_lengths_and_angles:
            return computeLengthsAndAngles(box_fix)

        boxflat      = np.zeros(6, dtype=float)
        boxflat[0]   = box_fix[0,0]
        boxflat[1:3] = box_fix[1,:][:2]
        boxflat[3:]  = box_fix[2,:]

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

        box_fix = self.fix_flat_box(boxflat)
        box_fix = self.reduce_flat_box(box_fix)
        if self.use_lengths_and_angles:
            box = computePeriodicBoxVectors(*box_fix)
            box = box.value_in_unit(unit.nanometer)
        else:
            
            box = np.zeros((3,3), dtype=float)
            box[0,0] = box_fix[0]
            box[1,0] = box_fix[1]
            box[1,1] = box_fix[2]
            box[2,:] = box_fix[3:]

        return box 
    
    def reduce_flat_box(self, boxflat):

        """
        Reduce flat box vectors openmm-style.
        """

        box_fix = self.fix_flat_box(boxflat)
        if self.use_lengths_and_angles:
            box = computePeriodicBoxVectors(*box_fix)
            box  = box.value_in_unit(unit.nanometer)
        else:
            box = np.zeros((3,3), dtype=float)
            box[0,0] = box_fix[0]
            box[1,0] = box_fix[1]
            box[1,1] = box_fix[2]
            box[2,:] = box_fix[3:]

        box     = self.reduce_box(box)
        boxflat = self.box_to_flat(box)

        return boxflat


    def reduce_box(self, box):

        """
        Reduce box vectors openmm-style.
        """

        box = app.internal.unitcell.reducePeriodicBoxVectors(box)
        box = box.value_in_unit(unit.nanometer)
        box = np.array(box)

        return box

    def scale_pos_to_box(self, pos, box_old, box_new):

        """
        Scales coordinates to new box.
        """

        box_old = self.reduce_box(box_old)
        box_new = self.reduce_box(box_new)

        pos     = pos.value_in_unit(unit.nanometer)

        box_old_inv = np.linalg.inv(box_old)
        pos_frac    = np.matmul(box_old_inv, pos.T).T
        pos_new     = np.matmul(box_new, pos_frac.T).T * unit.nanometer

        return pos_new

    def get_frac_pos(self, pos, box):

        """
        Get fractional coordinates of pos in box. pos must be real coordinates.
        """

        box = self.reduce_box(box)
        pos = pos.value_in_unit(unit.nanometer)

        box_inv  = np.linalg.inv(box)
        pos_frac = np.matmul(box_inv, pos.T).T

        return pos_frac

    def get_real_pos(self, pos, box):

        """
        Get real coordinates of pos in box. pos must be fractional coordinates.
        """

        box = self.reduce_box(box)
        pos_real = np.matmul(box, pos.T).T * unit.nanometer

        return pos_real

    def ene_box(self, boxflat):

        """
        Compute system energy from flat boxvector. Energy in kJ/mol.
        Note this also alters the positions.
        """
        
        box_new = self.flatbox_to_box(boxflat)
        box_new = self.reduce_box(box_new)

        state = self.context.getState(
            getPositions=True,
        )

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)

        pos_new = self.scale_pos_to_box(pos_old, box_old, box_new)

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_new
            )
        self.context.reinitialize()
        self.context.setPeriodicBoxVectors(*box_new)
        self.context.setPositions(pos_new)
        state = self.context.getState(
            getEnergy=True,
        )
        ene = state.getPotentialEnergy()

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_old
            )
        self.context.reinitialize()
        self.context.setPositions(pos_old)
        self.context.setPeriodicBoxVectors(*box_old)

        return ene._value
    
    def ene(self, pos_box_flat):

        """
        Compute system energy from flat position and boxvector. Energy in kJ/mol
        """
        
        pos_new, box_new = self.flat_to_pos_box(pos_box_flat)
        box_new = self.reduce_box(box_new)

        state = self.context.getState(
            getPositions=True,
        )

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)
        pos_new = self.scale_pos_to_box(pos_new, box_old, box_new)

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_new
            )
        self.context.reinitialize()
        self.context.setPeriodicBoxVectors(*box_new)
        self.context.setPositions(pos_new)
        state = self.context.getState(
            getEnergy=True,
        )
        ene = state.getPotentialEnergy()

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_old
            )
        self.context.reinitialize()
        self.context.setPeriodicBoxVectors(*box_old)
        self.context.setPositions(pos_old)

        return ene._value
    
    def grad_box(self, boxflat):

        """
        Compute energy gradient from flat box vector. Gradient in kJ/mol/nm
        """
        
        epsilon = self._epsilon
        
        box_new = self.flatbox_to_box(boxflat)
        box_new = self.reduce_box(box_new)

        state = self.context.getState(
            getEnergy=True,
            getPositions=True,
        )

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)

        pos_new  = self.scale_pos_to_box(pos_old, box_old, box_new)
        pos_frac = self.get_frac_pos(pos_new, box_new)

        self.context.setPositions(pos_old)
        self.context.setPeriodicBoxVectors(*box_old)
        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_new
            )
        self.context.reinitialize()
        self.context.setPositions(pos_new)
        self.context.setPeriodicBoxVectors(*box_new)

        boxgrad_flat = np.zeros(6, dtype=float)
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box_new[r,c] += epsilon
                box_plus      = self.reduce_box(box_new)
                pos_new       = self.get_real_pos(pos_frac, box_plus)
                x_new         = self.pos_box_to_flat(
                    pos_new, 
                    box_plus
                    )
                ene_plus  = self.ene(x_new)

                box_new[r,c] -= 2. * epsilon
                box_minus     = self.reduce_box(box_new)
                pos_new       = self.get_real_pos(pos_frac, box_minus)
                x_new         = self.pos_box_to_flat(
                    pos_new, 
                    box_minus
                    )
                ene_minus = self.ene(x_new)

                boxgrad_flat[grad_idx] = (ene_plus - ene_minus)/(epsilon * 2.)

                box_new[r,c] += epsilon

                grad_idx += 1

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_old
            )
        self.context.reinitialize()
        self.context.setPositions(pos_old)
        self.context.setPeriodicBoxVectors(*box_old)

        return boxgrad_flat

    def grad(self, pos_box_flat):

        """
        Compute energy gradient from flat position and boxvector. Gradient in kJ/mol/nm
        """        
        
        epsilon = self._epsilon
        
        pos_new, box_new = self.flat_to_pos_box(pos_box_flat)

        box_new = self.reduce_box(box_new)

        state = self.context.getState(
            getEnergy=True,
            getPositions=True,
        )

        pos_old = state.getPositions(asNumpy=True)
        box_old = state.getPeriodicBoxVectors(asNumpy=True)

        pos_new  = self.scale_pos_to_box(pos_new, box_old, box_new)
        pos_frac = self.get_frac_pos(pos_new, box_new)

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_new
            )
        self.context.reinitialize()
        self.context.setPositions(pos_new)
        self.context.setPeriodicBoxVectors(*box_new)

        boxgrad_flat = np.zeros(6, dtype=float)
        grad_idx = 0
        for r in [0,1,2]:
            ### Note that first vec must be parallel to x axis
            for c in range(r+1):
                box_new[r,c] += epsilon
                box_plus      = self.reduce_box(box_new)
                pos_new       = self.get_real_pos(pos_frac, box_plus)
                x_new         = self.pos_box_to_flat(
                    pos_new, 
                    box_plus
                    )
                ene_plus  = self.ene(x_new)

                box_new[r,c] -= 2. * epsilon
                box_minus     = self.reduce_box(box_new)
                pos_new       = self.get_real_pos(pos_frac, box_minus)
                x_new         = self.pos_box_to_flat(
                    pos_new, 
                    box_minus
                    )
                ene_minus = self.ene(x_new)

                boxgrad_flat[grad_idx] = (ene_plus - ene_minus)/(epsilon * 2.)

                box_new[r,c] += epsilon

                grad_idx += 1

        state   = self.context.getState(getForces=True)
        forces  = state.getForces(asNumpy=True)._value
        forces *= -1

        grad    = np.hstack((forces.flatten(), boxgrad_flat))

        system = self.context.getSystem()
        system.setDefaultPeriodicBoxVectors(
            *box_old
            )
        self.context.reinitialize()
        self.context.setPositions(pos_old)
        self.context.setPeriodicBoxVectors(*box_old)

        return grad


def run_xtal_min(
    xml_path, 
    pdb_path,
    steps = 100,
    epsilon = 1.e-5,
    alternating = False,
    use_lengths_and_angles = False,
    method = "Nelder-Mead",
    positions_only = False,
    platform_name = "CPU",
    property_dict = {
        "Threads" : "4"
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

    try:
        with open(xml_path, "r") as fopen:
            xml_str = fopen.read()
        system = openmm.XmlSerializer.deserialize(xml_str)
        system.setDefaultPeriodicBoxVectors(
            *topology.getPeriodicBoxVectors()
            )

        integrator = openmm.LangevinMiddleIntegrator(
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
    except Exception as e:
        return 0, e

    cw = ContextWrapper(
        context, 
        system.getNumParticles(),
        epsilon,
        use_lengths_and_angles
        )
    logfile = Logger(system, "")

    ### Only atomic positions
    if positions_only:
        try:
            openmm.LocalEnergyMinimizer.minimize(context)
            state = context.getState(
                    getEnergy=True,
                    getPositions=True)
            logfile.write(state)
            logfile.write_footer(
                    f"# Final Energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)} kJ/mol")
        except Exception as e:
            return 0, e

        return openmm.XmlSerializer.serialize(state), logfile.str

    ### Alternate between openmm native minimizer
    ### and xtal cell minimization.
    if alternating:
        try:
            best_ene = 999999999999999999999.
            best_x   = None
            for _ in range(steps):
                openmm.LocalEnergyMinimizer.minimize(context)
                pos_old, box_old = cw.flat_to_pos_box(
                    cw.pos_box_flat
                    )
                x0   = np.copy(cw.box_flat)
                ene1 = cw.ene_box(x0)
                result = optimize.minimize(
                    fun=cw.ene_box,
                    x0=x0, 
                    method=method,
                    jac=cw.grad_box,
                    hess='3-point',
                    options = {
                        "maxiter" : steps,
                        "disp"    : False,
                        "gtol"    : 1.e-5,
                    }

                )
                box_new = cw.flatbox_to_box(result.x)
                pos_new = cw.scale_pos_to_box(
                    pos_old, 
                    box_old, 
                    box_new
                    )

                system = context.getSystem()
                system.setDefaultPeriodicBoxVectors(
                    *box_new
                    )
                context.reinitialize()
                context.setPeriodicBoxVectors(*box_new)
                context.setPositions(pos_new)
                state = context.getState(
                    getEnergy=True,
                    getPositions=True
                )
                ene = state.getPotentialEnergy()._value
                if ene < best_ene:
                    best_ene = ene
                    best_x   = cw.pos_box_to_flat(pos_new, box_new)
                else:
                    ### We are done.
                    break
                logfile.write(state)

            logfile.write_footer(
                f"# Final Energy: {result.fun} kJ/mol"
                )

            pos, box = cw.flat_to_pos_box(best_x)
            system = context.getSystem()
            system.setDefaultPeriodicBoxVectors(
                *box
                )
            context.reinitialize()
            context.setPeriodicBoxVectors(*box)
            context.setPositions(pos)

            state = context.getState(
                getEnergy=True,
                getPositions=True
            )
        except Exception as e:
            return 0, e

    ### Minimize everything together.
    else:
        try:
            def callback(xk):
                state = context.getState(
                    getPositions=True
                )
                pos_old = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                box_old = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
                pos_new, box_new = cw.flat_to_pos_box(xk)
                system = context.getSystem()
                system.setDefaultPeriodicBoxVectors(
                    *box_new
                    )
                context.reinitialize()
                context.setPositions(pos_new)
                context.setPeriodicBoxVectors(*box_new)
                state = context.getState(
                    getEnergy=True,
                    getPositions=True
                )
                logfile.write(state)

                system = context.getSystem()
                system.setDefaultPeriodicBoxVectors(
                    *box_old
                    )
                context.reinitialize()
                context.setPositions(pos_old)
                context.setPeriodicBoxVectors(*box_old)

            openmm.LocalEnergyMinimizer.minimize(context)
            x0 = cw.pos_box_flat
            result = optimize.minimize(
                fun=cw.ene,
                x0=x0, 
                method=method,
                jac=cw.grad,
                hess='3-point',
                callback=callback,
                options = {
                    "maxiter" : steps,
                    "disp"    : True,
                    "gtol"    : 1.e-5,
                }
            )
            logfile.write_footer(
                f"# Final Energy: {result.fun} kJ/mol"
                )
            pos, box = cw.flat_to_pos_box(result.x)
            system = context.getSystem()
            system.setDefaultPeriodicBoxVectors(
                *box
                )
            context.reinitialize()
            context.setPositions(pos)
            context.setPeriodicBoxVectors(*box)
            state = context.getState(
                getEnergy=True,
                getPositions=True
            )
        except Exception as e:
            return 0, e

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
                        "input"                  : args.input,
                        "pdb"                    : args.pdb,
                        "prefix"                 : args.prefix,
                        "steps"                  : args.steps,
                        "epsilon"                : args.epsilon,
                        "alternating"            : args.alternating,
                        "use_lengths_and_angles" : args.use_lengths_and_angles,
                        "method"                 : args.method,
                        "positions_only"         : args.positions_only,
                        "platform_name"          : args.platform_name,
                    }
            }

        input_dict["num_cpus"] = 4

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

            if "platform_name" in input_dict:
                platform_name = input_dict["platform_name"]
                if platform_name == "CPU":
                    ray_dict = {
                            "num_cpus" : input_dict["num_cpus"],
                            "num_gpus" : 0
                            }
                elif platform_name in ["CUDA", "OpenCL"]:
                    ray_dict = {
                            "num_cpus" : 1,
                            "num_gpus" : 1
                            }
                else:
                    import warnings
                    warnings.warn(
                            f"Platform name {platform_name} not known.")
                    ray_dict = {}
            else:
                ray_dict = {}

            ### Wrapper around `run_xtal_min` function for ray
            @ray.remote(**ray_dict)
            def run_xtal_min_remote(
                xml_path, 
                pdb_path,
                steps = 100,
                epsilon = 1.e-5,
                alternating = False,
                use_lengths_and_angles = False,
                method = "Nelder-Mead",
                positions_only = False,
                platform_name = "CPU",
                property_dict = {"Threads" : "4"},
                prefix = "xtal_min"):

                return run_xtal_min(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    steps = steps,
                    epsilon = epsilon,
                    alternating = alternating,
                    use_lengths_and_angles  = use_lengths_and_angles,
                    method = method,
                    positions_only = positions_only,
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
        if output_dir == "prefix":
            continue
        if output_dir == "steps":
            continue
        if output_dir == "method":
            continue
        if output_dir == "alternating":
            continue
        if output_dir == "use_lengths_and_angles":
            continue
        if output_dir == "epsilon":
            continue
        if output_dir == "positions_only":
            continue
        if output_dir == "platform_name":
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

        ### Look for global options
        if "prefix" in input_dict[output_dir]:
            prefix = input_dict[output_dir]["prefix"]
        elif "prefix" in input_dict:
            prefix = input_dict["prefix"]
        else:
            prefix = args.prefix

        if os.path.exists(f"{output_dir}/{prefix}.xml"):
            warnings.warn(f"{output_dir}/{prefix}.xml already exists. Skipping.")
            continue
        if os.path.exists(f"{output_dir}/{prefix}.pdb"):
            warnings.warn(f"{output_dir}/{prefix}.pdb already exists. Skipping.")
            continue
        if os.path.exists(f"{output_dir}/{prefix}.csv"):
            warnings.warn(f"{output_dir}/{prefix}.csv already exists. Skipping.")
            continue

        if "steps" in input_dict[output_dir]:
            steps = input_dict[output_dir]["steps"]
        elif "steps" in input_dict:
            steps = input_dict["steps"]
        else:
            steps = args.steps

        if "method" in input_dict[output_dir]:
            method = input_dict[output_dir]["method"]
        elif "method" in input_dict:
            method = input_dict["method"]
        else:
            method = args.method

        if "alternating" in input_dict[output_dir]:
            alternating = input_dict[output_dir]["alternating"]
        elif "alternating" in input_dict:
            alternating = input_dict["alternating"]
        else:
            alternating = args.alternating

        if "use_lengths_and_angles" in input_dict[output_dir]:
            use_lengths_and_angles = input_dict[output_dir]["use_lengths_and_angles"]
        elif "use_lengths_and_angles" in input_dict:
            use_lengths_and_angles = input_dict["use_lengths_and_angles"]
        else:
            use_lengths_and_angles = args.use_lengths_and_angles

        if "epsilon" in input_dict[output_dir]:
            epsilon = input_dict[output_dir]["epsilon"]
        elif "epsilon" in input_dict:
            epsilon = input_dict["epsilon"]
        else:
            epsilon = args.epsilon

        if "positions_only" in input_dict[output_dir]:
            positions_only = input_dict[output_dir]["positions_only"]
        elif "positions_only" in input_dict:
            positions_only = input_dict["positions_only"]
        else:
            positions_only = args.positions_only

        if "platform_name" in input_dict[output_dir]:
            platform_name = input_dict[output_dir]["platform_name"]
        elif "platform_name" in input_dict:
            platform_name = input_dict["platform_name"]
        else:
            platform_name = args.platform_name

        if platform_name == "CPU":
            property_dict = {"Threads" : "4"}
        elif platform_name in ["CUDA", "OpenCL"]:
            property_dict = {}
        else:
            property_dict = {}

        worker_id = min_func(
            xml_path = input_dict[output_dir]["input"],
            pdb_path = input_dict[output_dir]["pdb"],
            steps = int(steps),
            epsilon = float(epsilon),
            alternating = bool(alternating),
            use_lengths_and_angles = bool(use_lengths_and_angles),
            method = method,
            positions_only = positions_only,
            platform_name = platform_name,
            property_dict = property_dict,
            prefix=prefix
        )
        worker_id_dict[worker_id] = output_dir

    worker_id_list = list(worker_id_dict.keys())
    while worker_id_list:
    
        if HAS_RAY:
            [worker_id], worker_id_list = ray.wait(worker_id_list)
            state, file_str = ray.get(worker_id)
        else:
            worker_id       = worker_id_list.pop(0)
            state, file_str = worker_id
        if state == 0:
            warnings.warn(f"No output for {output_dir}. Minimization failed with {file_str}.")
            continue
            
        output_dir = worker_id_dict[worker_id]

        state = openmm.XmlSerializer.deserialize(state)

        os.makedirs(output_dir, exist_ok=True)
        if "prefix" in input_dict[output_dir]:
            prefix = input_dict[output_dir]["prefix"]
        elif "prefix" in input_dict:
            prefix = input_dict["prefix"]
        else:
            prefix = "out"
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
