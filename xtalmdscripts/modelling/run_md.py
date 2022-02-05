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
        '--nvt', 
        dest='nvt', 
        action='store_true',
        default=False,
        required=False,
        help="Perform md in nvt only."
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
        help="Output prefix for csv and dcd files.", 
        default="xtal_md",
        required=False
        )

    xml_parse.add_argument(
        '--nanoseconds', 
        "-ns", 
        type=int, 
        help="Production length in nanoseconds.", 
        required=False,
        default=100
        )

    xml_parse.add_argument(
        '--temperature', 
        "-t", 
        type=float, 
        help="Target temperature in md run.", 
        required=False,
        default=298.15
        )

    return parser.parse_args()


class Logger(object):
    
    def __init__(self, system, filehandle):
        
        self.dof  = system.getNumParticles() * 3.
        self.dof -= system.getNumConstraints()
        if any(isinstance(force, openmm.CMMotionRemover) for force in system.getForces()):
            self.dof -= 3.
        self.total_mass = 0. * unit.dalton
        for i in range(system.getNumParticles()):
            self.total_mass += system.getParticleMass(i)
            
        self.filehandle = filehandle
        self.write_header()
            
    def write_header(self):
        
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

        self.filehandle.write(f"{time._value},{pot_ene._value},{kin_ene._value},")
        self.filehandle.write(f"{boxvectors[0][0]._value},{boxvectors[0][1]._value},{boxvectors[0][2]._value},")
        self.filehandle.write(f"{boxvectors[1][0]._value},{boxvectors[1][1]._value},{boxvectors[1][2]._value},")
        self.filehandle.write(f"{boxvectors[2][0]._value},{boxvectors[2][1]._value},{boxvectors[2][2]._value},")
        self.filehandle.write(f"{vol._value},{dens._value},{temp._value}")
        self.filehandle.write("\n")


def run_nvt_md(
    xml_path, 
    pdb_path,
    temperature,
    nanoseconds,
    time_step = 2. * unit.femtoseconds,
    platform_name = "CUDA",
    property_dict = {
        "Precision" : "mixed"
    },
    prefix = "xtal_md"
    ):
    
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters
    
    ### Initialize system things, create integrator, etc...
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()
    system = openmm.XmlSerializer.deserialize(xml_str)

    ### 1. Temperature equilibration
    ### ============================    
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )
    
    simulation.context.setPositions(pdbfile.positions)
    if topology.getPeriodicBoxVectors() != None:
        simulation.context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )

    ### Run for 100 picoseconds
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    simulation.step(500 * 100)
    simulation.saveState(f"./{prefix}_thermalization.xml")

    ### 2. Production run
    ### =================
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)

    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )
    
    simulation.context.setPositions(pdbfile.positions)
    if topology.getPeriodicBoxVectors() != None:
        simulation.context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )
    simulation.loadState(f"./{prefix}_thermalization.xml")

    ### Run for 1 nanoseconds (1000 picoseconds)
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    ### 1 nanosecond is 500000 steps
    write_at = 500 * 20  # Save every 20 picosends
    for i in range(nanoseconds):
        filehandle_dcd = open(f"./{prefix}_production_{i}.dcd", "wb")
        filehandle_logger = open(f"./{prefix}_production_{i}.csv", "w")
        dcdfile = app.dcdfile.DCDFile(
            filehandle_dcd, 
            topology, 
            time_step, 
            simulation.currentStep,
            write_at
            )
        logfile = Logger(system, filehandle_logger)
        N_iter  = int(500000/write_at)
        for _ in range(N_iter):
            simulation.step(write_at)
            state = simulation.context.getState(
                getEnergy=True,
                getPositions=True
            )
            dcdfile.writeModel(
                positions=state.getPositions(),
            )
            logfile.write(state)

        filehandle_dcd.close()
        filehandle_logger.close()
        with open(f"./{prefix}_production_{i}.xml", "w") as fopen:
            simulation.saveState(fopen)

    return 1


def run_xtal_md(
    xml_path, 
    pdb_path,
    temperature,
    nanoseconds,
    time_step = 2. * unit.femtoseconds,
    platform_name = "CUDA",
    property_dict = {
        "Precision" : "mixed"
    },
    prefix = "xtal_md"
    ):
    
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters
    
    ### Initialize system things, create integrator, etc...
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()
    system = openmm.XmlSerializer.deserialize(xml_str)

    ### 1. Temperature equilibration
    ### ============================    
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )
    
    simulation.context.setPositions(pdbfile.positions)
    simulation.context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )

    ### Run for 100 picoseconds
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    simulation.step(500 * 100)
    simulation.saveState(f"./{prefix}_thermalization.xml")
    
    ### 2. Pressure equilibration anisotropic
    ### =====================================
    
    barostat_aniso = openmm.MonteCarloAnisotropicBarostat(
        (1., 1., 1.) * unit.bar,
        temperature.value_in_unit_system(unit.md_unit_system),
    )
    ### Default is 25
    barostat_aniso.setFrequency(25)
    
    system.addForce(barostat_aniso)
    
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )
    
    simulation.context.setPositions(pdbfile.positions)
    simulation.context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )
    simulation.loadState(f"./{prefix}_thermalization.xml")
    ### Run for 500 picoseconds
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    simulation.step(500 * 500)
    simulation.saveState(f"./{prefix}_pressure1.xml")
    
    ### 3. Pressure equilibration flexible
    ### ==================================
    
    barostat_aniso = openmm.MonteCarloFlexibleBarostat(
        1. * unit.bar,
        temperature.value_in_unit_system(unit.md_unit_system),
    )
    ### Default is 25
    barostat_aniso.setFrequency(25)
    
    system.addForce(barostat_aniso)
    
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(0.00001)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )
    
    simulation.context.setPositions(pdbfile.positions)
    simulation.context.setPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )
    simulation.loadState(f"./{prefix}_pressure1.xml")
    ### Run for 500 picoseconds
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    simulation.step(500 * 500)
    simulation.saveState(f"./{prefix}_pressure2.xml")
    
    ### 4. Production run
    ### =================
    ### Run for 1 nanoseconds (1000 picoseconds)
    ### 1 step is 0.002 picoseconds
    ### 1 picosecond is 500 steps
    ### 1 nanosecond is 500000 steps
    write_at = 500 * 20  # Save every 20 picosends
    for i in range(nanoseconds):
        filehandle_dcd = open(f"./{prefix}_production_{i}.dcd", "wb")
        filehandle_logger = open(f"./{prefix}_production_{i}.csv", "w")
        dcdfile = app.dcdfile.DCDFile(
            filehandle_dcd, 
            topology, 
            time_step, 
            simulation.currentStep,
            write_at
            )
        logfile = Logger(system, filehandle_logger)
        N_iter  = int(500000/write_at)
        for _ in range(N_iter):
            simulation.step(write_at)
            state = simulation.context.getState(
                getEnergy=True,
                getPositions=True
            )
            dcdfile.writeModel(
                positions=state.getPositions(),
                periodicBoxVectors=state.getPeriodicBoxVectors(),
            )
            logfile.write(state)

        filehandle_dcd.close()
        filehandle_logger.close()
        with open(f"./{prefix}_production_{i}.xml", "w") as fopen:
            simulation.saveState(fopen)

    return 1


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
                        "temperature" : args.temperature,
                        "nanoseconds" : args.nanoseconds,
                        "nvt"         : args.nvt,
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

            ### Wrapper around `run_xtal_md` function for ray
            @ray.remote(num_cpus=1, num_gpus=1)
            def run_xtal_md_remote(
                xml_path, 
                pdb_path,
                temperature,
                nanoseconds,
                time_step = 2. * unit.femtoseconds,
                platform_name = "CUDA",
                property_dict = {
                    "Precision" : "mixed"
                },
                prefix = "xtal_md"
                ):
                return run_xtal_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    time_step = time_step,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix
                    )

            @ray.remote(num_cpus=1, num_gpus=1)
            def run_nvt_md_remote(
                xml_path, 
                pdb_path,
                temperature,
                nanoseconds,
                time_step = 2. * unit.femtoseconds,
                platform_name = "CUDA",
                property_dict = {
                    "Precision" : "mixed"
                },
                prefix = "xtal_md"
                ):
                return run_nvt_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    time_step = time_step,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix
                    )

    else:
        raise NotImplementedError(f"Command {args.command} not understood.")

    import os
    import warnings
    worker_id_dict = dict()
    for output_dir in input_dict:

        if not os.path.exists(input_dict[output_dir]["input"]):
            warnings.warn(f"{input_dict[output_dir]['input']} not found. Skipping.")
            continue
        if not os.path.exists(input_dict[output_dir]["pdb"]):
            warnings.warn(f"{input_dict[output_dir]['pdb']} not found. Skipping.")
            continue

        os.makedirs(output_dir, exist_ok=True)
        prefix = input_dict[output_dir]["prefix"]
        prefix = f"{output_dir}/{prefix}"

        if input_dict[output_dir]["nvt"]:
            if HAS_RAY:
                md_func = run_nvt_md_remote.remote
            else:
                md_func = run_nvt_md
        else:
            if HAS_RAY:
                md_func = run_xtal_md_remote.remote
            else:
                md_func = run_xtal_md

        worker_id = md_func(
            xml_path = input_dict[output_dir]["input"],
            pdb_path = input_dict[output_dir]["pdb"],
            temperature = float(input_dict[output_dir]["temperature"]) * unit.kelvin,
            nanoseconds = int(input_dict[output_dir]["nanoseconds"]),
            time_step = 0.002 * unit.picosecond,
            platform_name = "CUDA",
            property_dict = {
                "Precision" : "mixed"
            },
            prefix = prefix,
        )
        worker_id_dict[output_dir] = worker_id

    if HAS_RAY:
        ### Only if we have ray, wait for all jobs to finish
        for output_dir in worker_id_dict:
            output = ray.get(worker_id_dict[output_dir])


def entry_point():

    main()


if __name__ == "__main__":

    entry_point()