#!/usr/bin/env python

import openmm
from openmm import unit
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
        '--replicates', 
        "-r", 
        type=int, 
        help="Number of replicates to generate.", 
        required=False,
        default=10
        )

    xml_parse.add_argument(
        '--resume', 
        "-re", 
        help="Resume md run.", 
        action='store_true',
        default=False,
        required=False,
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


def run_nvt_md(
    xml_path, 
    pdb_path,
    temperature,
    nanoseconds,
    platform_name = "CUDA",
    property_dict = {
        "CudaPrecision" : "mixed"
    },
    prefix = "nvt_md",
    resume = False,
    ):
    
    import os
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters
    
    ### Initialize system things, create integrator, etc...

    ### better not go higher than 1 fs
    time_step = 1. * unit.femtoseconds
    ### The default 1.e-5 seems to be too high and will crash
    ### with the flexible barostat.
    constraint_tolerance = 1.e-6
    
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
    integrator.setConstraintTolerance(constraint_tolerance)

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
    if not resume:
        run_sim = True
    elif resume and os.path.exists(f"{prefix}_thermalization.xml"):
        run_sim = False
    else:
        run_sim = True

    if run_sim:
        ### Initial minimization.
        simulation.minimizeEnergy()

        ### Run for 100 picoseconds
        ### 1 step is 0.001 picoseconds
        ### 1 picosecond is 1000 steps
        try:
            simulation.step(1000 * 100)
        except Exception as e:
            return e
        simulation.saveState(f"{prefix}_thermalization.xml")

    ### 2. Production run
    ### =================
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        1./unit.picoseconds,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(constraint_tolerance)

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
    simulation.loadState(f"{prefix}_thermalization.xml")

    ### Run for 1 nanoseconds (1000 picoseconds)
    ### 1 step is 0.001 picoseconds
    ### 1 picosecond is 1000 steps
    ### 1 nanosecond is 1000000 steps
    write_at = 1000 * 20  # Save every 20 picosends
    N_steps_nanosecond = 1000000
    for i in range(nanoseconds):
        if not resume:
            run_sim = True
        elif resume and os.path.exists(f"{prefix}_production_{i}.xml"):
            run_sim = False
        else:
            run_sim = True
        if run_sim:
            filehandle_dcd = open(f"{prefix}_production_{i}.dcd", "wb")
            filehandle_logger = open(f"{prefix}_production_{i}.csv", "w")
            dcdfile = app.dcdfile.DCDFile(
                filehandle_dcd, 
                topology, 
                time_step, 
                simulation.currentStep,
                write_at
                )
            logfile = Logger(system, filehandle_logger)
            N_iter  = int(N_steps_nanosecond/write_at)
            for _ in range(N_iter):
                try:
                    simulation.step(write_at)
                except Exception as e:
                    filehandle_dcd.close()
                    filehandle_logger.close()
                    return e
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
            with open(f"{prefix}_production_{i}.xml", "w") as fopen:
                simulation.saveState(fopen)

    return 1


def run_xtal_md(
    xml_path, 
    pdb_path,
    temperature,
    nanoseconds,
    platform_name = "CUDA",
    property_dict = {
        "CudaPrecision" : "mixed"
    },
    prefix = "xtal_md",
    resume = False,
    ):
    
    import os
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters
    
    ### Initialize system things, create integrator, etc...

    ### better not go higher than 1 fs
    time_step = 1.0 * unit.femtoseconds
    ### The default 1.e-5 seems to be too high and will crash
    ### with the flexible barostat.
    constraint_tolerance = 1.e-6
    ### friction coefficient for Langevin Thermostat
    friction = 1.0/unit.picoseconds
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()

    ### 1. Temperature equilibration
    ### ============================    
    system = openmm.XmlSerializer.deserialize(xml_str)
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        friction,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(constraint_tolerance)

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
    if not resume:
        run_sim = True
    elif resume and os.path.exists(f"{prefix}_thermalization.xml"):
        run_sim = False
    else:
        run_sim = True

    if run_sim:
        ### Run energy minimization first.
        simulation.minimizeEnergy()
        ### Run for 100 picoseconds
        ### 1 step is 0.001 picoseconds
        ### 1 picosecond is 1000 steps
        try:
            simulation.step(1000 * 100)
        except Exception as e:
            return e
        simulation.saveState(f"{prefix}_thermalization.xml")
    
    ### 2. Pressure equilibration anisotropic
    ### =====================================
    system = openmm.XmlSerializer.deserialize(xml_str)    
    barostat_aniso = openmm.MonteCarloAnisotropicBarostat(
        (1., 1., 1.) * unit.bar,
        temperature.value_in_unit_system(unit.md_unit_system),
    )
    ### Default is 25
    barostat_aniso.setFrequency(25)

    system.addForce(barostat_aniso)
    
    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        friction,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )

    if not resume:
        run_sim = True
    elif resume and os.path.exists(f"{prefix}_pressure1.xml"):
        run_sim = False
    else:
        run_sim = True
    
    if run_sim:
        with open(f"{prefix}_thermalization.xml", "r") as fopen:
            state = openmm.XmlSerializer.deserialize(
                fopen.read()
                )
        simulation.context.setPositions(state.getPositions())
        simulation.context.setVelocities(state.getVelocities())
        simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())

        ### Run for 1 ns
        ### 1 step is 0.001 picoseconds
        ### 1 picosecond is 1000 steps
        try:
            simulation.step(1000 * 1000)
        except Exception as e:
            return e
        simulation.saveState(f"{prefix}_pressure1.xml")
    
    ### 3. Pressure equilibration flexible
    ### ==================================
    system = openmm.XmlSerializer.deserialize(xml_str)
    barostat_aniso = openmm.MonteCarloFlexibleBarostat(
        1. * unit.bar,
        temperature.value_in_unit_system(unit.md_unit_system),
    )
    ### Default is 25
    barostat_aniso.setFrequency(25)
    ### Default is True
    barostat_aniso.setScaleMoleculesAsRigid(True)
    system.addForce(barostat_aniso)

    integrator = openmm.LangevinIntegrator(
        temperature.value_in_unit_system(unit.md_unit_system),
        friction,
        time_step.value_in_unit_system(unit.md_unit_system))
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = openmm.Platform.getPlatformByName(platform_name)
    for property_name, property_value in property_dict.items():
        platform.setPropertyDefaultValue(property_name, property_value)
        
    simulation = app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )

    if not resume:
        run_sim = True
    elif resume and os.path.exists(f"{prefix}_pressure2.xml"):
        run_sim = False
    else:
        run_sim = True
    
    if run_sim:
        with open(f"{prefix}_pressure1.xml", "r") as fopen:
            state = openmm.XmlSerializer.deserialize(
                fopen.read()
                )
        simulation.context.setPositions(state.getPositions())
        simulation.context.setVelocities(state.getVelocities())
        simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        ### Run for 1 ns
        ### 1 step is 0.001 picoseconds
        ### 1 picosecond is 1000 steps
        try:
            simulation.step(1000 * 1000)
        except Exception as e:
            return e
        simulation.saveState(f"{prefix}_pressure2.xml")

    ### 4. Production run
    ### =================
    ### Run for 1 nanoseconds (1000 picoseconds)
    ### 1 step is 0.001 picoseconds
    ### 1 picosecond is 1000 steps
    ### 1 nanosecond is 1000000 steps
    write_at = 1000 * 20  # Save every 20 picosends
    N_steps_nanosecond = 1000000
    for i in range(nanoseconds):
        if not resume:
            run_sim = True
        elif resume and os.path.exists(f"{prefix}_production_{i}.xml"):
            run_sim = False
        else:
            run_sim = True
        if run_sim:
            filehandle_dcd = open(f"{prefix}_production_{i}.dcd", "wb")
            filehandle_logger = open(f"{prefix}_production_{i}.csv", "w")
            dcdfile = app.dcdfile.DCDFile(
                filehandle_dcd, 
                topology, 
                time_step, 
                simulation.currentStep,
                write_at
                )
            logfile = Logger(system, filehandle_logger)
            N_iter  = int(N_steps_nanosecond/write_at)
            for _ in range(N_iter):
                try:
                    simulation.step(write_at)
                except Exception as e:
                    filehandle_dcd.close()
                    filehandle_logger.close()
                    return e
                state = simulation.context.getState(
                    getEnergy=True,
                    getPositions=True,
                    enforcePeriodicBox=True
                )
                dcdfile.writeModel(
                    positions=state.getPositions(),
                    periodicBoxVectors=state.getPeriodicBoxVectors(),
                )
                logfile.write(state)
                simulation.saveState(f"{prefix}_production_latest.xml")

            filehandle_dcd.close()
            filehandle_logger.close()
            with open(f"{prefix}_production_{i}.xml", "w") as fopen:
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
                        "replicates"  : args.replicates,
                        "nvt"         : args.nvt,
                        "resume"      : args.resume,
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
                platform_name = "CUDA",
                property_dict = {
                    "CudaPrecision" : "mixed"
                },
                prefix = "xtal_md",
                resume = False,
                ):
                return run_xtal_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix,
                    resume = resume,
                    )

            @ray.remote(num_cpus=1, num_gpus=1)
            def run_nvt_md_remote(
                xml_path, 
                pdb_path,
                temperature,
                nanoseconds,
                platform_name = "CUDA",
                property_dict = {
                    "CudaPrecision" : "mixed"
                },
                prefix = "nvt_md",
                resume = False,
                ):
                return run_nvt_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix,
                    resume = resume,
                    )

    else:
        raise NotImplementedError(f"Command {args.command} not understood.")

    import os
    import warnings
    worker_id_dict = dict()
    for output_dir in input_dict:
        if output_dir == "ray_host":
            continue
        if not os.path.exists(input_dict[output_dir]["input"]):
            warnings.warn(f"{input_dict[output_dir]['input']} not found. Skipping.")
            continue
        if not os.path.exists(input_dict[output_dir]["pdb"]):
            warnings.warn(f"{input_dict[output_dir]['pdb']} not found. Skipping.")
            continue

        for replicate_idx in range(input_dict[output_dir]["replicates"]):

            output_dir_replicate = f"{output_dir}/run_{replicate_idx:d}"

            os.makedirs(output_dir_replicate, exist_ok=True)
            prefix = input_dict[output_dir]["prefix"]
            prefix = f"{output_dir_replicate}/{prefix}"

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
                platform_name = "CUDA",
                property_dict = {
                    "CudaPrecision" : "mixed"
                },
                prefix = prefix,
                resume = input_dict[output_dir]["resume"],
            )
            worker_id_dict[output_dir_replicate] = worker_id

    ### Only if we have ray, wait for all jobs to finish
    for output_dir_replicate in worker_id_dict:
        if HAS_RAY:
            result = ray.get(worker_id_dict[output_dir_replicate])
        else:
            result = worker_id_dict[output_dir_replicate]
        if result != 1:
            warnings.warn(f"No output for {output_dir_replicate}. MD run failed with {result}.")

def entry_point():

    main()


if __name__ == "__main__":

    entry_point()