#!/usr/bin/env python

import openmm
from openmm import unit
import numpy as np
from scipy import optimize
from .utils import Logger

DEBUG = False

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
        help="Perform md in nvt."
        )

    xml_parse.add_argument(
        '--restart', 
        dest='restart', 
        action='store_true',
        default=False,
        required=False,
        help="Try to restart md simulation from previously saved xml files."
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
        '--temperature', 
        "-t", 
        type=float, 
        help="Target temperature in md run.", 
        required=False,
        default=298.15
        )

    return parser.parse_args()


def add_CMMotionRemover(system):

    __doc__ = """
    Add CMMotionRemover to system *in place*.
    """

    import openmm

    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    if not "CMMotionRemover" in forces:
        cmm = openmm.CMMotionRemover()
        system.addForce(cmm)

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
    restart = False,
    ):
    
    import os
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters
    
    ### Initialize system things, create integrator, etc...

    ### better not go higher than 1 fs
    time_step = 1 * unit.femtoseconds
    ### The default 1.e-5 seems to be too high and will crash
    ### with the flexible barostat.
    constraint_tolerance = 1.e-6

    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()
    system = openmm.XmlSerializer.deserialize(xml_str)
    add_CMMotionRemover(system)

    ### 1. Temperature equilibration
    ### ============================    
    integrator = openmm.LangevinMiddleIntegrator(
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
    ### Run for 100 picoseconds
    ### 1 step is 0.001 picoseconds
    ### 1 picosecond is 1000 steps
    if os.path.exists(f"{prefix}_thermalization.xml") and restart:
        simulation.loadState(f"{prefix}_thermalization.xml")
    else:
        restart = False
        try:
            simulation.step(1000 * 100)
        except Exception as e:
            return e
    simulation.saveState(
        f"{prefix}_thermalization.xml"
        )

    ### 2. Production run
    ### =================
    integrator = openmm.LangevinMiddleIntegrator(
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
        if os.path.exists(f"{prefix}_production_{i}.xml") and restart:
            simulation.loadState(f"{prefix}_production_{i}.xml")
        else:
            restart = False
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
                    enforcePeriodicBox=False
                )
                dcdfile.writeModel(
                    positions=state.getPositions(),
                )
                logfile.write(state)

            filehandle_dcd.close()
            filehandle_logger.close()
        simulation.saveState(f"{prefix}_production_{i}.xml")

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
    restart = False,
    ):

    import os
    import openmm
    from openmm import unit
    from openmm import app
    import numpy as np
    from mdtraj import reporters

    ### Initialize system things, create integrator, etc...

    ### better not go higher than 1 fs
    time_step = 1 * unit.femtoseconds
    ### The default 1.e-5 seems to be too high and will crash
    ### with the flexible barostat.
    constraint_tolerance = 1.e-6
    ### friction coefficient for Langevin Thermostat
    friction = 1./unit.picoseconds
    
    pdbfile  = app.PDBFile(pdb_path)
    topology = pdbfile.topology

    with open(xml_path, "r") as fopen:
        xml_str = fopen.read()

    ### Pre OpenMM 8.1.2, we have to initialize
    ### the system with a triclinc box. Otherwise
    ### the MIC distance calculations are wrong.
    tric_box = [
            openmm.Vec3(3,0,0)*unit.nanometer,
            openmm.Vec3(0,3,0)*unit.nanometer,
            openmm.Vec3(0,1,3)*unit.nanometer]

    ### 1. Temperature equilibration
    ### ============================
    system = openmm.XmlSerializer.deserialize(xml_str)
    add_CMMotionRemover(system)
    system.setDefaultPeriodicBoxVectors(
        *tric_box)

    integrator = openmm.LangevinMiddleIntegrator(
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
    
    simulation.context.setPositions(
        pdbfile.positions
        )
    simulation.context.setPeriodicBoxVectors(
        *topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(
            temperature.value_in_unit_system(
                unit.md_unit_system
                )
            )

    ### Run for 100 picoseconds
    ### 1 step is 0.001 picoseconds
    ### 1 picosecond is 1000 steps
    if os.path.exists(f"{prefix}_thermalization.xml") and restart:
        simulation.loadState(f"{prefix}_thermalization.xml")
    else:
        restart = False
        try:
            simulation.step(1000 * 100)
        except Exception as e:
            return e
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        getForces=True,
        enforcePeriodicBox=False
        )
    simulation.saveState(
        f"{prefix}_thermalization.xml"
        )
    topology.setPeriodicBoxVectors(
        state.getPeriodicBoxVectors()
        )
    with open(f"{prefix}_thermalization.pdb", "w") as fopen:
        app.PDBFile.writeHeader(
            topology,
            fopen
            )
        app.PDBFile.writeModel(
            topology,
            state.getPositions(),
            fopen
            )
        app.PDBFile.writeFooter(
            topology,
            fopen
            )
    
    ### 2. Pressure equilibration anisotropic
    ### =====================================
    if os.path.exists(f"{prefix}_pressure1.xml") and restart:
        system = openmm.XmlSerializer.deserialize(xml_str)
        add_CMMotionRemover(system)
        system.setDefaultPeriodicBoxVectors(
            *tric_box)

        barostat_aniso = openmm.MonteCarloAnisotropicBarostat(
            (1., 1., 1.) * unit.bar,
            temperature.value_in_unit_system(unit.md_unit_system),
        )
        ### Default is 25
        barostat_aniso.setFrequency(25)
        system.addForce(barostat_aniso)

        integrator = openmm.LangevinMiddleIntegrator(
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
        simulation.loadState(f"{prefix}_pressure1.xml")
    else:
        restart = False
        for _ in range(100):
            system = openmm.XmlSerializer.deserialize(xml_str)
            add_CMMotionRemover(system)
            system.setDefaultPeriodicBoxVectors(
                *tric_box)

            barostat_aniso = openmm.MonteCarloAnisotropicBarostat(
                (1., 1., 1.) * unit.bar,
                temperature.value_in_unit_system(unit.md_unit_system),
            )
            ### Default is 25
            barostat_aniso.setFrequency(25)
            system.addForce(barostat_aniso)

            integrator = openmm.LangevinMiddleIntegrator(
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

            simulation.context.setPositions(
                state.getPositions()
                )
            simulation.context.setVelocities(
                state.getVelocities()
                )
            simulation.context.setPeriodicBoxVectors(
                *state.getPeriodicBoxVectors()
                )
            ### Run for 0.01 ns
            ### 1 step is 0.001 picoseconds
            ### 1 picosecond is 1000 steps
            try:
                simulation.step(1000 * 10)
            except Exception as e:
                return e
            state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getEnergy=True,
                getForces=True,
                enforcePeriodicBox=False
            )
            topology.setPeriodicBoxVectors(
                state.getPeriodicBoxVectors()
                )
    simulation.saveState(
        f"{prefix}_pressure1.xml"
        )
    with open(f"{prefix}_pressure1.pdb", "w") as fopen:
        state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getEnergy=True,
            getForces=True,
            enforcePeriodicBox=False
            )
        app.PDBFile.writeHeader(
            topology,
            fopen
            )
        app.PDBFile.writeModel(
            topology,
            state.getPositions(),
            fopen
            )
        app.PDBFile.writeFooter(
            topology,
            fopen
            )
    
    ### 3. Pressure equilibration flexible
    ### ==================================
    if os.path.exists(f"{prefix}_pressure2.xml") and restart:
        system = openmm.XmlSerializer.deserialize(xml_str)
        add_CMMotionRemover(system)
        system.setDefaultPeriodicBoxVectors(
            *tric_box)

        barostat_aniso = openmm.MonteCarloFlexibleBarostat(
            1.0 * unit.bar,
            temperature.value_in_unit_system(unit.md_unit_system),
        )

        ### Default is 25
        barostat_aniso.setFrequency(25)
        ### Default is True
        barostat_aniso.setScaleMoleculesAsRigid(False)
        system.addForce(barostat_aniso)

        integrator = openmm.LangevinMiddleIntegrator(
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
        simulation.loadState(f"{prefix}_pressure2.xml")
    else:
        restart = False
        for _ in range(100):
            system = openmm.XmlSerializer.deserialize(xml_str)
            add_CMMotionRemover(system)
            system.setDefaultPeriodicBoxVectors(
                *tric_box
                )

            barostat_aniso = openmm.MonteCarloFlexibleBarostat(
                1.0 * unit.bar,
                temperature.value_in_unit_system(unit.md_unit_system),
            )

            ### Default is 25
            barostat_aniso.setFrequency(25)
            ### Default is True
            barostat_aniso.setScaleMoleculesAsRigid(False)
            system.addForce(barostat_aniso)

            integrator = openmm.LangevinMiddleIntegrator(
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

            simulation.context.setPositions(
                state.getPositions()
                )
            simulation.context.setVelocities(
                state.getVelocities()
                )
            simulation.context.setPeriodicBoxVectors(
                *state.getPeriodicBoxVectors()
                )
            ### Run for 0.01 ns
            ### 1 step is 0.001 picoseconds
            ### 1 picosecond is 1000 steps
            try:
                simulation.step(1000 * 10)
            except Exception as e:
                return e
            state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getEnergy=True,
                getForces=True,
                enforcePeriodicBox=False
            )
            topology.setPeriodicBoxVectors(
                state.getPeriodicBoxVectors()
                )
    simulation.saveState(
        f"{prefix}_pressure2.xml"
        )
    with open(f"{prefix}_pressure2.pdb", "w") as fopen:
        state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getEnergy=True,
            getForces=True,
            enforcePeriodicBox=False
            )
        app.PDBFile.writeHeader(
            topology,
            fopen
            )
        app.PDBFile.writeModel(
            topology,
            state.getPositions(),
            fopen
            )
        app.PDBFile.writeFooter(
            topology,
            fopen
            )

    ### 4. Production run
    ### =================
    ### Run for 1 nanoseconds (1000 picoseconds)
    ### 1 step is 0.001 picoseconds
    ### 1 picosecond is 1000 steps
    ### 1 nanosecond is 1000000 steps
    write_at = 1000 * 20  # Save every 20 picosends
    N_steps_nanosecond = 1000000
    for i in range(nanoseconds):
        if os.path.exists(f"{prefix}_production_{i}.xml") and restart:
            simulation.loadState(
                f"{prefix}_production_{i}.xml"
                )
        else:
            restart = False
            if DEBUG:
                filehandle_dcd_vel = open(f"{prefix}_production_v_{i}.dcd", "wb")
                dcdfile_vel = app.dcdfile.DCDFile(
                    filehandle_dcd_vel, 
                    topology, 
                    time_step, 
                    simulation.currentStep,
                    write_at
                    )
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
                    getPositions=True,
                    getVelocities=True,
                    getEnergy=True,
                    getForces=True,
                    enforcePeriodicBox=True
                )
                dcdfile.writeModel(
                    positions=state.getPositions(),
                    periodicBoxVectors=state.getPeriodicBoxVectors(),
                )
                if DEBUG:
                    vel  = state.getVelocities()._value
                    vel *= unit.nanometer
                    dcdfile_vel.writeModel(
                        positions=vel,
                        )
                logfile.write(state)
                simulation.saveState(
                    f"{prefix}_production_latest.xml"
                    )

            filehandle_dcd.close()
            if DEBUG:
                filehandle_dcd_vel.close()
            filehandle_logger.close()
        simulation.saveState(
            f"{prefix}_production_{i}.xml"
            )

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
                        "restart"     : args.restart,
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
                restart = False,
                ):
                return run_xtal_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix,
                    restart = restart,
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
                restart = False,
                ):
                return run_nvt_md(
                    xml_path = xml_path, 
                    pdb_path = pdb_path,
                    temperature = temperature,
                    nanoseconds = nanoseconds,
                    platform_name = platform_name,
                    property_dict = property_dict,
                    prefix = prefix,
                    restart = restart,
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
                #platform_name = "CPU",
                #property_dict = {"Threads" : "2"},
                prefix = prefix,
                restart = bool(input_dict[output_dir]["restart"])
            )
            worker_id_dict[output_dir_replicate] = worker_id

    ### Only if we have ray, wait for all jobs to finish
    for output_dir_replicate in worker_id_dict:
        if not output_dir_replicate in worker_id_dict:
            continue
        if HAS_RAY:
            try:
                result = ray.get(worker_id_dict[output_dir_replicate])
            except Exception as e:
                result = e
        else:
            result = worker_id_dict[output_dir_replicate]
        if result != 1:
            warnings.warn(f"No output for {output_dir_replicate}. MD run failed with {result}.")


def entry_point():

    main()


if __name__ == "__main__":

    entry_point()
