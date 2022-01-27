from setuptools import setup, find_packages

__version__ = "0.1"

setup(
    name='xtalmdscripts',
    author='Tobias Huefner',
    author_email='thuefner@health.ucsd.edu',
    description='A collection of scripts useful for the xtalmd project \
in the group of Mike Gilson at UC San Diego.',
    version=__version__,
    license='MIT',
    platforms=['Linux'],
    zip_safe=False,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'make_supercell = xtalmdscripts.supercellbuilding.make_supercell:entry_point',
            'make_tables    = xtalmdscripts.tableformatting.make_tables:entry_point',
            'build_system   = xtalmdscripts.buildsystem.build_system:entry_point',
            ]
        },
    )