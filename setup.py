from setuptools import setup, find_packages

setup(
    name='FMOkit',
    version='0.4.0',
    py_modules=['FMOkit'],
    install_requires=[
        'click',
        'gemmi',
        # If you are using Maestro or MOE and mmcifprep is not required, 
        # you may safely comment out the following libraries.
        'rdkit',
        'pdbfixer',
        'mdtraj',
        'openmm',
        'openmmforcefields',
        'openff-toolkit'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        mmcifprep=utils.mmcifprep:cli
        mkfmoinp=utils.mkfmoinp:cli
        mkfmodftbinp=utils.mkfmodftbinp:cli
        gamoutparser=utils.gamoutparser:cli
    ''',
)
