from setuptools import setup, find_packages

setup(
    name='FMOkit',
    version='0.3.0',
    py_modules=['FMOkit'],
    install_requires=[
        'click',
        'gemmi',
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
        gamoutparser=utils.gamoutparser:cli
    ''',
)
