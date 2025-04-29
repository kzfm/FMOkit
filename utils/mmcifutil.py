#!/usr/bin/env python

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import gemmi
import click
import tempfile

@click.command(help="Add partial charges and hydrogen atoms to a mmCIF file")
@click.argument("ifile")
@click.argument("ofile")
@click.option("--minimize", is_flag=True, show_default=True, default=False, help="Minimize the system after adding hydrogens")
def cli(ifile, ofile, minimize):   
    pdb = PDBxFile(ifile)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    residues=modeller.addHydrogens(forcefield)
    system = forcefield.createSystem(
        modeller.topology, 
        nonbondedMethod=PME, 
        nonbondedCutoff=1.0*nanometer, 
        constraints=HBonds)

    ### minimize option does not work so far
    #if minimize:
    #    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    #    simulation = Simulation(modeller.topology, system, integrator)
    #    simulation.context.setPositions(modeller.positions)
    #    simulation.minimizeEnergy()
    
    partial_charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for atom_index in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(atom_index)
                partial_charges.append(str(charge)[:-2])

    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp_file:
        PDBxFile.writeFile(modeller.topology, modeller.positions, tmp_file)

        doc = gemmi.cif.read(tmp_file.name)
        block = doc.sole_block()
        loop = block.find_loop('_atom_site.type_symbol').get_loop()
        loop.add_columns(['_atom_site.partial_charge'], value='?')
        pcs = block.find_values('_atom_site.partial_charge')
        for n, x in enumerate(pcs):
            pcs[n] = partial_charges[n]
        doc.write_file(ofile)

if __name__ == "__main__":
    cli()