#!/usr/bin/env python

try:
    from FMOkit import System
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from FMOkit import System

import click

@click.command(help="FMO input generator for GAMESS")
@click.argument("input")
@click.option("--output", "-o", default="fmo.inp", help="FMO inputfile name")
@click.option("--nodes", "-n", default=1, help="num of nodes")
@click.option("--cores", "-c", default=8, help="num of cores")
@click.option("--memory", "-m", default=14000, help="memory in MB")
@click.option("--basissets", "-b", default="6-31G", help="Basis Sets")
@click.option("--charge", "-C", default="partial_charge", help="Basis Sets")
def cli(input, output, nodes, cores, memory, basissets, charge):
    """
    Command line interface for the script.
    """
    mkfmoinp(input, output, nodes, cores, memory, basissets, charge)

def mkfmoinp(input, output, nodes, cores, memory, basissets, charge):
    """
    Generate FMO input file for GAMESS.
    """
    s = System(nodes=nodes, cores=cores, memory=memory, basissets=basissets, charge=charge)
    s.read_file(input)
    s.prepare_fragments()
    with open(output, "w") as f:
        f.write(s.print_fmoinput())
          
if __name__ == "__main__":
    cli()