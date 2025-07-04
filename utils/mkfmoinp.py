#!/usr/bin/env python

try:
    from FMOkit import System
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from FMOkit import System
import os
import click

default_config = os.path.join(os.path.dirname(__file__), "..", "FMOkit", "config", "config.toml")

@click.command(help="FMO input generator for GAMESS")
@click.argument("input")
@click.option("--output", "-o", default="[basename].inp", help="FMO inputfile name")
@click.option("--nodes", "-n", default=1, help="num of nodes")
@click.option("--cores", "-c", default=8, help="num of cores")
@click.option("--memory", "-m", default=14000, help="memory in MB")
@click.option("--basissets", "-b", default="6-31G*", help="Basis Sets")
@click.option("--charge", "-C", default="partial_charge", help="Charge Selection: Partial / Formal")
@click.option("--toml", "-T", default=default_config, help="Path to TOML configuration file")
@click.option("--asym_id", "-A", default="label_asym_id", help="Asym ID Selection: label_asym_id / auth_asym_id")
@click.option("--maestro_cif",  "-M", is_flag=True, help="Use Maestro CIF format (default: mmCIF)")
@click.option("--pcm",  "-P", is_flag=True, help="If True, PCM computation is performed (default: False)")
def cli(input, output, nodes, cores, memory, basissets, charge, asym_id, maestro_cif, pcm, toml):
    """
    Command line interface for the script.
    """
    if maestro_cif:
        s = System(nodes=nodes, cores=cores, memory=memory, basissets=basissets, charge="pdbx_formal_charge", asym_id="auth_asym_id", pcm=pcm, toml=toml)
    else:
        s = System(nodes=nodes, cores=cores, memory=memory, basissets=basissets, charge=charge, asym_id=asym_id, pcm=pcm, toml=toml)
    try:
        s.read_file(input)
    except ValueError as e:
        print(f"Error reading input file: {e}")
        exit(1)
    
    s.prepare_fragments()
    if output.startswith("[basename]"):
        basename, suffix = input.rsplit(".", 1)
        output = f"{basename}.inp"
    with open(output, "w") as f:
        f.write(s.print_fmoinput())
          
if __name__ == "__main__":
    cli()