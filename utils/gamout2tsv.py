#!/usr/bin/env python

#  Two-body FMO properties.
#    I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ    total     Ees      Eex    Ect+mix    Edisp    Gsol
# ---------------------------------------------------------------------------------------------------------
#    2    1 C1  0   0.00 -0.0199 -9581.007    1.730 -9579.277 -9310.380 -167.795  -73.445  -27.657    0.000
#    3    1 C1 -1   1.27  0.0002   -32.497   -3.140   -35.637   -34.494   -0.008   -0.546   -0.589    0.000
#    3    2 C1  0   0.00 -0.0257 -9582.357    0.764 -9581.593 -9320.896 -155.280  -77.719  -27.698    0.000

import re
import click

@click.command(help="gamout parser")
@click.argument("gamout")
@click.option("--output", "-o", default="pieda.tsv", help="Output file name")
def cli(gamout, output):
    """
    Command line interface for the script.
    """
    parse_gamout(gamout, output)

def parse_gamout(gamout, output):
    """
    Parse the GAMESS output file and extract fragment information and two-body properties.
    :param gamout: The GAMESS output file.
    :param output: The output file name.
    """
    # Regular expressions to match the relevant sections in the GAMESS output
    fragment_re = re.compile("CONV\n ={70,90}\n(.*?)\n\n Close fragment pairs", re.DOTALL)
    pieda_re = re.compile(" -{105}\n(.*?)\n\n Total energy", re.DOTALL)

    with open(output, "w") as wf:
        wf.write(f"I\tIFRG\tJ\tJFRG\tCOMPONENT\tENERGY\tTOTAL\n") # header
        frgs = {}
        gamout_str = open(gamout, "r").read()
        for m in fragment_re.finditer(gamout_str):
            for l in (m.group(1).split("\n")):
                els = l.split()
                frgs[int(els[0])] = els[1]

        for m in pieda_re.finditer(gamout_str):
            for l in (m.group(1).split("\n")):
                I, J, DL, Z, R, QIJ, EIJ, dDIJVIJ, total, Ees, Eex, Ectmix, Edisp, Gsol = l.split()
                wf.write(f"{I}\t{frgs[int(I)]}\t{J}\t{frgs[int(J)]}\tES\t{Ees}\t{total}\n")
                wf.write(f"{I}\t{frgs[int(I)]}\t{J}\t{frgs[int(J)]}\tEX\t{Eex}\t{total}\n")
                wf.write(f"{I}\t{frgs[int(I)]}\t{J}\t{frgs[int(J)]}\tCT\t{Ectmix}\t{total}\n")
                wf.write(f"{I}\t{frgs[int(I)]}\t{J}\t{frgs[int(J)]}\tDI\t{Edisp}\t{total}\n")
                wf.write(f"{I}\t{frgs[int(I)]}\t{J}\t{frgs[int(J)]}\tSOL\t{Gsol}\t{total}\n")  
          
if __name__ == "__main__":
    cli()