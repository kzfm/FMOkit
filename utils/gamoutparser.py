#!/usr/bin/env python

#  Two-body FMO properties.
#    I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ    total     Ees      Eex    Ect+mix    Edisp    Gsol
# ---------------------------------------------------------------------------------------------------------
#    2    1 C1  0   0.00 -0.0199 -9581.007    1.730 -9579.277 -9310.380 -167.795  -73.445  -27.657    0.000
#    3    1 C1 -1   1.27  0.0002   -32.497   -3.140   -35.637   -34.494   -0.008   -0.546   -0.589    0.000
#    3    2 C1  0   0.00 -0.0257 -9582.357    0.764 -9581.593 -9320.896 -155.280  -77.719  -27.698    0.000

import re
import csv
import click

@click.command(help="gamout parser")
@click.argument("gamout")
@click.option("--output", "-o", default="[basename]_pieda.csv", help="Output file name")
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
    if output.startswith("[basename]"):
        basename, suffix = gamout.rsplit(".", 1)
        output = f"{basename}_pieda.csv"
    else:
        basename, suffix = output.rsplit(".")
    # Regular expressions to match the relevant sections in the GAMESS output
    fragment_re = re.compile("CONV\n ={70,90}\n(.*?)\n\n Close fragment pairs", re.DOTALL)
    pieda_re = re.compile(" -{105}\n(.*?)\n\n Total energy", re.DOTALL)

    with open(output, 'w') as f:
        if suffix == "tsv":
            writer = csv.writer(f, delimiter="\t")
        else:
            writer = csv.writer(f)
        writer.writerow(["I", "IFRG", "J", "JFRG", "COMPONENT", "ENERGY", "TOTAL"])
        frgs = {}
        gamout_str = open(gamout, "r").read()

        for m in fragment_re.finditer(gamout_str):
            for l in (m.group(1).split("\n")):
                els = l.split()
                frgs[int(els[0])] = els[1]

        for m in pieda_re.finditer(gamout_str):
            for l in m.group(1).split("\n"):
                fields = {
                    'I': l[:5], 'J': l[5:10], 'DL': l[10:13], 'Z': l[13:16],
                    'R': l[16:23], 'QIJ': l[23:31], 'EIJ': l[31:41], 'dDIJVIJ': l[41:50],
                    'total': l[50:60], 'Ees': l[60:70], 'Eex': l[70:79],
                    'Ectmix': l[79:88], 'Edisp': l[88:97], 'Gsol': l[97:106]
                }

                I, J = int(fields['I']), int(fields['J'])
                components = [("ES", fields['Ees']), ("EX", fields['Eex']), ("CT", fields['Ectmix']),
                              ("DI", fields['Edisp']), ("SOL", fields['Gsol'])]

                for tag, energy in components:
                    for a, b in [(I, J), (J, I)]:
                        writer.writerow([str(a), frgs[a], str(b), frgs[b], tag, energy, fields['total']])   
          
if __name__ == "__main__":
    cli()