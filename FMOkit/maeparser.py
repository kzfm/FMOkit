#!/usr/bin/env python

import sys
from sys import argv
import re
import shlex
import gzip

ANUM2ATOM = {
    1:  "H",
    3:  "LI",
    5:  "B",
    6:  "C",
    7:  "N",
    8:  "O",
    9:  "F",
    11: "NA",
    12: "MG",
    14: "SI",
    15: "P",
    16: "S",
    17: "CL",
    19: "K",
    20: "CA",
    26: "FE",
    30: "ZN",
    35: "BR",
    53: "I"
}

def maeparse(maefile):
    """
    Parses a MAE file and returns a list of atom records as lists.
    """
    pattern = re.compile(r"m_atom\[(\d+)\].{0,70}{ {0,5}\n(.*?)\n {0,5}:::\n {0,5}}", re.DOTALL)

    field_map = {
        "atomic_number": "anum", "pdb_atom_name": "aid",
        "m_x_coord": "x", "m_y_coord": "y", "m_z_coord": "z",
        "charge1": "charge", "pdb_residue_name": "compid",
        "chain_name": "asymid", "residue_number": "seqid"
    }
    if maefile.endswith("maegz"):
        with gzip.open(maefile, "rt") as f:
            content = f.read()
    elif maefile.endswith("mae"):
        with open(maefile, "r") as f:
            content = f.read()

    table = []
    for match in pattern.finditer(content):
        atom_num = int(match.group(1))
        header, data = match.group(2).split(":::", 1)
        header = header.strip()
        headers = header.split("\n")
        hdict = {v: i for i, h in enumerate(headers) for k, v in field_map.items() if h.endswith(k)}
        lines = data.strip().split("\n")

        def parse_line(line_pair):
            els = sum((shlex.split(l) for l in line_pair), [])
            if not int(els[hdict["anum"]]) in ANUM2ATOM:
                    print(f"Atom number:{els[hdict["anum"]]} is not in fmodata")
                    exit()

            return [
                int(els[0]),
                ANUM2ATOM[int(els[hdict["anum"]])],
                els[hdict["aid"]].strip(),
                float(els[hdict["x"]]),
                float(els[hdict["y"]]),
                float(els[hdict["z"]]),
                float(els[hdict["charge"]]),
                els[hdict["compid"]].strip(),
                els[hdict["asymid"]],
                int(els[hdict["seqid"]])
            ]

        if len(lines) == atom_num:
            for line in lines:
                table.append(parse_line([line]))
        elif len(lines) == atom_num * 2:
            for i in range(0, len(lines), 2):
                table.append(parse_line(lines[i:i+2]))

    return table

if __name__ == "__main__":
    maeparse(argv[1])