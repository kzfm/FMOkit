#!/usr/bin/env python

import sys
from sys import argv
import re
import shlex

ANUM2ATOM = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S", 20: "CA",
             9: "F", 15: "P", 17: "CL"}

def maeparse(maefile):
    """
    Parses a MAE file and returns a list of atom records as lists.
    """
    pattern = re.compile(r"m_atom\[(\d+)\].{0,70}{\n(.*?)\n:::\n}", re.DOTALL)
    field_map = {
        "atomic_number": "anum", "pdb_atom_name": "aid",
        "x_coord": "x", "y_coord": "y", "z_coord": "z",
        "charge1": "charge", "pdb_residue_name": "compid",
        "chain_name": "asymid", "residue_number": "seqid"
    }

    with open(maefile, "r") as f:
        content = f.read()

    table = []
    for match in pattern.finditer(content):
        atom_num = int(match.group(1))
        header, data = match.group(2).split("\n:::\n", 1)
        headers = header.split("\n")
        hdict = {v: i for i, h in enumerate(headers) for k, v in field_map.items() if h.endswith(k)}
        lines = data.strip().split("\n")

        def parse_line(line_pair):
            els = sum((shlex.split(l) for l in line_pair), [])
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