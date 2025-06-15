#!/usr/bin/env python

import sys
from sys import argv
import re
import shlex

ANUM2ATOM = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S", 20: "CA",
             9: "F", 15: "P", 17: "CL"}

def maeparse(maefile):
    """
    Parses a MAE file and returns a list of dictionaries containing the data.
    
    Args:
        maefile (str): Path to the MAE file.
        
    Returns:
        list: A list of dictionaries with parsed data.
    """
    atomtable_re = re.compile(r"m_atom\[(\d+)\].{0,70}{\n(.*?)\n:::\n}", re.DOTALL)

    table = []
    with open(maefile, 'r') as file:
        mae_str = open(maefile, "r").read()

    for m in atomtable_re.finditer(mae_str):
        atom_num = int(m.group(1))
        header, data = m.group(2).split("\n:::\n", 1)
        header_list = header.split("\n")
        hdict = {"id":0}
        for i, h in enumerate(header_list):
            #["id", "type_symbol", "label_atom_id", "Cartn_x",
            #"Cartn_y", "Cartn_z", self.charge, "label_comp_id", self.asym_id, "label_seq_id"]
            if h.endswith("atomic_number"):
                hdict["anum"] = i
            elif h.endswith("pdb_atom_name"):
                hdict["aid"] = i
            elif h.endswith("x_coord"):
                hdict["x"] = i
            elif h.endswith("y_coord"):
                hdict["y"] = i
            elif h.endswith("z_coord"):
                hdict["z"] = i
            elif h.endswith("charge1"):
                hdict["charge"] = i
            elif h.endswith("pdb_residue_name"):
                hdict["compid"] = i
            elif h.endswith("chain_name"):
                hdict["asymid"] = i
            elif h.endswith("residue_number"):
                hdict["seqid"] = i
        line_num = len(data.split("\n"))
        if line_num == atom_num:
            for line in data.split("\n"):
                els = shlex.split(line)
                table.append([
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
                ])
        elif line_num == atom_num*2:
            lines = data.split("\n")
            for i in range(0, len(lines), 2):
                els = shlex.split(lines[i])
                els2 = shlex.split(lines[i+1])
                els.extend(els2)
                table.append([
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
                ])
    return table

if __name__ == "__main__":
    maeparse(argv[1])