import gemmi
from .Fragment import Fragment
from .Atom import Atom

class System:
    def __init__(self, **kwargs):
        self.nodes: int = kwargs["nodes"]
        self.cores: int = kwargs["cores"]
        self.memory: int = kwargs["memory"]
        self.basissets: str = kwargs["basissets"]
        self.fragments: List[Fragment] = []
    
    def print_header(self):
        """
        Generate the header for the system.
        :return: The header string.
        """
        return f""" $contrl runtyp=energy nprint=-5 maxit=200 $end
 $system mwords={int(self.memory * 1000 / (self.cores * 8))} memddi=0 $end
 $gddi ngroup={self.nodes} $end
 $scf dirscf=.t. npunch=0 $end"""

    def read_file(self, structure_file: str):
        """
        Read the structure file and populate the fragments.
        :param structure_file: The path to the structure file.
        """
        if structure_file.endswith(".cif"):
            self.from_cif(structure_file)
        elif structure_file.endswith(".pdb"):
            print("pdb format is not supported")
            raise NotImplementedError("Only cif format is supported at the moment.")
        else:
            print("not implemented yet")
            raise NotImplementedError("Only cif format is supported at the moment.")
    
    def from_cif(self, structure_file: str):
        """
        Read a CIF file and populate the fragments.
        :param structure_file: The path to the CIF file.
        """
        doc = gemmi.cif.read(structure_file)
        block = doc.sole_block()

        table = block.find("_atom_site.", ["id", "type_symbol", "label_atom_id", "Cartn_x",
             "Cartn_y", "Cartn_z", "partial_charge", "label_comp_id", "label_asym_id", "label_seq_id"])

        fragment = None
        for row in table:
            atom = Atom(
                id=int(row[0]),
                type_symbol=row[1],
                atom_id=row[2],
                x=float(row[3]),
                y=float(row[4]),
                z=float(row[5]),
                charge=float(row[6])
            )
            
            comp_id = row[7]
            asym_id = row[8]
            seq_id = int(row[9])
            if fragment is None or seq_id != fragment.seq_id or asym_id != fragment.asym_id:
                fragment = Fragment(comp_id=comp_id, asym_id=asym_id, seq_id=seq_id)
                self.fragments.append(fragment)
            fragment.atoms.append(atom)

        return table