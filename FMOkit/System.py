import gemmi
from .Fragment import Fragment
from .Atom import Atom
from .hyb_carbon import hybrid_orbitals
from .maeparser import maeparse
from math import dist
import tomllib

ANUMBERS = {"h": 1, "c": 6, "n": 7, "o": 8, "s": 16, "ca": 20,
           "f": 9, "p": 15, "cl": 17}

def coef_format(coef):
    """
    Format the coefficient for FMO.
    :param coef: The coefficient value.
    :return: The formatted coefficient string.
    """
    formatted = ""
    for l in coef.split("\n"):
        els = l.split()
        if len(els) > 0:
            if els[0] in ["0", "1"]:
                formatted +=  f"    {els[0]} {els[1]}" + "".join(f"{float(v):11.6f}" for v in els[2:]) + "\n"
            else:
                formatted += f"       " + "".join(f"{float(v):11.6f}" for v in els) + "\n"
    return formatted.rstrip('\n')

def atom_dist(atom1, atom2):
    """
    Calculate the distance between two atoms.
    :param atom1: The first atom.
    :param atom2: The second atom.
    :return: The distance between the two atoms.
    """
    return dist(
        (atom1.x, atom1.y, atom1.z),
        (atom2.x, atom2.y, atom2.z)
        )

class System:
    def __init__(self, **kwargs):
        self.nodes: int = kwargs["nodes"]
        self.cores: int = kwargs["cores"]
        self.memory: int = kwargs["memory"]
        self.basissets: str = kwargs["basissets"]
        self.fragments: List[Fragment] = []
        self.title: str = "Structure from PDB" # todo: add title to the structure file
        self.fmobnd_list: List[(Atom, Atom)] = []
        self.charge: str = kwargs["charge"]
        self.asym_id: str = kwargs["asym_id"]
        self.pcm: bool = kwargs.get("pcm", False)
        self.runtyp: str = kwargs.get("runtyp", "energy")  # default run type
        self.config: str = kwargs["toml"]

        
        with open(self.config, "rb") as f:
            config = tomllib.load(f)
        self.AAs = config["AminoAcids"]
        self.NTs = config["Nucleotides"]

    def read_file(self, structure_file: str):
        """
        Read the structure file and populate the fragments.
        :param structure_file: The path to the structure file.
        """
        if structure_file.endswith(".cif"):
            self.from_cif(structure_file)
            self.find_fmobnd()
        elif structure_file.endswith(".mae") or structure_file.endswith(".maegz"):
            self.from_mae(structure_file)
            self.find_fmobnd()
        elif structure_file.endswith(".pdb"):
            print("pdb format is not supported")
            raise NotImplementedError("Only cif/mae/maegz formats are supported at the moment.")
        else:
            print("not implemented yet")
            raise NotImplementedError("Only cif/mae/maegz formats are supported at the moment.")
    
    def from_mae(self, structure_file: str):
        """
        Read a Maestro file and populate the fragments.
        :param structure_file: The path to the Maestro file.
        """
        table = maeparse(structure_file)
        
        fragment_dict = {}
        for row in table:
            atom = Atom(
                id=row[0],
                type_symbol=row[1],
                atom_id=row[2],
                x=row[3],
                y=row[4],
                z=row[5],
                charge=row[6]
            )
            
            comp_id = row[7]
            asym_id = row[8]
            seq_id = row[9]

            fragment_key = (comp_id, asym_id, seq_id)
            if fragment_key not in fragment_dict:
                fragment = Fragment(comp_id=comp_id, asym_id=asym_id, seq_id=seq_id)
                fragment_dict[fragment_key] = fragment
                self.fragments.append(fragment)
            else:
                fragment = fragment_dict[fragment_key]

            fragment.atoms.append(atom)
    
    def from_cif(self, structure_file: str):
        """
        Read a CIF file and populate the fragments.
        :param structure_file: The path to the CIF file.
        """
        doc = gemmi.cif.read(structure_file)
        block = doc.sole_block()

        table = block.find("_atom_site.", ["id", "type_symbol", "label_atom_id", "Cartn_x",
             "Cartn_y", "Cartn_z", self.charge, "label_comp_id", self.asym_id, "label_seq_id"])
        if table.width() == 0:
            raise ValueError(f"The specified column is not valid in {structure_file}")
        
        fragment_dict = {}
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

            fragment_key = (comp_id, asym_id, seq_id)
            if fragment_key not in fragment_dict:
                fragment = Fragment(comp_id=comp_id, asym_id=asym_id, seq_id=seq_id)
                fragment_dict[fragment_key] = fragment
                self.fragments.append(fragment)
            else:
                fragment = fragment_dict[fragment_key]

            fragment.atoms.append(atom)

    @property
    def fmoheader(self):
        """
        Generate the header for the system.
        :return: The header string.
        """
        if self.runtyp == "energy":
            if self.basissets == "dftb":
                header = f""" $contrl runtyp={self.runtyp} nprint=-5 ispher=-1 maxit=200 $end
 $system mwords={int(self.memory / (self.cores * 8))} $end
 $gddi ngroup={self.nodes} $end
 $scf dirscf=.true. npunch=0 $end
 $basis gbasis=dftb $end
 $dftb
   scc=.true. dftb3=.true. dampxh=.true. dampex=4.00
   param=3ob-3-1
 $end
 $dft dc=.true. idcver=3 $end
 $pcm solvnt=water ief=-10 icomp=2 icav=1 idisp=1 modpar=73 ifmo=-1 $end
 $pcmcav radii=suahf $end
 $tescav ntsall=60 $end"""
            else:
                header = f""" $contrl runtyp={self.runtyp} nprint=-5 maxit=200 $end
 $system mwords={int(self.memory / (self.cores * 8))} memddi=0 $end
 $gddi ngroup={self.nodes} $end
 $scf dirscf=.t. npunch=0 $end"""
                if self.pcm:
                    header += """\n $pcm solvnt=WATER icomp=2 icav=1 idisp=1 ifmo=-1 mxts=1000000 $end
 $pcmcav radii=vandw $end
 $tescav ntsall=60 mthall=2 $end"""
        elif self.runtyp == "optfmo":
            header = f""" $contrl runtyp=optfmo nprint=-5 ispher=-1 maxit=200 $end
 $system mwords={self.memory} $end
 $gddi ngroup={self.nodes} $end
 $scf dirscf=.true. npunch=0 $end
 $basis gbasis={self.basissets} $end
 $dftb
   scc=.true. dftb3=.true. dampxh=.true. dampex=4.00
   param=3ob-3-1
 $end
 $optfmo
         method=hssupd nstep=2000 opttol=1e-4
 $end
 $pcm solvnt=water ief=-10 icomp=0 icav=1 idisp=1 modpar=65 ifmo=-1 $end
 $pcmcav radii=suahf $end
 $tescav ntsall=60 $end"""

        return header

    @property
    def fmoprp(self):
        """
        Generate the FMO property string.
        :return: The FMO property string.
        """
        if self.runtyp == "energy":
            if self.basissets == "dftb":
                fmoprp = f""" $fmoprp
    modpar=8205
    naodir=210
    ngrfmo(1)={self.nodes}, {self.nodes}, 0, 0, 0,   0, 0, 0, 0, 0
    ipieda=1
    nprint=9
 $end"""                
            else:
                fmoprp = f""" $fmoprp
    ngrfmo(1)={self.nodes}, {self.nodes}, 0, 0, 0,   0, 0, 0, 0, 0
    ipieda=1
    naodir=220
    nprint=9
    maxit=100
 $end"""
        elif self.runtyp == "optfmo":
            fmoprp = f""" $fmoprp
    modpar=8205
    naodir=210
    ngrfmo(1)={self.nodes}, {self.nodes}, 0, 0, 0,   0, 0, 0, 0, 0
    nprint=9
 $end"""

        return fmoprp

    @property
    def icharge(self):
        charges = [f"{f.charge:> d}" for f in self.fragments]
        lines = [
                ("      icharg(1)=" if i == 0 else "                ") +
                ",".join(charges[i:i+10]) +
                ("\n" if i + 10 >= len(charges) else ",\n")
                for i in range(0, len(charges), 10)
                ]
        return ''.join(lines)

    @property
    def fmofragnam(self):
        """
        Generate the FMO fragment names string.
        :return: The FMO fragment names string.
        """
        fragment_names = [f.fragment_name for f in self.fragments]
        lines = [
            ("      frgnam(1)= " if i == 0 else "                 ") +
            ", ".join(fragment_names[i:i+5]) +
            (",\n" if i + 5 < len(fragment_names) else "\n")
            for i in range(0, len(fragment_names), 5)
        ]
        return "".join(lines)

    @property
    def indat(self):
        return "      indat(1)= 0\n" + ''.join(frag.indat for frag in self.fragments)

    @property
    def fmo(self):
        """
        Generate the FMO section string.
        :return: The FMO section string.
        """
        if self.runtyp == "energy":
            if self.basissets == "dftb":
                fmo_str = (
                f" $fmo\n"
                f"      scftyp(1)=rhf\n"
                f"      modgrd=10\n"
                f"      nlayer=1\n"
                f"      nfrag={len(self.fragments)}\n"
                f"{self.icharge}"
                f"{self.fmofragnam}"
                f"{self.indat}"
                f" $end"
                )
            else:
                fmo_str = (
                f" $fmo\n"
                f"      nlayer=1\n"
                f"      mplevl(1)=2\n"
                f"      nfrag={len(self.fragments)}\n"
                f"{self.icharge}"
                f"{self.fmofragnam}"
                f"{self.indat}"
                f" $end"
                )
        elif self.runtyp == "optfmo":
            fmo_str = (
                f" $fmo\n"
                f"      scftyp(1)=rhf\n"
                f"      modgrd=42 modesp=0\n"
                f"      nlayer=1\n"
                f"      nfrag={len(self.fragments)}\n"
                f"{self.icharge}"
                f"{self.fmofragnam}"
                f"{self.indat}"
                f" $end"
                )
        return fmo_str

    @property
    def fmohyb(self):
        """
        Generate the FMO hybrid orbital.
        :return: The FMO hybrid orbitarl string.
        """
        if self.basissets == "dftb":
            fmohyb = """ $fmolmo
 dftb-c 4 4
   1 0  0.558756  0.000000  0.000000  0.829332
   0 1  0.558757  0.781901  0.000000 -0.276445
   0 1  0.558756 -0.390951  0.677146 -0.276445
   0 1  0.558756 -0.390951 -0.677146 -0.276445
 $end"""
        else:
            if self.basissets   not in hybrid_orbitals:
                raise ValueError(f"Basis sets {self.basissets} is not supported.")
            else:
                fmohyb = f""" $fmohyb
  {hybrid_orbitals[self.basissets]["name"]:<10}{hybrid_orbitals[self.basissets]["bda"]:>4}{hybrid_orbitals[self.basissets]["baa"]:>4}
{coef_format(hybrid_orbitals[self.basissets]["coef"])}
  {hybrid_orbitals["MINI"]["name"]:<10}{hybrid_orbitals["MINI"]["bda"]:>4}{hybrid_orbitals["MINI"]["baa"]:>4}
{coef_format(hybrid_orbitals["MINI"]["coef"])}
 $end"""
        return fmohyb

    def find_fmobnd(self):
        ligands = []
        for frg1, frg2 in zip(self.fragments, self.fragments[1:]):
            if frg1.asym_id == frg2.asym_id:
                if (
                    frg1.comp_id in self.AAs and
                    frg2.comp_id in self.AAs and
                    1.2 < atom_dist(frg1.find_atom("C"), frg2.find_atom("N")) < 1.5
                    ):
                    ca_atom, c_atom = frg1.find_atom("CA"), frg1.find_atom("C")
                    self.fmobnd_list.append((ca_atom, c_atom))
                elif frg1.comp_id in self.NTs and frg2.comp_id in self.NTs:
                    # Todo: check for phosphodiester bond
                    cs_atom, c4_atom = frg2.find_atom("C5'"), frg2.find_atom("C4'")
                    self.fmobnd_list.append((cs_atom, c4_atom))
    
            # check for ligands
            if frg2.comp_id not in self.AAs and frg2.comp_id not in self.NTs:
                ligands.append(frg2)

        # check for ligand-ligand interactions
        if len(ligands) > 1:
            for i, lig1 in enumerate(ligands[:-1]):
                for lig2 in ligands[i+1:]:
                    for lig1_c in [lig1_atom for lig1_atom in lig1.atoms if lig1_atom.type_symbol == "C"]:
                        for lig2_c in [lig2_atom for lig2_atom in lig2.atoms if lig2_atom.type_symbol == "C"]:
                            if 1.4 < atom_dist(lig1_c, lig2_c) < 1.6:
                                    self.fmobnd_list.append((lig1_c, lig2_c))

    @property
    def fmobnd(self):
        """
        Generate the FMO bond string.
        :return: The FMO bond string.
        """
        ligands = []
        lines = [" $fmobnd"]
        for atom1, atom2 in self.fmobnd_list:
            if self.basissets == "dftb":
                lines.append(f"{-atom1.id:>10d}{atom2.id:>10d}  {'dftb-c':<10}")
            else:
                lines.append(f"{-atom1.id:>8d}{atom2.id:>6d}  {self.basissets:<10}  {'MINI':<10}")
        return '\n'.join(lines) + "\n $end"

    @property
    def fmodata(self):
        """
        Generate the FMO data string.
        :return: The FMO data string.
        """
        atoms = [atom.type_symbol.lower() for fragment in self.fragments for atom in fragment.atoms]
        atoms = list(set(atoms))
        lines = [f" $data\n {self.title}\n C1"]
        if self.basissets == "dftb":
            for a in atoms:
                if not a in ANUMBERS:
                    print(f"{a} is not in fmodata")
                    exit()                
                lines.append(f" {a}{ANUMBERS[a]:>5}")
        else:
            for a in atoms:
                if not a in ANUMBERS:
                    print(f"{a} is not in fmodata")
                    exit()
                
                lines.append(f" {a}.1-1  {ANUMBERS[a]:>4}")

                if self.basissets == "STO-3G":
                    lines.append("       sto 3\n")
                elif self.basissets == "6-31G":
                    lines.append("       n31 6\n")
                elif self.basissets.startswith("6-31G*"): # 6-31G* or 6-31G**
                    lines.append("       n31 6")                
                    if a in ["h"]:
                        if self.basissets == "6-31G**":
                            lines.append("       d 1\n       1 1.100 1.0\n")
                        else:
                            lines.append("")
                    elif a in ["c", "n", "o", "f"]:
                        lines.append("       d 1\n       1 0.800 1.0\n")
                    elif a in ["cl"]:
                        lines.append("       d 1\n       1 0.750 1.0\n")
                    elif a in ["s"]:
                        lines.append("       d 1\n       1 0.650 1.0\n")
                    elif a in ["p"]:
                        lines.append("       d 1\n       1 0.550 1.0\n")
                    elif a in ["ca"]:
                        lines.append("       d 1\n       1 0.200 1.0\n")
                else:
                    print(f"basis sets({self.basissets}) is not implemented yet")
                    exit()
        lines.append(" $end\n")
        return "\n".join(lines).rstrip('\n')

    @property
    def fmoxyz(self):
        """
        Generate the FMO XYZ format for the system.
        :return: The FMO XYZ format string.
        """
        atoms = sorted(
                    (atom for fragment in self.fragments for atom in fragment.atoms),
                    key=lambda atom: atom.id
                    )
        return  " $fmoxyz\n" + "\n".join(atom.fmoxyz for atom in atoms) + "\n $end"
#        return  " $fmoxyz\n" + "".join(fragment.fmoxyz for fragment in self.fragments) + " $end"

### process fragments

    def prepare_fragments(self):
        self.process_peptide_bond()
        self.process_phosphodiester_bond()
        self.process_cys()
        self.process_nterminal()
        self.process_cterminal()

    def process_peptide_bond(self):
        """
        Process the peptide bond between fragments.
        This method modifies the fragments by moving the C and O atoms from one fragment to the next.
        """
        for frg1, frg2 in zip(self.fragments, self.fragments[1:]):
            if (
                frg1.asym_id == frg2.asym_id and
                frg1.comp_id in self.AAs and
                frg2.comp_id in self.AAs and
                1.2 < atom_dist(frg1.find_atom("C"), frg2.find_atom("N")) < 1.5
            ):
                frg2.atoms.append(frg1.atoms.pop(frg1.find_atom_index("C")))
                frg2.atoms.append(frg1.atoms.pop(frg1.find_atom_index("O")))

    def process_phosphodiester_bond(self):
        """
        Process the phosphodiester bond between fragments.
        This method modifies the fragments by moving the phosphate and oxygen atoms from one fragment to the next.
        """
        for frg1, frg2 in zip(self.fragments, self.fragments[1:]):
            if (
                frg1.asym_id == frg2.asym_id and
                frg1.comp_id in self.NTs and
                frg2.comp_id in self.NTs and
                1.6 < atom_dist(frg1.find_atom("O3'"), frg2.find_atom("P")) < 1.8
            ):
                for atom_id in ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2"]:
                    ai =  frg2.find_atom_index(atom_id)
                    if ai is not None:
                        frg1.atoms.append(frg2.atoms.pop(ai))

    def process_cys(self):
        """
        Process the system by merging fragments and searching for disulfide bonds.
        """
        cys_pairs =self.search_disulfied_bonds()
        for cysname1, cysname2 in cys_pairs:
            self.merge_fragments(cysname1, cysname2, remove_bonds=False)
    
    def process_nterminal(self):
        """
        Process the N-terminal of the system.
        This method modifies the fragments by moving the N-terminal atom to the first fragment.
        """
        for frg in self.fragments:
            if frg.comp_id == "ACE":
                nfrg = self.find_fragment_by_seqid(frg.asym_id, frg.seq_id+1)
                if nfrg is not None:
                    self.merge_fragments(nfrg.fragment_name, frg.fragment_name)

    def process_cterminal(self):
        """
        Process the C-terminal of the system.
        This method modifies the fragments by moving the C-terminal atom to the last fragment.
        """
        for frg in self.fragments:
            if frg.comp_id == "NMA":
                cfrg = self.find_fragment_by_seqid(frg.asym_id, frg.seq_id) # NME/NMA is always at the end
                if cfrg is not None:
                    self.merge_fragments(cfrg.fragment_name, frg.fragment_name)
            elif frg.comp_id == "NME":
                cfrg = self.find_fragment_by_seqid(frg.asym_id, frg.seq_id-1) # NME/NMA is always at the end
                if cfrg is not None:
                    self.merge_fragments(cfrg.fragment_name, frg.fragment_name)
                    
    def merge_fragments(self, frgnam1, frgnam2, remove_bonds=True):
        """
        Merge two fragments by their names.
        :param frgnam1: The name of the first fragment.     
        :param frgnam2: The name of the second fragment.
        """
        frg2 = self.fragments.pop(self.find_index(frgnam2))
        frg1 = self.fragments[self.find_index(frgnam1)]
        if remove_bonds:
            for i, (atom1, atom2) in enumerate(self.fmobnd_list):
                if (atom1 in frg1.atoms and atom2 in frg2.atoms) or (atom1 in frg2.atoms and atom2 in frg1.atoms):
                    self.fmobnd_list.pop(i)
                    break
        # merge the fragments
        for a in frg2.atoms:
            frg1.atoms.append(a)

    
    def find_index(self, fragment_name):
        """
        Find the index of a fragment by its name.
        :param fragment_name: The name of the fragment to find.
        :return: The index of the fragment if found, otherwise None.
        """
        for i, fragment in enumerate(self.fragments):
            if fragment.fragment_name == fragment_name:
                return i
        return None

    def find_fragment_by_seqid(self, asym_id, seq_id):
        """
        Find a fragment by its asym_id and seq_id.
        :param asym_id: The asym_id of the fragment.
        :param seq_id: The seq_id of the fragment.
        :return: The Fragment object if found, otherwise None.
        """
        for fragment in self.fragments:
            if fragment.asym_id == asym_id and fragment.seq_id == seq_id:
                return fragment
        return None

    def search_disulfied_bonds(self):
        cyss = [f for f in self.fragments if f.comp_id == "CYS"]
        sspairs = []
        sspairs.extend(
        (cys1.fragment_name, cys2.fragment_name)
        for i, cys1 in enumerate(cyss[:-1])
        for cys2 in cyss[i+1:]
        if 1.9 < atom_dist(cys1.find_atom("SG"), cys2.find_atom("SG")) < 2.1
        )
        return sspairs
 
    def print_fmoinput(self):
        """
        Generate the FMO input for the system.
        :return: The FMO input string.
        """
        return (
                f"{self.fmoheader}\n"
                f"{self.fmoprp}\n"
                f"{self.fmo}\n"
                f"{self.fmohyb}\n"
                f"{self.fmobnd}\n"
                f"{self.fmodata}\n"
                f"{self.fmoxyz} \n"
                )   