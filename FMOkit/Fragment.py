class Fragment:
    def __init__(self, **kwargs):
        """
        Initialize a Fragment object with the given parameters.
        :param kwargs: Dictionary containing the parameters for the Fragment.
            - comp_id (str): The component ID of the fragment.
            - asim_id (str): The asymmetric ID of the fragment.
            - seq_id (int): The sequence ID of the fragment.
        """
        self.comp_id: str = kwargs["comp_id"]
        self.asim_id: str = kwargs["asim_id"]
        self.seq_id: int = kwargs["seq_id"]
        self.atoms: List[Atom] = []

    def __repr__(self):
        return f"{self.comp_id}: {self.asim_id} {self.seq_id} ({len(self.atoms)} atoms)"
    
    @property
    def charge(self):
        """
        Calculate the total charge of the fragment.
        :return: The total charge of the fragment.
        """
        return sum(atom.charge for atom in self.atoms)
    
    @property
    def fragment_name(self):
        """
        Get the fragment name in the format "comp_id:seq_id".
        :return: The fragment name.
        """
        return f"{self.comp_id}{self.seq_id:03d}"
    
    @property
    def fmoxyz(self):
        """
        Generate the FMO XYZ format for the fragment.
        :return: The FMO XYZ format string.
        """
        return "\n".join(atom.fmoxyz for atom in self.atoms) + "\n"
    
    @property
    def fmoindat(self):
        """
        1.	Extracts and sorts the seq attribute from all atoms.
	    2.	Compresses consecutive integers into ranges, using a format like:
	        •	1, 2, 3 → 1, -3 (start, negative end of range)
	        •	4 → 4 (single value)
	        •	5, 6 → 5, 6 (just two values are left as-is)
	    3.	Formats the result into lines, with:
	        •	6 values per line, each right-aligned in a 7-character-wide field.
	        •	A leading 10-space indentation ("          ").
    	    •	A trailing 0 added at the end of the final line.
        """
        seqs = sorted(a.id for a in self.atoms) + [10000000]
        ranges = []
        s = t = None

        for i in seqs:
            if s is None:
                s = t = i
            elif i == t + 1:
                t = i
            else:
                ranges.append(s)
                if t > s + 1:
                    ranges.append(-t)
                elif t != s:
                    ranges.append(t)
                s = t = i

        lines = []
        for i in range(0, len(ranges), 6):
            line = ''.join(f"{n:7d}" for n in ranges[i:i+6])
            line += f"{0:7d}" if i + 6 >= len(ranges) else ''
            lines.append("          " + line)

        return '\n'.join(lines) + '\n'