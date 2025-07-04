# Ref: https://pdbj.org/cms-data/workshop/20200916/mmCIF_20200915.pdf
class Atom:
    def __init__(self, **kwargs):
        """
        Initialize an Atom object with the given parameters.
        :param kwargs: Dictionary containing the parameters for the Atom.
            - id (int): The ID of the atom.
            - type_symbol (str): The type symbol of the atom.
            - atom_id (str): The ID of the atom in the system.
            - charge (float): The charge of the atom.
            - x (float): The x-coordinate of the atom.
            - y (float): The y-coordinate of the atom.
            - z (float): The z-coordinate of the atom.
        """
        self.id: int = kwargs["id"]
        self.type_symbol: str = kwargs["type_symbol"]
        self.atom_id: str = kwargs["atom_id"]
        self.charge: float = kwargs["charge"]
        self.x: float = kwargs["x"]
        self.y: float = kwargs["y"]
        self.z: float = kwargs["z"]

    def __repr__(self):
        return f"<{self.id:3d} ({self.atom_id:<3}) {self.x:7.3f} {self.y:7.3f} {self.z:7.3f}>"

    @property
    def fmoxyz(self):
        return f"{self.id:>7d} {self.type_symbol:>6}    {self.x:>18.8f}{self.y:>18.8f}{self.z:>18.8f}"