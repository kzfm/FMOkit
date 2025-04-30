import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from FMOkit import Atom
a = Atom(x=25.8240, y=21.6710, z=10.2380, id=1, type_symbol="N", atom_id="N", charge=0.194)
print(a.fmoxyz)
