#!/usr/bin/python3

from compchem import *
from compchem.pdb import *
from compchem.view import *

a = load_pdbqt('out-methyl-L-ph.pdbqt')
frame0 = a[0]
frame0.rotate(.57, 0., 0.)

e = View(frame0)
e.show()
