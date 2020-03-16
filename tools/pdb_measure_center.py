#!/usr/bin/python3

from compchem import *
from compchem.pdb import *
import sys

a = load_pdb(sys.argv[1])
frame0 = a[0]

print(frame0.measure_center())
