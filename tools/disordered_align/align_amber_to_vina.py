#!/usr/bin/python3

import sys
import os, os.path
from compchem import *

# important variables
pdb_fn = ''
pdbqt_fn = ''
pdbqt_frame = 0

# print usage if needed
def usage():
  print(sys.argv[0], '[pdb target filename, filename.pdb]', '[pdbqt mobile filename, pdbqt]', '[pdbqt frame number, e.g. 0]', )

# check arguments
if len(sys.argv) != 4:
  print("Error: Improper arguments passed to script.")
  usage()
  sys.exit()

# grab PDB argument, make sure it exists
pdb_fn = sys.argv[1]
if not os.path.exists(pdb_fn):
  print('Error: PDB file "%s" does not exist. Quitting.' % pdb_fn)
  sys.exit()

# grab PDBQT argument, make sure it exists
pdbqt_fn = sys.argv[2]
if not os.path.exists(pdbqt_fn):
  print('Error: PDBQT file "%s" does not exist. Quitting.' % pdbqt_fn)
  sys.exit()

# grab PDBQT FRAME argument
try:
  pdbqt_frame = int(sys.argv[3])
except ValueError:
  print('Error: PDBQT Frame number "%s" invalid. Can\'t convert to integer.' % sys.argv[3])
  sys.exit()
  

# create our PDB objects, pull out relevant PDBFrame
pdbqt = load_pdbqt(pdbqt_fn)
target = pdbqt[pdbqt_frame]

pdb = load_pdb(pdb_fn)
mobile = pdb[0]

# calculate center of target, center the mobile molecule there
c = target.measure_center()
target.center()
mobile.center()
mcrd_noh = mobile.coordinates[mobile.types != 'H']
tcrd_noh = target.coordinates[target.types != 'H']
mtyp_noh = mobile.types[mobile.types != 'H']
ttyp_noh = target.types[target.types != 'H']


# use network data (bonding information) to find unique "ids" for each atom
def calculate_degree(crd, typ):
  # first calculate the degree of each atom
  # also, collect indices of the connected atoms
  degree_t = []
  connections_t = []
  for i,crdi in enumerate(crd):
    degree = 0
    connections = []
    for j,crdj in enumerate(crd):
      if i == j: continue
      dcrd = crdi - crdj
      dist = np.sqrt(np.square(dcrd[0]) + np.square(dcrd[1]) + np.square(dcrd[2]))
      if dist < 1.8:
        degree += 1
        connections.append(j)
    degree_t.append(degree)
    connections_t.append(connections)
  # get second degree connections
  seconds_t = []
  for i,crdi in enumerate(crd):
    seconds = []
    for c in connections_t[i]:
      seconds.append(connections_t[c])
    seconds_t.append(seconds)
  # get third degree connections
  thirds_t = []
  for i,crdi in enumerate(crd):
    thirds = []
    for sec in seconds_t[i]:
      for c in sec:
        if c != i:
          for tc in connections_t[c]:
            thirds.append(tc)
    thirds_t.append(thirds)
  # print identifier for each atom
  descriptions = []
  for i,crdi in enumerate(crd):
    seconds_flat = [item for sublist in seconds_t[i] for item in sublist]
    thirds_flat = [item for item in thirds_t[i]]
    desc = '%s %d %d %d %d %d %d' % (typ[i], degree_t[i], len(connections_t[i]), len(seconds_flat), len(set(seconds_flat)), len(thirds_flat), len(set(thirds_flat)))
    descriptions.append(desc)

  return descriptions


# determine index translation map
itrans = {}
d1 = calculate_degree(tcrd_noh, ttyp_noh)
d2 = calculate_degree(mcrd_noh, mtyp_noh)

for i in range(len(d1)):
  itj = 0
  cnt = 0
  for j in range(len(d2)):
    if d1[i] == d2[j]:
      itj = j
      cnt += 1
  if cnt == 1:
    itrans[i] = itj

## create matrices of mapped atoms
# create array of bools to subset only the atoms that were mapped
tbool = [False for x in range(len(tcrd_noh))]
for it in itrans.keys():
  tbool[it] = True
tbool = np.array(tbool)

# subset target coordinates to only atoms that got mapped
tcrd_ss = tcrd_noh[tbool]
X = np.matrix(tcrd_ss)

# copy the target coordinate array, and fill in values using atom map
mcrd_noh_ordered = mcrd_noh.copy()
for i in range(len(mcrd_noh)):
  try:
    mcrd_noh_ordered[i] = mcrd_noh[itrans[i]]
  except KeyError:
    mcrd_noh_ordered[i] = np.array([0., 0., 0.])
    pass #we'll hit this for unmapped atoms

# subset target coordinates to only atoms that got mapped
mcrd_ss = mcrd_noh_ordered[tbool]
Y = np.matrix(mcrd_ss)

## finally align based on index translation map
# calculate covariance matrix
C = X.T.dot(Y)

u,s,vh = np.linalg.svd(C)

d = int(np.sign(np.linalg.det(C)))
T = np.matrix('1 0 0; 0 1 0; 0 0 %d' % d)

U = vh.T.dot(T).dot(u.T)
mobile.rotate_matrix(U, at_origin=True)
mobile.translate(c[0], c[1], c[2])
target.translate(c[0], c[1], c[2])

print("REMARK  Initial RMSD:", rmsd(X, Y), 'A')
print("REMARK   Final RMSD:", rmsd(X, Y.dot(U)), 'A')

print(mobile)
