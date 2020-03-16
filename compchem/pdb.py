
import os, os.path
import sys
import struct #for hex->dec
import re # regex
import numpy as np
import math
from compchem import *

## Function for loading PDB, PQR, PDBQT
def load_pdb(filename):
  if filename.lower().endswith('pdbqt'):
    return PDB(filename, PDBQT=True)
  elif filename.lower().endswith('pqr'):
    return PDB(filename, PQR=True)
  #else:
  return PDB(filename)


def load_pqr(filename):
  return PDB(filename, PQR=True)


def load_pdbqt(filename):
  return PDB(filename, PDBQT=True)



# Single frame from PDB-type file
class PDBFrame(MolecularFrame):
  def __init__(self):
    pass

  def __getitem__(self, key):
    adat = { 'index': self.indices[key],
             'name':  self.names[key],
             'type':  self.types[key],
             'chain':  self.chains[key],
             'resname':  self.resnames[key],
             'resid':  self.resids[key],
             'coordinates':  self.coordinates[key,:],
             'occupancy':  self.occupancies[key],
             'temp_factor':  self.temp_factors[key],
             'charge':  self.charges[key],
             'radius':  self.radii[key], }
    return adat

  def __repr__(self):
    pdb_atoms = ''
    for i in range(len(self.indices)):
      pdb_atoms += 'ATOM  %s  %3s %3s %s %s    %8.3f%8.3f%8.3f  1.00  0.00\n' % (str(self.indices[i]).rjust(5), self.names[i], self.resnames[i], self.chains[i], str(self.resids[i]).rjust(3), self.coordinates[i,0], self.coordinates[i,1], self.coordinates[i,2])
    pdb_atoms += "TER\n"
    return pdb_atoms



## Core PDB-style Format Parser
class PDB:
  def __init__(self, filename, PQR=False, PDBQT=False):
    # check to make sure file exists
    assert (os.path.isfile(filename)), "File does not exist (%s)!" % filename
    # assign object variables
    self.frames = []
    self.current_frame = PDBFrame()
    self.filename = filename
    # data storage
    self.title = ''
    self.authors = ''
    self.journal = {}
    self.journal_string = ''
    # assign type
    self.PQR = PQR
    self.PDBQT = PDBQT
    # parse the file
    self.parse()


  def __repr__(self):
    filetype = 'PDB'
    if self.PQR: filetype = 'PQR'
    if self.PDBQT: filetype = 'PDBQT'
    rep = '* %s-style object (%s).' % (filetype, self.filename)
    rep += '\n  + %d frames' % len(self.frames)
    if self.title != '': rep += '\n  + ' + self.title.replace("\n", "\n    ")
    if self.journal_string != '': rep += '\n  + ' + self.journal_string.replace("\n", "\n    ")
    return rep


  def __getitem__(self, indices):
    return self.frames[indices]


  def parse_index(self, val):
    # When PDB atom indices or residue IDs get over 9999, it switches to hex
    # This function figures that out, and returns an int
    try:
      if re.search('[a-zA-Z]', val):
        return int(val, 16)
      else:
        return int(val)
    except ValueError:
      print('!! compchem.pdb.parse_index could not resolve either hex or int from string "%s"' % val)
      return None
      

  def parse(self):
    # temp vars for empty records
    _index = 1
    _findex = False
    _chain = 'A'
    _fchain = False
    _coords = []

    for i,line in enumerate(open(self.filename)):
      line = line.rstrip()
      if len(line) >= 3:
        ## RECORD
        record = line[:6].strip()
        if record == 'TITLE':
          if self.title != '': self.title += '\n'
          self.title += line[10:]
        elif record == 'AUTHOR':
          if self.title != '': self.title += '\n'
          self.authors += line[10:]
        elif record == 'JRNL':
          key = line[12:17].strip()
          val = line[19:]
          try:
            self.journal[key] += '\n' + val
          except KeyError:
            self.journal[key] = val
        elif record == 'TER':
          _chain = chr(ord(_chain)+1)
        elif record == 'END' or record == 'ENDMDL':
          self.finalize_frame(_coords)
          _coords = []
          self.frames.append(self.current_frame)
          self.current_frame = PDBFrame()
        elif record == 'ATOM' or record == 'HETATM':
          ## ATOM INDEX
          index = line[6:11].strip()
          if index == '':
            _findex = True
            index = _index
            _index += 1
          else:
            index = self.parse_index(index)
          self.current_frame.indices.append(index)
          ## ATOM NAME
          atomn = line[12:16].strip()
          if atomn == '': atomn = 'X'
          self.current_frame.names.append(atomn)
          type_from_name = ''.join([c for c in atomn if not c.isdigit()])
          # RESIDUE NAME
          resna = line[17:21].strip()
          if resna == '': resna = 'X'
          self.current_frame.resnames.append(resna)
          ## CHAIN ID
          chain = line[21:22].strip()
          if chain == '':
            _fchain = True
            chain = _chain
          self.current_frame.chains.append(chain)
          ## RESIDUE NUMBER
          resnu = line[22:26].strip()
          if resnu == '':
            resnu = None
          else:
            resnu = self.parse_index(resnu)
          self.current_frame.resids.append(resnu)
          ## X, Y, Z Coordinates
          crdx = line[30:38]
          assert(crdx != ''), "No coordinate X in PDB (strict parsing) line# %d in %s\n%s" % (i, self.filename, line)
          crdy = line[38:46]
          assert(crdy != ''), "No coordinate Y in PDB (strict parsing) line# %d in %s\n%s" % (i, self.filename, line)
          crdz = line[46:52]
          assert(crdz != ''), "No coordinate Z in PDB (strict parsing) line# %d in %s\n%s" % (i, self.filename, line)
          try:
            crdx = float(crdx)
            crdy = float(crdy)
            crdz = float(crdz)
            extra = line[54:]
            atom_crd = [crdx, crdy, crdz]
          except:
            block = line[30:]
            d1 = block.find('.')
            d2 = block.find('.', d1+1)
            d3 = block.find('.', d2+1)
            s = int((d2 - d1 - 1) / 2)
            crdx = float(block[d1-s:d1+s])
            crdy = float(block[d2-s:d2+s])
            crdz = float(block[d3-s:d3+s])
            extra = block[d3+s:]
            atom_crd = [crdx, crdy, crdz]
          _coords.append(atom_crd)
          # extra columns beyond coordinates
          if self.PQR or self.PDBQT:
            cols = []
            for c in extra.split():
              if c != '':
                cols.append(c)
            if self.PQR:
              self.current_frame.charges.append(float(cols[-2]))
              self.current_frame.radii.append(float(cols[-1]))
              self.current_frame.types.append(type_from_name)
            if self.PDBQT:
              self.current_frame.charges.append(float(cols[-2]))
              self.current_frame.radii.append(0.)
              atomt = cols[-1].strip()
              if atomt == 'A': atomt = 'C'  # autodock-specific types
              if atomt == 'OA': atomt = 'O'
              if atomt == 'NA': atomt = 'N'
              if atomt == 'HD': atomt = 'H'
              self.current_frame.types.append(atomt)
            self.current_frame.occupancies.append(1.)
            self.current_frame.temp_factors.append(0.)
          else:
            try:
              O = float(line[54:60])
            except:
              O = 1.
            self.current_frame.occupancies.append(O)
            try:
              B = float(line[60:66])
            except:
              B = 0.
            self.current_frame.temp_factors.append(B)
            try:
              qs = line[78:80]
              if qs[1] == '+': qs = qs[0]
              if qs[1] == '-':
                qs = '-%s' % qs[0]
              q = float(qs)
            except:
              q = 0.
            self.current_frame.charges.append(q)
            # try to determine generic type for PDB atom
            atype = line[77:80].strip()
            atomt = ''
            if type_from_name == '' and atype == '': atomt = None
            if type_from_name == '' and atype != '': atomt = atype
            if type_from_name != '' and atype == '': atomt = type_from_name
            if type_from_name != '' and atype != '': atomt = atype
            self.current_frame.types.append(atomt)
            # no radii, just append a zero (lookup table for generic vdw radii?)
            self.current_frame.radii.append(0.)
    if self.current_frame.coordinates == None and len(_coords) > 0:
      self.finalize_frame(_coords)
      _coords = []
      self.frames.append(self.current_frame)

    # make some data more accessible
    if len(self.authors) > 0:
      self.authors = [name.strip() for name in self.authors.split(',')]

  def finalize_frame(self, coords):
    Natoms = len(coords)
    if Natoms > 0:
      self.current_frame.coordinates = np.zeros((Natoms, 3))
      for i,atom in enumerate(coords):
        # store x, y, z data
        self.current_frame.coordinates[i,0] = coords[i][0]
        self.current_frame.coordinates[i,1] = coords[i][1]
        self.current_frame.coordinates[i,2] = coords[i][2]
    try: self.journal_string  = self.journal['TITL']
    except KeyError: pass
    try: self.journal_string += '\n' + self.journal['AUTH']
    except KeyError: pass
    try: self.journal_string += '\n' + self.journal['REF']
    except KeyError: pass
    # convert lists to arrays
    self.current_frame.chains =       np.array(self.current_frame.chains)
    self.current_frame.indices =      np.array(self.current_frame.indices)
    self.current_frame.names =        np.array(self.current_frame.names)
    self.current_frame.types =        np.array(self.current_frame.types)
    self.current_frame.resnames =     np.array(self.current_frame.resnames)
    self.current_frame.resids =       np.array(self.current_frame.resids)
    self.current_frame.occupancies =  np.array(self.current_frame.occupancies )
    self.current_frame.temp_factors = np.array(self.current_frame.temp_factors)
    self.current_frame.charges =      np.array(self.current_frame.charges)
    self.current_frame.radii =        np.array(self.current_frame.radii)




#################################
#for fn in os.listdir('../../PDBdb/'):
#  pdb = load_pdb('%s%s' % ("../../PDBdb/",fn))
#  print(pdb)

