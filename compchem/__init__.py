
import os, os.path
import sys
import struct #for hex->dec
import re # regex
import numpy as np
import math




# calculations & manipualtion functions 
def rmsd(frame1, frame2):
  diff = frame1 - frame2
  if len(frame1.shape) == 1 and len(frame2.shape) == 1:
    rms = np.sqrt(np.mean(np.square(diff[0]) + np.square(diff[1]) + np.square(diff[2])))
  else:
    rms = np.sqrt(np.mean(np.square(diff[:,0]) + np.square(diff[:,1]) + np.square(diff[:,2])))
  return round(rms, 3)


# Single frame representing a single molecular system's 
#   conformation and properties
class MolecularFrame:
  chains = []
  indices = []
  names = []
  types = []
  resnames = []
  resids = []
  occupancies = []
  temp_factors = []
  charges = []
  radii = []
  coordinates = None
  bonds = []

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

  def __len__(self):
    return len(self.names)

  def measure_center(self):
    # find mean of each column of coordinates
    x = np.mean(self.coordinates[:,0])
    y = np.mean(self.coordinates[:,1])
    z = np.mean(self.coordinates[:,2])
    mean = np.array( (x,y,z) )
    return mean

  def measure_dimensions(self):
    # find mean of each column of coordinates
    maxx = np.max(self.coordinates[:,0])
    maxy = np.max(self.coordinates[:,1])
    maxz = np.max(self.coordinates[:,2])
    minx = np.min(self.coordinates[:,0])
    miny = np.min(self.coordinates[:,1])
    minz = np.min(self.coordinates[:,2])
    return np.array( (maxx-minx, maxy-miny, maxz-minz) )

  def translate(self, dx, dy, dz):
    # add differential to every row of each column
    self.coordinates[:,0] += dx
    self.coordinates[:,1] += dy
    self.coordinates[:,2] += dz

  def center(self):
    c = self.measure_center()
    self.translate(-c[0], -c[1], -c[2])

  def rotate(self, d0x, d0y, d0z):
    ## rotate around center of coordinates, so we need to center
    ## the molecule on the origin before rotation, then translate back
    c = self.measure_center()
    self.translate(-c[0], -c[1], -c[2])
    ## NOTE: d0's must be radians; add check TODO
    # rotate around x-axis
    cos0 = math.cos(d0x)
    sin0 = math.sin(d0x)
    Rx = np.matrix([[1.,0.,0.],[0., cos0, -sin0],[0., sin0, cos0]])
    # rotate around y-axis
    cos0 = math.cos(d0y)
    sin0 = math.sin(d0y)
    Ry = np.matrix([[cos0, 0., sin0],[0.,1.,0.],[-sin0, 0., cos0]])
    # rotate around z-axis
    cos0 = math.cos(d0z)
    sin0 = math.sin(d0z)
    Rz = np.matrix([[cos0, -sin0, 0.],[sin0, cos0, 0.],[0.,0.,1.]])
    # combine to a single rotation matrix
    Rt = Rx*Ry*Rz
    # apply rotation matrix to coordinates
    self.coordinates = np.array(self.coordinates * Rt)
    # translate molecule back to its original center
    self.translate(c[0], c[1], c[2])

  def rotate_matrix(self, Rt, at_origin=False): #rotation matrix
    # if needed, center on origin
    c = [0., 0., 0.]
    if not at_origin:
      c = self.measure_center()
      self.center()
    # rotate
    self.coordinates = self.coordinates.dot(Rt)
    # translate back to position if needed
    if not at_origin:
      self.translate(c[0], c[1], c[2])

  def rmsd(self, other_frame):
    return rmsd(self.coordinates, other_frame.coordinates)

  def align(self, other_frame):
    X = self.coordinates
    Y = other_frame.coordinates

    ## calculate covariance matrix
    C = X.T.dot(Y)

    #u,s,vh = np.linalg.svd(C)

    #d = int(np.sign(np.linalg.det(C)))
    #T = np.matrix('1 0 0; 0 1 0; 0 0 %d' % d)

    #U = vh.T.dot(T).dot(u.T)
    #mobile.rotate_matrix(U, at_origin=True)
    #mobile.translate(c[0], c[1], c[2])
    #target.translate(c[0], c[1], c[2])
