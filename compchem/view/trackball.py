'''
/*
 * Trackball.h
 * A virtual trackball implementation
 * Written by Gavin Bell for Silicon Graphics, November 1988.

 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 * Vector manip code:
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */
'''

import numpy as np
import math
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

#def build_rotmatrix(float m[4][4], float q[4]);
#def axis_to_quat(float a[3], float phi, float q[4]);

class Trackball:
  size = 400
  quat = np.array([0., 0., 0., 1.])

  def __init__(self):
    pass

  def axis_to_quat(self, a, phi): # axis, angle, quaternion
    a = a / np.linalg.norm(a)
    q = np.zeros(4)
    q[0] = a[0]
    q[1] = a[1]
    q[2] = a[2]
    q *= math.sin(phi / 2.0)
    q[3] = math.cos(phi / 2.0)
    return q

  def project_to_sphere(self, r, x, y):
    d = math.sqrt(x*x + y*y)
    if d < r * 0.70710678118654752440: #inside sphere
      z = math.sqrt(r*r - d*d)
    else:
      t = r / 1.41421356237309504880
      z = t*t / d
    return z

  def calculate_rotation(self, p1x, p1y, p2x, p2y):
    a = np.zeros(3) # axis of rotation
    phi = 0.
    p1 = np.zeros(3)
    p2 = np.zeros(3)
    d = np.zeros(3)
    t = 0.

    if p1x == p2x and p1y == p2y:
      Rm = self.quat_to_matrix()
      return Rm

    # figure out z-coords for projection of p1 and p2 to deformed sphere
    p1 = np.array([p1x, p1y, self.project_to_sphere(self.size, p1x, p1y)])
    p2 = np.array([p2x, p2y, self.project_to_sphere(self.size, p2x, p2y)])
   
    # P2 x P1
    a = np.cross(p2, p1)

    # rotate how much around axis, a?
    d = p1 - p2
    t = np.linalg.norm(d) / (2. * self.size)

    if t >  1.: t =  1.
    if t < -1.: t = -1.
    phi = 2.0 * math.asin(t)

    self.quat = self.axis_to_quat(a, phi)
    #d_quat = self.axis_to_quat(a, phi)
    #self.quat = self.quat + d_quat
    self.quat /= np.linalg.norm(self.quat)

    '''
    q1 = self.quat
    q2 = d_quat
    t1 = q1.copy()
    t1 *= q2[3]
    t2 = q2.copy()
    t2 *= q1[3]
    t3 = np.cross(q2, q1)
    tf = t1 + t2
    tf += t3
    tf[3] = q1[3] * q2[3] - q1.dot(q2)

    self.quat = tf
    '''

    Rm = self.quat_to_matrix()
    return Rm

  def quat_to_matrix(self):
    m = np.zeros((3,3))
    q = self.quat

    m[0,0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
    m[0,1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
    m[0,2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);

    m[1,0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
    m[1,1]= 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
    m[1,2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);

    m[2,0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
    m[2,1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
    m[2,2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);

    return m

  def quat_to_glmatrix(self):
    m = self.quat_to_matrix()
    a = (GLfloat * 16)()
    n = 0
    for i in range(3):
      for j in range(3):
        a[n] = m[i][j]
        n += 1
    return a
