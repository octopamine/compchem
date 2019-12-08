
import sys
from compchem.view.trackball import *
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *


class Mouse:
  def __init__(self):
    self.x = 0.
    self.y = 0.
    self.buttons = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

class View:
  Rm = np.zeros((3,3))
  mouse = Mouse()
  trackball = Trackball()
  quat = np.array([0., 0., 0., 1.])

  def __init__(self, pdbframe, width=800, height=600):
    self.frame = pdbframe
    self.origin = pdbframe.measure_center()
    self.Rm = self.trackball.calculate_rotation(0, 0, 0, 0)

    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(width, height)
    glutCreateWindow('compchem.view')
    glutDisplayFunc(self.draw)
    glutMouseFunc(self.mouse_buttons)
    glutMotionFunc(self.mouse_motions)
    glutPassiveMotionFunc(None)

    glClearColor(0.,0.,0.,1.)
    glShadeModel(GL_SMOOTH)
    glEnable(GL_CULL_FACE)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    lightZeroPosition = [10.,4.,10.,1.]
    lightZeroColor = [0.8,1.0,0.8,1.0] #green tinged
    glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
    glEnable(GL_LIGHT0)

    glMatrixMode(GL_PROJECTION)
    gluPerspective(40.,1., 1.,1000.)
    glMatrixMode(GL_MODELVIEW)
    gluLookAt(0, 0, 0,
              self.origin[0],self.origin[1],self.origin[2],
              0,1,0)
    glPushMatrix()
    glutMainLoop()

  #def show(self):
  #  self.window = glutCreateWindow("compchem view")
  #  glutDisplayFunc(self.draw)
  #  glutIdleFunc(self.draw)
  #  glutMainLoop()

  def mouse_buttons(self, button, state, x, y):
    self.mouse.x = x
    self.mouse.y = y
    self.mouse.buttons[button] = state

  def mouse_motions(self, x, y):
    # trackball rotation
    if self.mouse.buttons[0] == 0:
      self.Rm = self.trackball.calculate_rotation(self.mouse.x, self.mouse.y, x, y)
      self.frame.rotate_matrix(self.Rm)
      glutPostRedisplay();

    self.mouse.x = x
    self.mouse.y = y

  def draw(self):
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glPushMatrix()
    glScalef(0.1, 0.1, 0.1)

    #glTranslatef(-self.origin[0], -self.origin[1], -self.origin[2])

    for i in range(len(self.frame.names)):
      if self.frame.types[i] == 'O': color = [1.,0.,0.,1.]
      if self.frame.types[i] == 'C': color = [0.,1.,1.,1.]
      if self.frame.types[i] == 'H': color = [1.,1.,1.,1.]
      if self.frame.types[i] == 'N': color = [0.,0.,1.,1.]
      glMaterialfv(GL_FRONT, GL_DIFFUSE, color)
      glPushMatrix()
      crd = self.frame.coordinates[i,:]
      glTranslatef(crd[0], crd[1], crd[2])
      glutSolidSphere(0.8, 20, 20)
      glPopMatrix()

    glPopMatrix()
    glutSwapBuffers()

