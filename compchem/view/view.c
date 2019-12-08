#include "numb_visual.h"
#include "../c/model.h"
#include "../c/atom.h"
#include "../c/bond.h"


// Garbage collection
int view_traverse(view *self, visitproc visit, void *arg) {
  py_traverse(self->mdl);
  return 0;
}


int view_clear(view *self) {
  py_clear(self->mdl);
  return 0;
}


// Destructor
void view_dealloc(view* self) {
  if(SDL_JoystickOpened(0) && self->joystick) SDL_JoystickClose(self->joystick);
  SDL_Quit();
  view_clear(self);
  Py_TYPE(self)->tp_free((PyObject*)self);
}


// Constructor
PyObject* view_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  view *self = (view*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->mdl = (PyObject*)Py_None;
    self->exit = 0;
    self->tick_next = 0;
  }
  return (PyObject*)self;
}


// Initialization
int view_init(view *self, PyObject *args, PyObject *kwds) {
  PyObject *mdl = NULL;
  static char *kwlist[] = { "model", NULL };

  if(! PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &mdl)) return -1; 

  py_init(mdl);

  self->quaternion[0] = 0.;
  self->quaternion[1] = 0.;
  self->quaternion[2] = 0.;
  self->quaternion[3] = 1.;
  self->scale = 1.;

  self->rotating = 0;
  self->zooming = 0;

  // find center of molecular system                                                                                                                                                     
  self->origin[0] = self->origin[1] = self->origin[2] = 0.;
  for(int i=0; i < PyObject_Length(((model*)self->mdl)->atoms); i++) {
    self->origin[0] += numb_array_double_get(((model*)self->mdl)->coordinates, i, 0);
    self->origin[1] += numb_array_double_get(((model*)self->mdl)->coordinates, i, 1);
    self->origin[2] += numb_array_double_get(((model*)self->mdl)->coordinates, i, 2);
  }
  self->origin[0] /= PyObject_Length(((model*)self->mdl)->atoms);
  self->origin[1] /= PyObject_Length(((model*)self->mdl)->atoms);
  self->origin[2] /= PyObject_Length(((model*)self->mdl)->atoms);


  // initialize sdl
  if(SDL_Init(SDL_INIT_VIDEO|SDL_INIT_JOYSTICK) < 0) printf("well shit on sdl!\n");
  if (SDL_NumJoysticks()) {
    self->joystick = SDL_JoystickOpen(0);
    if(self->joystick) {
      SDL_JoystickEventState(SDL_ENABLE);
      printf("- Opened joystick\n");
      printf("  + Name: %s\n", SDL_JoystickName(0));
      printf("  + Number of axes: %d\n", SDL_JoystickNumAxes(self->joystick));
      printf("  + Number of buttons: %d\n", SDL_JoystickNumButtons(self->joystick));
      printf("  + Number of balls: %d\n", SDL_JoystickNumBalls(self->joystick));
    }
  }
  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
  const SDL_VideoInfo* screen_info = SDL_GetVideoInfo();
  self->screen = SDL_SetVideoMode(1024., 768., 32, SDL_OPENGL);
  SDL_WM_SetCaption("numb.visual.view", "numb.visual.view");
  SDL_ShowCursor(1);

  // initialize scene
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat ambientLight[] = { .2, .2, .2, 1. };
  GLfloat diffuseLight[] = { .2, .2, .2, 1. };
  GLfloat specularLight[] = { .2, .2, .2, 1. };
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60., self->screen->w / self->screen->h, .1, 10000.);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  return 0;
}


PyObject* view_update(view *self, PyObject *args, PyObject *kwds) {
  SDL_Event event;
  float d_quat[4];

  // handle events
  while(SDL_PollEvent(&event)) {
    switch(event.type) {
      case SDL_MOUSEMOTION:
        if(self->rotating) {
          trackball(d_quat, (2.0 * self->mouse_x - self->screen->w) / self->screen->w, (self->screen->h - 2.0 * self->mouse_y) / self->screen->h, (2.0 * event.button.x - self->screen->w) / self->screen->w, (self->screen->h - 2.0 * event.button.y) / self->screen->h);
          add_quats(d_quat, self->quaternion, self->quaternion);
        }
        if(self->zooming) {
          self->scale = self->scale * (1.0 + (self->mouse_y - event.button.y) / self->screen->h);
          if(self->scale < .001) self->scale = .001;
          if(self->scale > 10.) self->scale = 10.;
        }
        self->mouse_x = event.button.x;
        self->mouse_y = event.button.y;
        break;
      case SDL_MOUSEBUTTONDOWN:
        if(event.button.button == 1) self->rotating = 1;
        if(event.button.button == 3) self->zooming = 1;
        break;
      case SDL_MOUSEBUTTONUP:
        if(event.button.button == 1) self->rotating = 0;
        if(event.button.button == 3) self->zooming = 0;
        break;
      case SDL_JOYAXISMOTION:
        break;
      case SDL_KEYDOWN:
        if(event.key.keysym.sym == SDLK_ESCAPE)
          self->exit = 1;
        //else
        //  self->key_down(event.key.keysym.sym);
        break;
      case SDL_KEYUP:
        //self->key_up(event.key.keysym.sym);
        break;
      case SDL_QUIT:
        self->exit = 1;
        break;
      default:
        break;
    }
  }

  view_redraw(self, NULL, NULL);

  // tick?
  //if(self->tick_next > SDL_GetTicks()) SDL_Delay(self->tick_next - SDL_GetTicks());
  //self->tick_next = SDL_GetTicks() + ((1./30.) * 1000.);

  Py_INCREF(Py_None);
  return Py_None;
}


PyObject* view_redraw(view *self, PyObject *args, PyObject *kwds) {
  // draw
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(0., 0., 0., 1.);

  glPushMatrix();
    // Take a step back
    glTranslatef(0., 0., -70.);
    // Scale scene
    glScalef(self->scale, self->scale, self->scale);
    // Adjust the light position
    GLfloat positionLight[] = { 0., 0., 100., 1. };
    glLightfv(GL_LIGHT0, GL_POSITION, positionLight);
    // Rotate the scene
    float m[4][4];
    build_rotmatrix(m, self->quaternion);
    glMultMatrixf(&m[0][0]);
    // Adjust translation
    glTranslatef(-self->origin[0], -self->origin[1], -self->origin[2]);
    // Draw molecular system
    for(int i=0; i < PyObject_Length(((model*)self->mdl)->atoms); i++) {
      glPushMatrix();
      glTranslatef(numb_array_double_get(((model*)self->mdl)->coordinates, i, 0), numb_array_double_get(((model*)self->mdl)->coordinates, i, 1), numb_array_double_get(((model*)self->mdl)->coordinates, i, 2));
      atom *a = (atom*)PyList_GetItem(((model*)self->mdl)->atoms, i);
      GLfloat materialColor[] = { 0., 0., 0., 1. };
      switch(a->element) {
        case 1:
          materialColor[0] = materialColor[1] = materialColor[2] = .9;
          break;
        case 6:
          materialColor[0] = .1;
          materialColor[1] = .5;
        materialColor[2] = .5;
          break;
        case 7:
          materialColor[0] = .1;
          materialColor[1] = .3;
          materialColor[2] = .9;
          break;
        case 8:
          materialColor[0] = .9;
          materialColor[1] = .2;
          materialColor[2] = .1;
          break;
        case 16:
          materialColor[0] = 1.;
          materialColor[1] = .83;
          materialColor[2] = .36;
          break;
        default:
          materialColor[0] = materialColor[1] = materialColor[2] = .1;
          break;
      }
      glMaterialfv(GL_FRONT, GL_AMBIENT, materialColor);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, materialColor);
      glMaterialf(GL_FRONT, GL_SHININESS, 10.);
      GLUquadricObj *quadric=gluNewQuadric();
      gluQuadricNormals(quadric, GLU_SMOOTH);
      gluSphere(quadric, a->radius, 30, 30);
      gluDeleteQuadric(quadric);
      glPopMatrix();
    }
    glBegin(GL_LINES);
      for(int i=0; i < PyObject_Length(((model*)self->mdl)->bonds); i++) {
        bond *b = (bond*)PyList_GetItem(((model*)self->mdl)->bonds, i);
        glVertex3f(numb_array_double_get(((model*)self->mdl)->coordinates, b->atom1, 0), numb_array_double_get(((model*)self->mdl)->coordinates, b->atom1, 1), numb_array_double_get(((model*)self->mdl)->coordinates, b->atom1, 2));
        glVertex3f(numb_array_double_get(((model*)self->mdl)->coordinates, b->atom2, 0), numb_array_double_get(((model*)self->mdl)->coordinates, b->atom2, 1), numb_array_double_get(((model*)self->mdl)->coordinates, b->atom2, 2));
      }
    glEnd();
  glPopMatrix();
  SDL_GL_SwapBuffers(); 

  Py_INCREF(Py_None);
  return Py_None;
}
