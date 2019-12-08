#ifndef _view_h_
#define _view_h_
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_NUMB_API
#include <Python.h>
#include <structmember.h>
#include <SDL/SDL.h>                                                            
#include <GL/gl.h>
#include <GL/glu.h>


// view class
typedef struct {
  PyObject_HEAD

  // media interface
  SDL_Surface *screen;
  SDL_Joystick *joystick;
  double tick_next;
  int exit;
  int rotating;
  int zooming;
  float origin[3];
  float quaternion[4];
  float scale;
  float mouse_x;
  float mouse_y;

  // molecular interface
  PyObject *mdl;

} view;


// view properties
static PyMemberDef view_members[] = {
  { "model", T_OBJECT_EX, offsetof(view, mdl), 0, "model" },
  { "tick_next", T_DOUBLE, offsetof(view, tick_next), 0, "tick_next" },
  { "exit", T_INT, offsetof(view, exit), 0, "exit" },
  { NULL }
};


// view methods
extern PyObject* view_update(view *self, PyObject *args, PyObject *kwds);
extern PyObject* view_redraw(view *self, PyObject *args, PyObject *kwds);

static PyMethodDef view_methods[] = {
  { "update", (PyCFunction)view_update, METH_VARARGS|METH_KEYWORDS, "update and redraw the view" },
  { "redraw", (PyCFunction)view_update, METH_VARARGS|METH_KEYWORDS, "redraw the view" },
  { NULL } 
};

extern PyObject*  view_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
extern int        view_init(view *self, PyObject *args, PyObject *kwds);
extern void       view_dealloc(view* self);
extern int        view_clear(view *self);
extern int        view_traverse(view *self, visitproc visit, void *arg);


// view type
static PyTypeObject viewType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "numb.c.view",                       /* tp_name */
  sizeof(view),                      /* tp_basicsize */
  0,                                  /* tp_itemsize */
  (destructor)view_dealloc,          /* tp_dealloc */
  0,                                  /* tp_print */
  0,                                  /* tp_getattr */
  0,                                  /* tp_setattr */
  0,                                  /* tp_reserved */
  0,                                  /* tp_repr */
  0,                                  /* tp_as_number */
  0,                                  /* tp_as_sequence */
  0,                                  /* tp_as_mapping */
  0,                                  /* tp_hash */
  0,                                  /* tp_call */
  0,                                  /* tp_str */
  0,                                  /* tp_getattro */
  0,                                  /* tp_setattro */
  0,                                  /* tp_as_buffer */
  ( Py_TPFLAGS_DEFAULT  |
    Py_TPFLAGS_BASETYPE |
    Py_TPFLAGS_HAVE_GC ),             /* tp_flags */
  "numb view object",                  /* tp_doc */
  (traverseproc)view_traverse,       /* tp_traverse */
  (inquiry)view_clear,               /* tp_clear */
  0,                                  /* tp_richcompare */
  0,                                  /* tp_weaklistoffset */
  0,                                  /* tp_iter */
  0,                                  /* tp_iternext */
  view_methods,                      /* tp_methods */
  view_members,                      /* tp_members */
  0,                                  /* tp_getset */
  0,                                  /* tp_base */
  0,                                  /* tp_dict */
  0,                                  /* tp_descr_get */
  0,                                  /* tp_descr_set */
  0,                                  /* tp_dictoffset */
  (initproc)view_init,               /* tp_init */
  0,                                  /* tp_alloc */
  view_new,                          /* tp_new */
};


#endif
