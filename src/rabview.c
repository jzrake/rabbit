#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#define ESCAPE_KEY 27


static void GLUTDisplayFunc();
static void GLUTIdleFunc();
static void GLUTReshapeFunc(int width, int height);
static void GLUTKeyboardFunc(unsigned char key, int x, int y);
static void GLUTSpecialFunc(int key, int x, int y);

static double TranslateZ = -4.0;
static double RotationX = 0.0;
static double RotationY = 0.0;
static double RotationZ = 0.0;

#define RABBIT_INTERNAL
#include "rabbit.h"
static rabbit_mesh *Mesh;


int main(int argc, char **argv)
{

  if (argc == 1) {
    printf("usage: rabiew infile.mesh\n");
    return 0;
  }
  else {
    Mesh = rabbit_mesh_load(argv[1]);
  }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowSize(768, 768);
  glutInitWindowPosition(512, 0);
  glutCreateWindow("rabview");
  glutDisplayFunc(GLUTDisplayFunc);
  glutIdleFunc(GLUTIdleFunc);
  glutReshapeFunc(GLUTReshapeFunc);
  glutKeyboardFunc(GLUTKeyboardFunc);
  glutSpecialFunc(GLUTSpecialFunc);
  glutMainLoop();
  return 0;
}

/* --------------------------------------------------------------------------
 * GLUT function callbacks
 * --------------------------------------------------------------------------
 */
void GLUTDisplayFunc()
{
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  /* reset and configure camera */
  glLoadIdentity();

  glTranslated(0.0, 0.0, TranslateZ);
  glRotated(RotationX, 1.0, 0.0, 0.0);
  glRotated(RotationY, 0.0, 1.0, 0.0);
  glRotated(RotationZ, 0.0, 0.0, 1.0);
  glRotated(-90.0, 1.0, 0.0, 0.0);

  /* actually draw stuff */
  glEnable(GL_DEPTH_TEST);

  if (1) {
    double C = (double) (1 << Mesh->config.max_depth);
    rabbit_edge *edge, *tmp_edge;
    rabbit_geom geom;

    glBegin(GL_LINES);

    HASH_ITER(hh, Mesh->edges, edge, tmp_edge) {

      geom = rabbit_mesh_geom(Mesh, edge->rnp);

      double x0 = geom.vertices[0] / C - 0.5;
      double y0 = geom.vertices[1] / C - 0.5;
      double z0 = geom.vertices[2] / C - 0.5;
      double x1 = geom.vertices[3] / C - 0.5;
      double y1 = geom.vertices[4] / C - 0.5;
      double z1 = geom.vertices[5] / C - 0.5;
      /*
      double r2 = sqrt((x0 + x1) * (x0 + x1) +
		       (y0 + y1) * (y0 + y1) +
		       (z0 + z1) * (z0 + z1)) * 0.5;

      if (r2 < 0.5)
      */
      glVertex3d(x0, y0, z0);
      glVertex3d(x1, y1, z1);
    }
    glEnd();
  }

  /* end drawing */
  glutSwapBuffers();
}

void GLUTIdleFunc()
{
  glutPostRedisplay();
}

void GLUTReshapeFunc(int width, int height)
{
  if (height == 0) {
    height = 1;
  }
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)width/(GLfloat)height, 0.1f, 1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

void GLUTKeyboardFunc(unsigned char key, int x, int y)
{
  switch (key) {
  case ESCAPE_KEY: exit(0);
  case '-': TranslateZ *= 1.1; break;
  case '=': TranslateZ /= 1.1; break;
  }
}

void GLUTSpecialFunc(int key, int x, int y)
{
  double a = 4.0;
  switch (key) {
  case GLUT_KEY_RIGHT: RotationY -= a; break;
  case GLUT_KEY_LEFT:  RotationY += a; break;
  case GLUT_KEY_UP:    RotationX -= a; break;
  case GLUT_KEY_DOWN:  RotationX += a; break;
  }
  glutPostRedisplay();
}
