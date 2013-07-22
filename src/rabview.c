#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#define ESCAPE_KEY 27


static void GLUTDisplayFunc();
static void GLUTIdleFunc();
static void GLUTReshapeFunc(int width, int height);
static void GLUTKeyboardFunc(unsigned char key, int x, int y);
static void GLUTSpecialFunc(int key, int x, int y);

static double TranslateZ = -3.0;
static double RotationX = 0.0;
static double RotationY = 0.0;
static double RotationZ = 0.0;


int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowSize(768, 768);
  glutInitWindowPosition(512, 0);
  glutCreateWindow("GLUT Window");
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
  glTranslated(0.0,-1.0, 0.0);
  glTranslated(0.0, 0.0, TranslateZ);
  glRotated(RotationX, 1.0, 0.0, 0.0);
  glRotated(RotationY, 0.0, 1.0, 0.0);
  glRotated(RotationZ, 0.0, 0.0, 1.0);
  glRotated(-90.0, 1.0, 0.0, 0.0);

  /* actually draw stuff */
  glEnable(GL_DEPTH_TEST);
  glutWireTeapot(1.0);

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
  gluPerspective(45.0f, (GLfloat)width/(GLfloat)height, 0.1f, 100.0f);
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
  case GLUT_KEY_RIGHT: RotationZ -= a; break;
  case GLUT_KEY_LEFT:  RotationZ += a; break;
  case GLUT_KEY_UP:    RotationX -= a; break;
  case GLUT_KEY_DOWN:  RotationX += a; break;
  }
  glutPostRedisplay();
}
