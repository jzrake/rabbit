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

static double TranslateZ = -2.0;
static double RotationX = 90.0;
static double RotationY = 0.0;
static double RotationZ = 0.0;

#define RABBIT_INTERNAL
#include "rabbit.h"

static rabbit_mesh *Mesh;
static double MeshCentroid[3];

static int NumDrawModes = 2;
static int DrawMode = 0;


int main(int argc, char **argv)
{

  if (argc == 1) {
    printf("usage: rabiew infile.mesh\n");
    return 0;
  }
  else {
    Mesh = rabbit_mesh_load(argv[1]);

    rabbit_geom geom;
    rabbit_node *node, *tmp_node;
    double tot_vol = 0.0;
    int num_nodes = 0;
    int D = Mesh->config.max_depth;

    HASH_ITER(hh, Mesh->nodes, node, tmp_node) {

      geom = rabbit_mesh_geom(Mesh, node->rnp);

      double dx = 1.0 / (1 << geom.index[0]);
      double dy = 1.0 / (1 << geom.index[0]);
      double dz = 1.0 / (1 << geom.index[0]);
      double vol = dx * dy * dz;

      MeshCentroid[0] += node->rnp[0] * vol / (1 << D);
      MeshCentroid[1] += node->rnp[1] * vol / (1 << D);
      MeshCentroid[2] += node->rnp[2] * vol / (1 << D);

      tot_vol += vol;
      num_nodes += 1;
    }

    MeshCentroid[0] /= tot_vol;
    MeshCentroid[1] /= tot_vol;
    MeshCentroid[2] /= tot_vol;
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

  if (DrawMode == 00) {
    double C = (double) (1 << Mesh->config.max_depth);
    rabbit_edge *edge, *tmp_edge;
    rabbit_geom geom;

    glBegin(GL_LINES);

    HASH_ITER(hh, Mesh->edges, edge, tmp_edge) {

      geom = rabbit_mesh_geom(Mesh, edge->rnp);

      double x0 = geom.vertices[0] / C - MeshCentroid[0];
      double y0 = geom.vertices[1] / C - MeshCentroid[1];
      double z0 = geom.vertices[2] / C - MeshCentroid[2];
      double x1 = geom.vertices[3] / C - MeshCentroid[0];
      double y1 = geom.vertices[4] / C - MeshCentroid[1];
      double z1 = geom.vertices[5] / C - MeshCentroid[2];

      glVertex3d(x0, y0, z0);
      glVertex3d(x1, y1, z1);
    }
    glEnd();
  }
  if (DrawMode == 1) {
    double C = (double) (1 << Mesh->config.max_depth);
    rabbit_geom geom;
    rabbit_face *face, *tmp_face;
    rabbit_node *node;
    int h, node_rnp[3];

    glBegin(GL_QUADS);

    HASH_ITER(hh, Mesh->faces, face, tmp_face) {

      if (face->rnp[2] != 0) continue;

      geom = rabbit_mesh_geom(Mesh, face->rnp);

      double x0 = geom.vertices[0 ] / C - MeshCentroid[0];
      double y0 = geom.vertices[1 ] / C - MeshCentroid[1];
      double z0 = geom.vertices[2 ] / C - MeshCentroid[2];
      double x1 = geom.vertices[6 ] / C - MeshCentroid[0];
      double y1 = geom.vertices[7 ] / C - MeshCentroid[1];
      double z1 = geom.vertices[8 ] / C - MeshCentroid[2];
      double x2 = geom.vertices[9 ] / C - MeshCentroid[0];
      double y2 = geom.vertices[10] / C - MeshCentroid[1];
      double z2 = geom.vertices[11] / C - MeshCentroid[2];
      double x3 = geom.vertices[3 ] / C - MeshCentroid[0];
      double y3 = geom.vertices[4 ] / C - MeshCentroid[1];
      double z3 = geom.vertices[5 ] / C - MeshCentroid[2];

      h = face->mesh->config.max_depth - geom.index[0] - 1;

      node_rnp[0] = face->rnp[0];
      node_rnp[1] = face->rnp[1];
      node_rnp[2] = face->rnp[2] + 1 << h;

      node = rabbit_mesh_containing(Mesh, node_rnp, RABBIT_RNP);

      double r = node->data[0];
      double g = node->data[0];
      double b = node->data[0];

      glColor3d(r, g, b);
      glVertex3d(x0, y0, z0);
      glVertex3d(x1, y1, z1);
      glVertex3d(x2, y2, z2);
      glVertex3d(x3, y3, z3);
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
  case 'm': DrawMode = (DrawMode + 1) % NumDrawModes;
  }
}

void GLUTSpecialFunc(int key, int x, int y)
{
  double a = 2.0;
  switch (key) {
  case GLUT_KEY_RIGHT: RotationY -= a; break;
  case GLUT_KEY_LEFT:  RotationY += a; break;
  case GLUT_KEY_UP:    RotationX -= a; break;
  case GLUT_KEY_DOWN:  RotationX += a; break;
  }
  glutPostRedisplay();
}
