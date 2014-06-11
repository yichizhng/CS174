#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"

#include "uistate.hh"
#include "ArcBall.hh"

#include "main_tet.hh"
#include "geodesic_tet.hh"

// Handles keypresses.  Translates them into mode changes, etc.
void keyboard(unsigned char key, int x, int y);
// Handles reshaping of the program window.
void reshape(const int width, const int height);
// Handles motion of the mouse when a button is being held.
void motion(const int x, const int y);
// Handles mouse clicks and releases.
void mouse(int button, int state, int x, int y);
// Initializes the UI.
void initUI();

int xRes = 800;
int yRes = 800;

UIState *ui;

int ntets, nedges, nverts;
double *pos;
Tet **tets;

static double *dists;

using namespace std;

void redraw() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ui->ApplyViewingTransformation();

  /*
  glDisable( GL_TEXTURE_1D );
  glBegin(GL_POINTS);
  for (int i = 0; i < nverts; ++i) {
    glColor3f(1.0f - 7*dists[i], 1.0f - 5*dists[i], 1.0f - 2*dists[i]);
    glVertex3f(pos[3*i], pos[3*i+1], pos[3*i+2]);
  }
  glEnd();

  glEnable( GL_TEXTURE_1D ); */
  
  
  glBegin(GL_TRIANGLES);
  double m = 10;
  for (int i = 0; i < ntets; ++i) {
    glTexCoord1f(m*dists[tets[i]->verts[0]]);
    glVertex3f(pos[3*tets[i]->verts[0]],
	       pos[3*tets[i]->verts[0]+1],
	       pos[3*tets[i]->verts[0]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[1]]);
    glVertex3f(pos[3*tets[i]->verts[1]],
	       pos[3*tets[i]->verts[1]+1],
	       pos[3*tets[i]->verts[1]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[2]]);
    glVertex3f(pos[3*tets[i]->verts[2]],
	       pos[3*tets[i]->verts[2]+1],
	       pos[3*tets[i]->verts[2]+2]);

    glTexCoord1f(m*dists[tets[i]->verts[0]]);
    glVertex3f(pos[3*tets[i]->verts[0]],
	       pos[3*tets[i]->verts[0]+1],
	       pos[3*tets[i]->verts[0]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[2]]);
    glVertex3f(pos[3*tets[i]->verts[2]],
	       pos[3*tets[i]->verts[2]+1],
	       pos[3*tets[i]->verts[2]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[3]]);
    glVertex3f(pos[3*tets[i]->verts[3]],
	       pos[3*tets[i]->verts[3]+1],
	       pos[3*tets[i]->verts[3]+2]);

    glTexCoord1f(m*dists[tets[i]->verts[0]]);
    glVertex3f(pos[3*tets[i]->verts[0]],
	       pos[3*tets[i]->verts[0]+1],
	       pos[3*tets[i]->verts[0]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[3]]);
    glVertex3f(pos[3*tets[i]->verts[3]],
	       pos[3*tets[i]->verts[3]+1],
	       pos[3*tets[i]->verts[3]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[1]]);
    glVertex3f(pos[3*tets[i]->verts[1]],
	       pos[3*tets[i]->verts[1]+1],
	       pos[3*tets[i]->verts[1]+2]);

    glTexCoord1f(m*dists[tets[i]->verts[1]]);
    glVertex3f(pos[3*tets[i]->verts[1]],
	       pos[3*tets[i]->verts[1]+1],
	       pos[3*tets[i]->verts[1]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[2]]);
    glVertex3f(pos[3*tets[i]->verts[2]],
	       pos[3*tets[i]->verts[2]+1],
	       pos[3*tets[i]->verts[2]+2]);
    glTexCoord1f(m*dists[tets[i]->verts[3]]);
    glVertex3f(pos[3*tets[i]->verts[3]],
	       pos[3*tets[i]->verts[3]+1],
	       pos[3*tets[i]->verts[3]+2]);
  }
  glEnd();
  glutSwapBuffers();
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(xRes, yRes); 
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " meshfile" << endl;
    return -1;
  }
  
  // Open the file as a stream
  ifstream fin(argv[1]);
  if(!fin.good()) {
    cerr << "mesh couldn't be opened" << endl;
    return -1;
  }

  // Mesh format used here is a bastardization of the PLY format, mainly
  // so that I don't have to rewrite my parser much

  string buf;
  do {
    fin >> buf;
    if (buf == "vertex")
      fin >> nverts;
    if (buf == "tet")
      fin >> ntets;
  } while (buf != "end_header");

  // Read in vertex positions
  pos = (double *) malloc(3 * nverts * sizeof(double));
  for (int i = 0; i < nverts; ++i) {
    fin >> pos[3*i] >> pos[3*i+1] >> pos[3*i+2];
    fin.ignore(numeric_limits<streamsize>::max(),'\n');
  }

  int v1, v2, v3, v4;
  tets = (Tet **) malloc(ntets * sizeof(Tet *));
  // Read in tets
  for (int i = 0; i < ntets; ++i) {
    fin >> buf; // Discard the first field
    fin >> v1 >> v2 >> v3 >> v4;
    tets[i] = new Tet(v1, v2, v3, v4);
  }
  cerr << "Done reading" << endl;

  dists = (double *) malloc(nverts * sizeof(double));

  init();
  cerr << "Done initializing" << endl;
  step0();
  cerr << "Done factoring" << endl;
  vector<int> blah;
  blah.push_back(0);
  setverts(blah);
  step1(dists);

  cerr << dists[0] << endl;
  cerr << dists[nverts-1] << endl;

  // Setting up a texture to use

  int stripeImageWidth = 32;
  GLubyte   stripeImage[3*stripeImageWidth];
  stripeImage[0] = 255;
  stripeImage[1] = 255;
  stripeImage[2] = 255;
  for (int j = 1; j < stripeImageWidth; j++) {
    stripeImage[3*j] = 0;
    stripeImage[3*j+1] = 155;
    stripeImage[3*j+2] = 0;
  }
  
  GLuint texID;
  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_1D, texID);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage1D(GL_TEXTURE_1D, 0, 3, stripeImageWidth, 0, GL_RGB, GL_UNSIGNED_BYTE, stripeImage);

  // We want the texture to wrap, so that values outside the range [0, 1] 
  // are mapped into a gradient sawtooth
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glDisable( GL_TEXTURE_GEN_S ); 
  glDisable(GL_TEXTURE_2D);
  glEnable( GL_TEXTURE_1D );
  glBindTexture(GL_TEXTURE_1D, texID);

  // Print out all distances; I'll display stuff later (displaying is relatively
  // easy by comparison -.-;;)
  /*
  for (int i = 0; i < nverts; ++i) {
    cout << dists[i] << endl;
    } */

  glDisable(GL_CULL_FACE);
  //glCullFace(GL_BACK);
  glDepthMask(GL_TRUE);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  
    // Set up projection and modelview matrices
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  // Disable lighting; we really aren't going to use it

  glDisable(GL_LIGHTING);
  initUI();

  //  fprintf(stderr, "%d\n", glGetError());

  glutIdleFunc(glutPostRedisplay);
  glutKeyboardFunc(keyboard);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  glutDisplayFunc(redraw);
  glutMainLoop();
  return 0;
}



void keyboard(unsigned char key, int x, int y)
{
  vector<int> blah;
    switch (key) {
    case 27: // ESC
    case 'q':
    case 'Q':
      exit(0);
      break;
    case 'r':
    case 'R':
      // Reset the source to some other point
      blah.push_back(rand()%nverts);
      setverts(blah);
      step1(dists);
      break;
    }
}

//--------------------------------------------------------------------------
// Handles reshaping of the program window
//--------------------------------------------------------------------------
void reshape(const int width, const int height)
{
    xRes = width;
    yRes = height;
    
    if( width <= 0 || height <= 0 ) return;
    
    ui->WindowX() = width;
    ui->WindowY() = height;
    
    ui->Aspect() = double( width ) / height;
    ui->SetupViewport();
    ui->SetupViewingFrustum();
}

//--------------------------------------------------------------------------
// Handles motion of the mouse when a button is being held
//--------------------------------------------------------------------------
void motion(const int x, const int y)
{
    // Just pass it on to the ui controller.
    ui->MotionFunction(x, y);
}

//--------------------------------------------------------------------------
// Handles mouse clicks and releases
//--------------------------------------------------------------------------
void mouse(const int button, const int state, const int x, const int y)
{
    // Just pass it on to the ui controller.
    ui->MouseFunction(button, state, x, y);
}

//--------------------------------------------------------------------------
// Initializes the UI
//--------------------------------------------------------------------------
void initUI()
{
    ui = new UIState;
    ui->Trans() = Vector3(0, 0, 0);
    ui->Radius() = 2;
    ui->Near() = 0.01;
    ui->Far() = 10;
    ui->CTrans() = Vector3(0,0,-0.5);
    reshape(xRes, yRes);
}
