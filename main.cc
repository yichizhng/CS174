// Parser-ish thing for the Stanford Bunny.
// Requires CHOLMOD; this is available as an Ubuntu package, or if
// you don't have Ubuntu, you can build it from source yourself.

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

#include "main.hh"
#include "geodesic.hh"
#include "marching.hh"

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

// For now we're just working with bun_zipper.ply (the Stanford Bunny)

// Global variables! Programmers love global variables, don't they?
int xRes = 800;
int yRes = 800;

UIState *ui;

int nvertices;
int ntriangles;
double *pos; // Array of length 3 * nvertices
Triangle **tris;

bool running = false;

double *temps;

using namespace std;

void redraw() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ui->ApplyViewingTransformation();
  glEnable(GL_BLEND);
  glBlendFunc(GL_ZERO, GL_SRC_COLOR);
  //glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  // TODO: replace with indexed list, this is pretty much an ideal application
  // for those
  for (int i = 0; i < ntriangles; ++i) {
    //glColor3f(temps[tris[i]->vertices[0]], 0, 1.0 - temps[tris[i]->vertices[0]]);
    glTexCoord1f(128*temps[tris[i]->vertices[0]]);
    glVertex3f(pos[3*tris[i]->vertices[0]],
	       pos[3*tris[i]->vertices[0]+1],
	       pos[3*tris[i]->vertices[0]+2]);
    //glColor3f(temps[tris[i]->vertices[1]], 0, 1.0 - temps[tris[i]->vertices[1]]);
    glTexCoord1f(128*temps[tris[i]->vertices[1]]);
    glVertex3f(pos[3*tris[i]->vertices[1]],
	       pos[3*tris[i]->vertices[1]+1],
	       pos[3*tris[i]->vertices[1]+2]);
    //glColor3f(temps[tris[i]->vertices[2]], 0, 1.0 - temps[tris[i]->vertices[2]]);
    glTexCoord1f(128*temps[tris[i]->vertices[2]]);
    glVertex3f(pos[3*tris[i]->vertices[2]],
	       pos[3*tris[i]->vertices[2]+1],
	       pos[3*tris[i]->vertices[2]+2]);
  }
  glEnd();
  glDisable(GL_BLEND);
  glutSwapBuffers();
  /*
  if(running) {
    // Run several timesteps; these are faster than drawing anyway
    for (int i = 0; i < 10; ++i)
      timestep(temps);
      }*/
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(xRes, yRes); 
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);
  // Open the bunny if no file specified
  ifstream fin((argc==1)?"bun_zipper.ply":argv[1]);
  if (!fin.good())
    return -1;
  
  // Technically speaking there's a lot of things we *should* check for
  // (like the actual number of vertices and triangles)
  // but I'm not going to; instead just ignore everything until end_header

  string buf;
  do {
    fin >> buf;
    if (buf.compare("vertex") == 0)
      fin >> nvertices;
    if (buf.compare("face") == 0)
      fin >> ntriangles;
  } while (buf.compare("end_header"));

  //nvertices = 35947;
  //ntriangles = 69451;

  // Read in vertices

  pos = (double *) malloc(3 * nvertices * sizeof(double));
  for (int i = 0; i < nvertices; ++i) {
    // Read in the coordinates (ignore the other fields)
    fin >> pos[3*i] >> pos[3*i+1] >> pos[3*i+2];
    fin.ignore(numeric_limits<streamsize>::max(),'\n');
  }

  int v1, v2, v3;
  tris = (Triangle **) malloc(ntriangles * sizeof(Triangle *));
  // Read in triangles
  for (int i = 0; i < ntriangles; ++i) {
    fin >> buf; // Discard the first field
    fin >> v1 >> v2 >> v3;
    tris[i] = new Triangle(v1, v2, v3);
  }
  
  temps = (double *) malloc(nvertices* sizeof(double));

  edgesetup();
  findboundary();
  step0();
  vector<int> blah;
  blah.push_back(0);
  setverts(blah);
  step1(temps);

  // Setting up a texture to use

  int stripeImageWidth = 32;
  GLubyte   stripeImage[3*stripeImageWidth];
  stripeImage[0] = 255;
  stripeImage[1] = 255;
  stripeImage[2] = 255;
  for (int j = 1; j < stripeImageWidth; j++) {
    stripeImage[3*j] = 0; // use a gradient instead of a line
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


  // Irrelevant OpenGL Crap

  //glEnable(GL_CULL_FACE);
  //glCullFace(GL_FRONT_AND_BACK);
  //glDepthMask(GL_TRUE);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  glShadeModel(GL_SMOOTH);

  // Disable depth test (for transparency)
  glDisable(GL_DEPTH_TEST);

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
      blah.push_back(rand()%nvertices);
      setverts(blah);
      step1(temps);
      break;
    case 'b':
    case 'B':
      boundset();
      break;
    }
}

// Edges are represented by vertex/vertex adjacency lists
// Of course they only need to be calculated once, but it's kind of a pain
vector<int> *edgeset = NULL;
void edgesetup() {
  if (edgeset) {
    delete[] edgeset;
    edgeset = NULL;
  }
  edgeset = new vector<int>[nvertices];
  // Iterate through the triangles
  for (int i = 0; i < ntriangles; ++i) {
    Triangle *tri = tris[i];
    // Add the six (!) edges
    edgeset[tri->vertices[0]].push_back(tri->vertices[2]);
    edgeset[tri->vertices[0]].push_back(tri->vertices[1]);
    edgeset[tri->vertices[1]].push_back(tri->vertices[0]);
    edgeset[tri->vertices[1]].push_back(tri->vertices[2]);
    edgeset[tri->vertices[2]].push_back(tri->vertices[1]);
    edgeset[tri->vertices[2]].push_back(tri->vertices[0]);
  }
}

vector<bool> isbound;

void findboundary() {
  isbound = vector<bool>(nvertices);
  for (int i = 0; i < nvertices; ++i) {
    // We can identify a boundary vertex as one which only has one
    // edge to some other vertex
    if (edgeset[i].size()%2) {
      // An odd number of edges guarantees a boundary
      isbound[i] = true;
      continue;
    }
    sort(edgeset[i].begin(), edgeset[i].end());
    for (int j = 0; j < edgeset[i].size(); j += 2) {
      if (edgeset[i][j] != edgeset[i][j+1]) {
	isbound[i] = true;
	break;
      }
    }
  }
}

// Solves heat problem with boundary set
void boundset() {
  if (edgeset == NULL)
    edgesetup();
  vector<int> b;
  for (int i = 0; i < nvertices; ++i) {
    // We can identify a boundary vertex as one which only has one
    // edge to some other vertex
    if (edgeset[i].size()%2) {
      // An odd number of edges guarantees a boundary
      b.push_back(i);
      continue;
    }
    sort(edgeset[i].begin(), edgeset[i].end());
    for (int j = 0; j < edgeset[i].size(); j += 2) {
      if (edgeset[i][j] != edgeset[i][j+1]) {
	b.push_back(i);
	break;
      }
    }
  }
  setverts(b);
  step1(temps);
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
