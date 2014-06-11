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

int ntets, nedges, nverts;
double *pos;
Tet **tets;

static double *dists;

using namespace std;

int main(int argc, char **argv) {
  glutInit(&argc, argv);
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

  dists = (double *) malloc(nverts * sizeof(double));

  init();
  step0();
  vector<int> blah;
  blah.push_back(0);
  setverts(blah);
  step1(dists);

  // Print out all distances; I'll display stuff later (displaying is relatively
  // easy by comparison -.-;;)

  for (int i = 0; i < nverts; ++i) {
    cout << dists[i] << endl;
  }

  return 0;
}
