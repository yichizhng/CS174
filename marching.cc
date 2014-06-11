// (Crude) implementation of fast marching
#include <queue>
#include <vector>
#include <limits>
#include <cmath>
#include <cfloat>
#include "main.hh"

#include <iostream>

using namespace std;

extern int nvertices, ntriangles;
extern double *pos;
extern Triangle **tris;

extern vector<int> vertset;

class vert_dist {
public:
  int v; // Corresponding vertex
  double dist;
  bool operator<(const vert_dist &rhs) const {
    return dist>rhs.dist;
    // Yes, I know that's reversed; we want the smallest
  }
  vert_dist(int x, double d) : v(x),dist(d) {}
};

vector<vert_dist> *edges;
vector<int> *triadjs;
vector<bool> alive;

static inline double len(double x, double y, double z) {
  return sqrt( (x*x) + (y*y) + (z*z) );
}

void setupmarch() {
  // Sets up the edge vector
  alive = vector<bool>(nvertices, false);
  edges = new vector<vert_dist>[nvertices];
  triadjs = new vector<int>[nvertices];
  for (int i = 0; i < ntriangles; ++i) {
    Triangle *tri = tris[i];

    triadjs[tri->vertices[0]].push_back(i);
    triadjs[tri->vertices[1]].push_back(i);
    triadjs[tri->vertices[2]].push_back(i);

    double dx1, dx2, dy1, dy2, dz1, dz2;
    dx1 = pos[3*tri->vertices[1]] - pos[3*tri->vertices[0]];
    dy1 = pos[3*tri->vertices[1]+1] - pos[3*tri->vertices[0]+1];
    dz1 = pos[3*tri->vertices[1]+2] - pos[3*tri->vertices[0]+2];
    dx2 = pos[3*tri->vertices[2]] - pos[3*tri->vertices[0]];
    dy2 = pos[3*tri->vertices[2]+1] - pos[3*tri->vertices[0]+1];
    dz2 = pos[3*tri->vertices[2]+2] - pos[3*tri->vertices[0]+2];

    // There are plenty of ways of calculating the cotangents from here,
    // mostly involving trigonometric identities. Since we're going to need
    // side lengths anyway, let's use those.
    double a = len(dx1-dx2, dy1-dy2, dz1-dz2);
    double b = len(dx2, dy2, dz2);
    double c = len(dx1, dy1, dz1);

    // Add the six (!) edges
    edges[tri->vertices[0]].push_back(vert_dist(tri->vertices[2],b));
    edges[tri->vertices[0]].push_back(vert_dist(tri->vertices[1],c));
    edges[tri->vertices[1]].push_back(vert_dist(tri->vertices[0],c));
    edges[tri->vertices[1]].push_back(vert_dist(tri->vertices[2],a));
    edges[tri->vertices[2]].push_back(vert_dist(tri->vertices[1],a));
    edges[tri->vertices[2]].push_back(vert_dist(tri->vertices[0],b));
  }
}

priority_queue<vert_dist> q;

void march(double *dists) {
  // Uses the fast marching algorithm to calculate distances
  // outputting into dists
  for (int i = 0; i < nvertices; ++i) {
    dists[i] = numeric_limits<double>::infinity();
  }
  for (int i = 0; i < vertset.size(); ++i) {
    dists[vertset[i]] = 0;
    // And update all adjacent vertices; these are now "close"
    for (int j = 0; j < edges[i].size(); ++j) {
      if (edges[i][j].dist < dists[edges[i][j].v])
	dists[edges[i][j].v] = edges[i][j].dist;
    }
    alive[vertset[i]] = true;
  }
  for (int i = 0; i < nvertices; ++i) {
    if ((dists[i] != 0) && (dists[i] != numeric_limits<double>::infinity())) {
      q.push(vert_dist(i,dists[i]));
      //cerr << i << endl;
    }
  }

  while(!q.empty()) {
    vert_dist x = q.top();
    q.pop();
    alive[x.v] = true;
    cerr << x.v << ' ' <<  x.dist << endl; 
    // Get all adjacent _triangles_ and update through them
    for (int i = 0; i < triadjs[x.v].size(); ++i) {
      Triangle *tri = tris[triadjs[x.v][i]];
      if (alive[tri->vertices[0]] && alive[tri->vertices[1]] && alive[tri->vertices[2]])
	continue; // This triangle is already done
      double dx1, dx2, dy1, dy2, dz1, dz2;
      dx1 = pos[3*tri->vertices[1]] - pos[3*tri->vertices[0]];
      dy1 = pos[3*tri->vertices[1]+1] - pos[3*tri->vertices[0]+1];
      dz1 = pos[3*tri->vertices[1]+2] - pos[3*tri->vertices[0]+2];
      dx2 = pos[3*tri->vertices[2]] - pos[3*tri->vertices[0]];
      dy2 = pos[3*tri->vertices[2]+1] - pos[3*tri->vertices[0]+1];
      dz2 = pos[3*tri->vertices[2]+2] - pos[3*tri->vertices[0]+2];

      // There are plenty of ways of calculating the cotangents from here,
      // mostly involving trigonometric identities. Since we're going to need
      // side lengths anyway, let's use those.
      double a = len(dx1-dx2, dy1-dy2, dz1-dz2);
      double b = len(dx2, dy2, dz2);
      double c = len(dx1, dy1, dz1);

      // There's some casework here which we eliminate by sorting everything
      // by distance
      int cpy[3];
      cpy[0] = tri->vertices[0];
      cpy[1] = tri->vertices[1];
      cpy[2] = tri->vertices[2];
      if (cpy[0] < cpy[1]) {
	if (cpy[1] < cpy[2]) {
	  // Already sorted
	}
	else if (cpy[0] < cpy[2]) {
	  // 0,2,1
	  double tmp = c;
	  c = b;
	  b = tmp;
	  tri->vertices[1] = cpy[2];
	  tri->vertices[2] = cpy[0];
	}
	else {
	  // 2,0,1
	  double tmp = a;
	  a = c;
	  c = b;
	  b = tmp;
	  tri->vertices[0] = cpy[2];
	  tri->vertices[1] = cpy[0];
	  tri->vertices[2] = cpy[1];
	}
      }
      else {
        // 1,0
	if (cpy[0] < cpy[2]) {
	  // 1,0,2
	  double tmp = a;
	  a = b;
	  b = tmp;
	  tri->vertices[1] = cpy[0];
	  tri->vertices[0] = cpy[1];
	}
	// In the other cases our new point is the farthest and there's no
	// point in updating :)
	else {
	  continue;
	}
      }
      
      if ( (!alive[tri->vertices[2]]) &&
dists[tri->vertices[1]] != numeric_limits<double>::infinity()) {
	double costheta = (a*a + b*b - c*c)/(2*a*b);
	if (costheta>1)
	  costheta=1;
	if (costheta<0) // Fudging for right triangles; if it's _actually_
	  // negative you'll get what you deserve
	  costheta = DBL_MIN;
	double sintheta = sqrt(1-costheta*costheta);
	bool topush = (dists[tri->vertices[2]] == numeric_limits<double>::infinity());
	// Using the terminology from Sethian and Kimmel, we now
	// have that A is vertices[0], B is vertices[1], and C
	// is vertices[2]
	
	// Theta is the angle at 2
	  
	double u = dists[tri->vertices[1]] - dists[tri->vertices[0]];
	double t;
	  
	double aa = a*a + b*b - 2*a*b*costheta;
	double bb = 2*b*u*(a*costheta-b);
	double cc = b*b*(u*u-a*a*sintheta*sintheta);
	if (bb*bb >= 4*aa*cc) {
	  t = (sqrt(bb*bb - 4*aa*cc) - bb) / (2*aa);
	    
	  // Check that our conditions for t hold (that the update
	  // occurs from inside the triangle)
	  if ((u<t) && (a*costheta < (b*(t-u)/t)) && ((b*(t-u)/t) < (a/costheta))) {
	    // Update is valid
	    if (t < dists[tri->vertices[2]])
	      dists[tri->vertices[2]] = t;
	  }
	}
	// Try the other two updates too while we're at it
	if (dists[0] + b < dists[tri->vertices[2]]) {
	  dists[tri->vertices[2]] = dists[0] + b;
	}
	if (dists[1] + a < dists[tri->vertices[2]]) {
	  dists[tri->vertices[2]] = dists[1] + a;
	}
	if (topush)
	  q.push(vert_dist(tri->vertices[2],dists[tri->vertices[2]]));
      }
      else {
	// Since we have two infinities what we need to do is obvious
	dists[tri->vertices[1]] = dists[tri->vertices[0]] + c;
	dists[tri->vertices[2]] = dists[tri->vertices[0]] + b;
	q.push(vert_dist(tri->vertices[1], dists[tri->vertices[1]]));
	q.push(vert_dist(tri->vertices[2], dists[tri->vertices[2]]));
      }
      tri->vertices[0] = cpy[0];
      tri->vertices[1] = cpy[1];
      tri->vertices[2] = cpy[2];
    }
  } /*
  for (int i = 0; i < nvertices; ++i) {
    cerr << i << dists[i] << endl;
    } */
}
