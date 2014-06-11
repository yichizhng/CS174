// The actual calculationy bits of the geodesic distance finder

#include "main.hh"
#include "geodesic.hh"
#include <cmath>
#include <cstring>
#include <cfloat>
#include <ctime>

#include <iostream>
#include <vector>

#define SMALL_FLOAT 1000000*DBL_EPSILON

using namespace std;

#include "suitesparse/cholmod.h"

double m = 100;

extern int nvertices, ntriangles;
extern double *pos;
extern Triangle **tris;
extern vector<int> edgeset;
extern vector<bool> isbound;

static cholmod_common *workspace;
// Workspace for all cholmod operations; note that this is never cleaned up

double *A;

double t; // Timestep

cholmod_sparse *mat1, *mat2, *mat3;
// mat1 is A-tL_c, mat2 is -L_c
cholmod_factor *L1, *L2, *L3;
// Factoriztions for mat1 and mat2 respectively
cholmod_dense *vec1, *vec2, *vec3;
// vec1 is delta, vec2 is phi, vec3 is b

static inline double len(double x, double y, double z) {
  return sqrt( (x*x) + (y*y) + (z*z) );
}

static inline double cot_asin(double x) {
  // Calculates cotangent of arcsin of x; this has a simple formula, but
  // there is a nasty rounding error to watch out for.
  if (x > 1)
    x = 1;
  return sqrt(1 - x * x) / x;
}

void step0() {
  clock_t tx;
  tx = clock();
  // Step 0 is to calculate relevant triangle data; the cotangents of
  // angles and the areas of triangles and vertices.

  cholmod_triplet *trip1;
  
  workspace = (cholmod_common *)malloc(sizeof(cholmod_common));

  cholmod_start(workspace);
  
  // We build the sparse matrices A - tL_c and -L_c (mat1 and mat2 respectively)
  // by first building the triplet matrices. The number of entries is bounded
  // (the cotan operator has up to three (well, six, but we take advantage of
  // symmetry throughout) entries for each triangle and one for each vertex,
  // and A - tL_c has the same number of entries)

  trip1 = cholmod_allocate_triplet(nvertices, nvertices, nvertices + 
				   3*ntriangles, 1, CHOLMOD_REAL, workspace);
  // stype>0 denotes a symmetric matrix stored in the upper half triangle
  // (i.e. i <= j)
  trip1->nnz = trip1->nzmax;

  // Initialize indices for the first nvertices terms; those are the
  // diagonal entries
  for (int i = 0; i < nvertices; ++i) {
    ((int *)trip1->i)[i] = i; // This is some pretty trippy code
    ((int *)trip1->j)[i] = i;
    ((double *)trip1->x)[i] = SMALL_FLOAT; // Pick your favorite (small) number
  }

  A = (double *) calloc(nvertices, sizeof(double));

  for(int i = 0; i < ntriangles; ++i) {
    // The area is simply one half the length of the
    // cross product of any two edges.
    Triangle *tri = tris[i];

    double dx1, dx2, dy1, dy2, dz1, dz2;
    dx1 = pos[3*tri->vertices[1]] - pos[3*tri->vertices[0]];
    dy1 = pos[3*tri->vertices[1]+1] - pos[3*tri->vertices[0]+1];
    dz1 = pos[3*tri->vertices[1]+2] - pos[3*tri->vertices[0]+2];
    dx2 = pos[3*tri->vertices[2]] - pos[3*tri->vertices[0]];
    dy2 = pos[3*tri->vertices[2]+1] - pos[3*tri->vertices[0]+1];
    dz2 = pos[3*tri->vertices[2]+2] - pos[3*tri->vertices[0]+2];

    // And the cross-product is simply...
    tri->area = 0.5 * len((dy1 * dz2) - (dz1 * dy2), 
			  (dz1 * dx2) - (dx1 * dz2),
			  (dx1 * dy2) - (dy1 * dx2));

    A[tri->vertices[0]] += tri->area / 3;
    A[tri->vertices[1]] += tri->area / 3;
    A[tri->vertices[2]] += tri->area / 3;

    // There are plenty of ways of calculating the cotangents from here,
    // mostly involving trigonometric identities. Since we're going to need
    // side lengths anyway, let's use those.
    double a = len(dx1-dx2, dy1-dy2, dz1-dz2);
    double b = len(dx2, dy2, dz2);
    double c = len(dx1, dy1, dz1);

    
    t += a+b+c;
    

    // We can get sin of each angle from the area and side lengths.
    tri->cotans[0] = cot_asin(2*tri->area / (b*c));
    tri->cotans[1] = cot_asin(2*tri->area / (a*c));
    tri->cotans[2] = cot_asin(2*tri->area / (a*b));

    // Armed with the cotangents we can set up the two positive (semi)definite
    // matrices that we care about, A - tL_c and -L_c. For A we just add
    // 1/3 the area (there are more complicated ideas like using the
    // area of the dual cell instead, but this is much easier), and for L_c,
    // well, we have the cotangents, so we can calculate that too.
    
    // Update diagonal entries
    ((double *)trip1->x)[tri->vertices[0]] += (tri->cotans[1] + tri->cotans[2])/2;
    ((double *)trip1->x)[tri->vertices[1]] += (tri->cotans[0] + tri->cotans[2])/2;
    ((double *)trip1->x)[tri->vertices[2]] += (tri->cotans[1] + tri->cotans[0])/2;

    // Update non-diagonal entries; there's a bit of annoyingness here
    if(tri->vertices[0] < tri->vertices[1]) {
      ((int *)trip1->i)[nvertices + 3*i] = tri->vertices[0];
      ((int *)trip1->j)[nvertices + 3*i] = tri->vertices[1];
    }
    else {
      ((int *)trip1->i)[nvertices + 3*i] = tri->vertices[1];
      ((int *)trip1->j)[nvertices + 3*i] = tri->vertices[0];
    }
    ((double *)trip1->x)[nvertices + 3*i] = -tri->cotans[2]/2;

    if(tri->vertices[1] < tri->vertices[2]) {
      ((int *)trip1->i)[nvertices + 3*i + 1] = tri->vertices[1];
      ((int *)trip1->j)[nvertices + 3*i + 1] = tri->vertices[2];
    }
    else {
      ((int *)trip1->i)[nvertices + 3*i + 1] = tri->vertices[2];
      ((int *)trip1->j)[nvertices + 3*i + 1] = tri->vertices[1];
    }
    ((double *)trip1->x)[nvertices + 3*i + 1] = -tri->cotans[0]/2;

    if(tri->vertices[2] < tri->vertices[0]) {
      ((int *)trip1->i)[nvertices + 3*i + 2] = tri->vertices[2];
      ((int *)trip1->j)[nvertices + 3*i + 2] = tri->vertices[0];
    }
    else {
      ((int *)trip1->i)[nvertices + 3*i + 2] = tri->vertices[0];
      ((int *)trip1->j)[nvertices + 3*i + 2] = tri->vertices[2];
    }
    ((double *)trip1->x)[nvertices + 3*i + 2] = -tri->cotans[1]/2;
    
  }

  mat2 = cholmod_triplet_to_sparse(trip1, 0, workspace);

  t /= 3*ntriangles;
  cerr << "Square of average edge length is " << t << endl;
  t *= m * t;

  if (L2) {
    cholmod_free_factor(&L2, workspace);
  }
  L2 = cholmod_analyze(mat2, workspace);
  cholmod_factorize(mat2, L2, workspace);
  //  cholmod_print_sparse(mat2, "mat2", workspace);
  //  cholmod_print_factor(L2, "fac2", workspace);
  
  // Modify trip1 to generate mat1 (A - tL_c) instead
  for (int i = 0; i < nvertices; ++i) {
    if (((double *)trip1->x)[i] < 2*SMALL_FLOAT) {
      ((double *)trip1->x)[i] = DBL_MIN;
    }
    ((double *)trip1->x)[i] *= t;
    ((double *)trip1->x)[i] += A[i];
  }
  for (int i = nvertices; i < nvertices+3*ntriangles; ++i) {
    ((double *)trip1->x)[i] *= t;
  }
  
  mat1 = cholmod_triplet_to_sparse(trip1, 0, workspace);
  
  if (L1) {
    cholmod_free_factor(&L1, workspace);
  }
  L1 = cholmod_analyze(mat1, workspace);
  cholmod_factorize(mat1, L1, workspace);
  for (int i = 0; i < nvertices; ++i) {
    if (isbound[i]) {
      // Dirichlet boundary conditions (sneaky, sneaky)
      ((double *)trip1->x)[i] = DBL_MAX;
    }
  }
  mat3 = cholmod_triplet_to_sparse(trip1, 0, workspace);
  if (L3) {
    cholmod_free_factor(&L3, workspace);
  }
  L3 = cholmod_analyze(mat3, workspace);
  cholmod_factorize(mat1, L3, workspace);

  cholmod_free_triplet(&trip1, workspace);
  tx = clock() - tx;
  cerr << "Finished prefactoring in " << (double)tx/CLOCKS_PER_SEC << " seconds";
  cerr << endl;
}

vector<int> vertset;

void setupHeat() {
  // Sets up the initial heat vector (a Kronecker delta) based on vertset
  if (vec1 != NULL)
    cholmod_free_dense(&vec1, workspace);
  vec1 = cholmod_zeros(nvertices, 1, CHOLMOD_REAL, workspace);
  for (vector<int>::iterator it = vertset.begin(); it != vertset.end(); ++it) {
    ((double *)vec1->x)[*it] = 1;
  }
}

void setverts(vector<int> &set) {
  vertset = set;
  setupHeat();
}


void step1(double *dists) {
  clock_t t;
  t = clock();
  // Assumes the heat vector has already been set up; outputs into dists,
  // which is assumed to be an array of size nvertices.
  double norm[3];
  if (vec2) {
    cholmod_free_dense(&vec2, workspace);
  }
  vec2 = cholmod_solve(CHOLMOD_A, L1, vec1, workspace);
  
  vec3 = cholmod_solve(CHOLMOD_A, L3, vec1, workspace);
  double sum1 = 0, sum2 = 0;
  
  for (int i = 0; i < nvertices; ++i) {
    sum1 += ((double *)vec2->x)[i];
    sum2 += ((double *)vec3->x)[i];
  }
  for (int i = 0; i < nvertices; ++i) {
    ((double *)vec2->x)[i] = ((((double *)vec3->x)[i])*sum1 + 
			      (((double *)vec2->x)[i])*sum2) / (sum1+sum2);
  }

  double *heats = ((double *)vec2->x); /*
  for (int i = 0; i < nvertices; ++i) {
    cerr << ((double *)vec2->x)[i] << endl;
  }
  cerr << endl; */
  vec3 = cholmod_zeros(nvertices, 1, CHOLMOD_REAL, workspace);
  double *vdiv = ((double *)vec3->x);
  // Solve for gradient over each triangle
  double trigrad[3];
  for (int i = 0; i < ntriangles; ++i) {
    Triangle *tri = tris[i];
    double dx1, dx2, dy1, dy2, dz1, dz2;
    dx1 = pos[3*tri->vertices[1]] - pos[3*tri->vertices[0]];
    dy1 = pos[3*tri->vertices[1]+1] - pos[3*tri->vertices[0]+1];
    dz1 = pos[3*tri->vertices[1]+2] - pos[3*tri->vertices[0]+2];
    dx2 = pos[3*tri->vertices[2]] - pos[3*tri->vertices[0]];
    dy2 = pos[3*tri->vertices[2]+1] - pos[3*tri->vertices[0]+1];
    dz2 = pos[3*tri->vertices[2]+2] - pos[3*tri->vertices[0]+2];

    norm[0] = (dy1*dz2 - dz1*dy2);
    norm[1] = (dz1*dx2 - dx1*dz2);
    norm[2] = (dx1*dy2 - dy1*dx2);

    // In fact there's not even any need to use a unit normal because
    // we're going to normalize the gradient anyway
    //double n = len(norm[0],norm[1],norm[2]);
    //norm[0] /= n;
    //norm[1] /= n;
    //norm[2] /= n;

    // The gradient is the sum of u_i times norm cross opposing edge
    trigrad[0] = 
      heats[tri->vertices[0]] * (norm[1]*(dz2-dz1) - norm[2]*(dy2-dy1)) +
      heats[tri->vertices[1]] * (norm[1]*(-dz2) - norm[2]*(-dy2)) +
      heats[tri->vertices[2]] * (norm[1]*(dz1) - norm[2]*(dy1));
    trigrad[1] = 
      heats[tri->vertices[0]] * (norm[2]*(dx2-dx1) - norm[0]*(dz2-dz1)) +
      heats[tri->vertices[1]] * (norm[2]*(-dx2) - norm[0]*(-dz2)) +
      heats[tri->vertices[2]] * (norm[2]*(dx1) - norm[0]*(dz1));
    trigrad[2] =
      heats[tri->vertices[0]] * (norm[0]*(dy2-dy1) - norm[1]*(dx2-dx1)) +
      heats[tri->vertices[1]] * (norm[0]*(-dy2) - norm[1]*(-dx2)) +
      heats[tri->vertices[2]] * (norm[0]*(dy1) - norm[1]*(dx1));

    double n = len(trigrad[0],trigrad[1],trigrad[2]);
    if (!n) continue; // The triangle has no gradient
    trigrad[0] /= n;
    trigrad[1] /= n;
    trigrad[2] /= n;

    // Update vertex divergences; there's some fancy linear algebra going on

    vdiv[tri->vertices[0]] += 
      (tri->cotans[1] * (dx2 * trigrad[0] + dy2 * trigrad[1] +
			 dz2 * trigrad[2]) +
       tri->cotans[2] * (dx1 * trigrad[0] + dy1 * trigrad[1] +
			 dz1 * trigrad[2]))/2;
    vdiv[tri->vertices[1]] += 
      (tri->cotans[2] * (-dx1 * trigrad[0] + (-dy1) * trigrad[1] +
			 (-dz1) * trigrad[2]) +
       tri->cotans[0] * ((dx2-dx1) * trigrad[0] + (dy2-dy1) * trigrad[1] +
			 (dz2-dz1) * trigrad[2]))/2;
    vdiv[tri->vertices[2]] += 
      (tri->cotans[0] * ((dx1-dx2) * trigrad[0] + (dy1-dy2) * trigrad[1] +
			 (dz1-dz2) * trigrad[2]) +
       tri->cotans[1] * ((-dx2) * trigrad[0] + (-dy2) * trigrad[1] +
			 (-dz2) * trigrad[2]))/2;
  }
  // Since that's all dealt with, solve for the last thing
  /* 
  for (int i = 0; i < nvertices; ++i) {
    cerr << vdiv[i] << endl;
  }
  cerr << endl;
  */ 
  cholmod_free_dense(&vec2, workspace);
  vec2 = cholmod_solve(CHOLMOD_A, L2, vec3, workspace);
  memcpy(dists, vec2->x, nvertices * sizeof(double));

  // Solution to poisson problem is only up to a constant; fortunately
  // this constant is not all that hard to establish
  double min = dists[0];
  for (int i = 1; i < nvertices; ++i) {
    if (dists[i] < min)
      min = dists[i];
  }
  for (int i = 0; i < nvertices; ++i) {
    //cerr << dists[i] << endl;
    dists[i] -= min;
    //cerr << dists[i] << endl;
  }

  cholmod_free_dense(&vec1, workspace);
  cholmod_free_dense(&vec2, workspace);
  cholmod_free_dense(&vec3, workspace);
  t = clock() - t;
  cerr << "Finished solving in " << (double)t/CLOCKS_PER_SEC << " seconds";
  cerr << endl;
}

/* Forward Euler simulation; this is why Step 0 exists separately from step
   1; it calculates some things that are needed for heat simulation */
/*
double dt = 0.00000005; // Time step used for simulation, basically found
// experimentally (I wanted something that flowed reasonably quickly
// and didn't become unstable too quickly)

// Runs Forward-Euler heat flow simulation
void timestep(double *u) {
  // We can calculate the standard discrete Laplacian (see section 3.2.1)
  // via iterating over all triangles in the mesh.

  // We actually do need to zero both of these arrays, since we're going
  // to update them by adding
  double *Lu = (double *) calloc(ntriangles, sizeof(double));

  for (int i = 0; i < ntriangles; ++i) {
    Triangle *tri = tris[i];
    
    // Add the appropriate cotangent terms; we need to do each edge twice,
    // since we're going to need to look at each vertex on the edge.

    // cotans[0] is the cotangent of the angle at vertices[0]
    // and so on. There is a pleasing symmetry between these operations.

    Lu[tri->vertices[0]] += 
      tri->cotans[1] * (u[tri->vertices[2]] - u[tri->vertices[0]])
      + tri->cotans[2] * (u[tri->vertices[1]] - u[tri->vertices[0]]);
    Lu[tri->vertices[1]] += 
      tri->cotans[0] * (u[tri->vertices[2]] - u[tri->vertices[1]])
      + tri->cotans[2] * (u[tri->vertices[0]] - u[tri->vertices[1]]);
    Lu[tri->vertices[2]] += 
      tri->cotans[0] * (u[tri->vertices[1]] - u[tri->vertices[2]])
      + tri->cotans[1] * (u[tri->vertices[0]] - u[tri->vertices[2]]);
  }
  // Now we iterate through vertices to update temperatures

  for (int i = 0; i < nvertices; ++i) {
    if (A[i] == 0)
      continue; 
    // Don't want to deal with this; these are
    // disconnected vertices so we can ignore heat flow to them
    u[i] += dt * Lu[i] / (2*A[i]);

    // Fun statistics for debugging:
    // the max cotangent is about 116. the minimum vertex area is about
    // 1e-08. This combines to give me a headache. My measurements say that
    // I would have to run it with a timestep on the order of magnitude of
    // 1e-8, but then the simulation would never go anywhere. Instead just
    // clamp heat values between 0 and 1. There's some visible weirdness
    // which I just refuse to deal with.
    if (u[i] < 0) {
      u[i] = 0;
    }
    if (u[i] > 1) {
      u[i] = 1;
    }
  }
    //  cerr << "Total heat is " << sum << endl;
  free(Lu);
}
*/
