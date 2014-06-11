// The implementation of the geodesic distance finder on a
// 3-simplicial manifold

// The steps are roughly analogous to those of the 2-simplicial manifold:

// The question here is what the Laplacian is; well, it's \delta_1 d_0. Neither
// of these is particularly _difficult_ to calculate, just very tedious.
// (the main annoyance is the question of what *_1 is; we approximate it as
// the diagonal matrix of size E x E with entries equal to the ratio of the area
// of the dual 2-form to the length of the edge. (Of course, we cheat when
// calculating the area of the dual 2-form, since this is an approximation
// anyway)

#include "main_tet.hh"
#include "suitesparse/cholmod.h"
#include <unordered_map>
#include <cmath>
#include <cfloat>
#include <vector>
#include <cstring>

#include <iostream>

// Data structures:

// Constants (well, per mesh)
extern int ntets;
extern int nedges;
extern int nverts;
extern Tet **tets;
extern double *pos;

#define SMALL_FLOAT 1000000*DBL_EPSILON

static cholmod_common *workspace;
cholmod_sparse *mat1, *mat2;

cholmod_factor *L1, *L2;

cholmod_dense *vec1, *vec2, *vec3;

static inline double len(double x, double y, double z) {
  return sqrt( (x*x) + (y*y) + (z*z) );
}

using namespace std;

// An unordered map, allowing us to retrieve an edge index from vertex pairs
unordered_map<Edge, int> vert_edge_map;
// Note that nedges == vert_edge_map.size()

// Returns the edge corresponding to the pair (i,j) should it exist; -1
// otherwise
int map_vert_to_edge(int i, int j) {
  Edge e(i,j);
  if (vert_edge_map.count(e)) {
    return vert_edge_map[e];
  }
  return -1;
}

// Inserts an edge into the vertex pair - edge map
void vem_insert_edge(int i, int j, int e) {
  Edge ee(i,j);
  vert_edge_map[ee] = e;
}

double timestep;

// Sum of incident tet volumes on an edge
double *edge_areas;

double *edge_lens;

// Sum of incident tet volumes on a vertex
double *vert_areas;

// tet volumes
double *tet_volumes;

void init() {
  // Set up things needed for calculations in this file; assumes all extern
  // variables have already been set, except nedges; we also calculate
  // the edge set here as well as the timestep t
  vert_areas = (double *) calloc(nverts, sizeof(double));
  int e = 0;
  for (int i = 0; i < ntets; ++i) {
    Tet *t = tets[i];
    if (map_vert_to_edge(t->verts[0], t->verts[1]) == -1)
      vem_insert_edge(t->verts[0], t->verts[1], e++);
    if (map_vert_to_edge(t->verts[0], t->verts[2]) == -1)
      vem_insert_edge(t->verts[0], t->verts[2], e++);
    if (map_vert_to_edge(t->verts[0], t->verts[3]) == -1)
      vem_insert_edge(t->verts[0], t->verts[3], e++);
    if (map_vert_to_edge(t->verts[1], t->verts[2]) == -1)
      vem_insert_edge(t->verts[1], t->verts[2], e++);
    if (map_vert_to_edge(t->verts[1], t->verts[3]) == -1)
      vem_insert_edge(t->verts[1], t->verts[3], e++);
    if (map_vert_to_edge(t->verts[2], t->verts[3]) == -1)
      vem_insert_edge(t->verts[2], t->verts[3], e++);
  }
  nedges = e;
  
  edge_areas = (double *) calloc(nedges, sizeof(double));
  edge_lens = (double *) calloc(nedges, sizeof(double));
  double edge_total = 0;
  for (auto it = vert_edge_map.begin(); it != vert_edge_map.end(); ++it) {
    double dx, dy, dz;
    dx = pos[3*it->first.v1] - pos[3*it->first.v2];
    dy = pos[3*it->first.v1+1] - pos[3*it->first.v2+1];
    dz = pos[3*it->first.v1+2] - pos[3*it->first.v2+2];
    edge_lens[it->second] = len(dx, dy, dz);
    edge_total += len(dx, dy, dz);
  }
  timestep = edge_total / nedges;
  timestep *= timestep;
  workspace = (cholmod_common *)malloc(sizeof(cholmod_common));
  cholmod_start(workspace);
  mat1 = 0;
  mat2 = 0;
  vec1 = 0;
  vec2 = 0;
  L1 = 0;
  L2 = 0;
}

// Calculate the Laplacian matrix L as well as the diagonal area matrix
// A
void step0() {
  for (int i = 0; i < ntets; ++i) {
    Tet *t = tets[i];
    // Calculate the volume of the tet
    // Call its vertices v0, v1, v2, and v3. Then the volume is simply
    // one sixth the triple product of the edges v1-v0, v2-v0, and v3-v0.
    // Of course I'm trying not to make the calculation orientation
    // dependent, so there will be an absolute value to make sure :)
    double x1, x2, x3;
    double y1, y2, y3;
    double z1, z2, z3;
    x1 = pos[3 * t->verts[1]] - pos[3 * t->verts[0]];
    x2 = pos[3 * t->verts[1]+1] - pos[3 * t->verts[0]+1];
    x3 = pos[3 * t->verts[1]+2] - pos[3 * t->verts[0]+2];
    y1 = pos[3 * t->verts[2]] - pos[3 * t->verts[0]];
    y2 = pos[3 * t->verts[2]+1] - pos[3 * t->verts[0]+1];
    y3 = pos[3 * t->verts[2]+2] - pos[3 * t->verts[0]+2];
    z1 = pos[3 * t->verts[3]] - pos[3 * t->verts[0]];
    z2 = pos[3 * t->verts[3]+1] - pos[3 * t->verts[0]+1];
    z3 = pos[3 * t->verts[3]+2] - pos[3 * t->verts[0]+2];

    t->volume = ( (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2) -
                     (x1 * y3 * z2 + x2 * y1 * z3 + x3 * y2 * z1) ) / 6;
    cerr << t->volume << endl;
    if (t->volume < 0)
      t->volume = -t->volume;

    // We actually need to iterate through the _edges_ of the
    // tet; 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
    // Update the volume of each edge (needed for *_1)
    for (int j = 0; j < 3; ++j)
      for (int k = j+1; k < 4; ++k) {
        int v1 = tets[i]->verts[j];
        int v2 = tets[i]->verts[k];
        int e = map_vert_to_edge(v1, v2);
        edge_areas[e] += t->volume;
      }

    // Update the volume of each vertex (needed for *_3)
    for (int j = 0; j < 4; ++j) {
      vert_areas[j] += t->volume;
    }
  }

  for (int i = 0; i < nedges; ++i) {
    cerr << edge_areas[i] << endl;
  }
  cerr << endl;

  // (i != j)
  // L_(i,j) = (sum over tetrahedra t incident on both i and j)
  // |e*_ij| / |e_ij|
  // |e_ij| is obvious (it's just the length of edge ij) but
  // |e*_ij| isn't; we'll approximate it as edge_areas[ij]
  // Meanwhile L_(i,i) is (sum over j adjacent to i) -L_(i,j)

  cholmod_triplet *trip1;
  // trip1 will be used to build both -L_c and A - tL_c
  
  trip1 = cholmod_allocate_triplet(nverts, nverts, nverts + 6 * ntets, 1,
                                   CHOLMOD_REAL, workspace);
  // 6, because each tet has 6 edges (12, really, but the matrix is symmetric)
  // (12 _directed_ edges, if you're confused)
  // And yes, there is overlap; it's okay, CHOLMOD will deal with it.

  trip1->nnz = trip1->nzmax; // We'll be filling the whole thing

  // Initialize the diagonal entries; they will be "zero" for now, and
  // will be updated every time we add a non-diagonal entry
  for (int i = 0; i < nverts; ++i) {
    ((int *)trip1->i)[i] = i;
    ((int *)trip1->j)[i] = i;
    ((double *)trip1->x)[i] = SMALL_FLOAT;
  }

  for (int i = 0; i < ntets; ++i) {
    Tet *t = tets[i];
    
    // Iterate through edges
    for (int j = 0; j < 3; ++j)
      for (int k = j+1; k < 4; ++k) {
        // I'm kind of faking this expression; draw it out yourself if
        // you'd like to see what's going on
        int v1 = tets[i]->verts[j];
        int v2 = tets[i]->verts[k];
        int e = map_vert_to_edge(v1, v2);
        // (*_1)_e is what we're using as the term
        // This time, find the edge length, because we actually
        // kind of need that
        double dx, dy, dz;
        dx = pos[3*v2] - pos[3*v1];
        dy = pos[3*v2+1] - pos[3*v1+1];
        dz = pos[3*v2+2] - pos[3*v1+2];
        double e_len = len(dx, dy, dz);
        // Update diagonal entries
        ((double *)trip1->x)[v1] += (0.75 * edge_areas[e] / (e_len * e_len));
        ((double *)trip1->x)[v2] += (0.75 * edge_areas[e] / (e_len * e_len));

        // Update non-diagonal entry (we're only doing one per time through
        // this loop, yes)

        int idx = j + k - (!j); 
        // Heh. Don't write code like this.
        if (v1 < v2) {
          ((int *)trip1->i)[nverts + 6 * i + idx] = v1;
          ((int *)trip1->j)[nverts + 6 * i + idx] = v2;
        }
        else {
          ((int *)trip1->i)[nverts + 6 * i + idx] = v2;
          ((int *)trip1->j)[nverts + 6 * i + idx] = v1;
        }
        ((double *)trip1->x)[nverts + 6 * i + idx] = -0.75 * edge_areas[e] /
          (e_len * e_len);
      }
    
  }

  mat2 = cholmod_triplet_to_sparse(trip1, 0, workspace);
  L2 = cholmod_analyze(mat2, workspace);
  cholmod_factorize(mat2, L2, workspace);
  
  cerr << endl;
  // Modify trip1 to generate A - tL
  for (int i = 0; i < nverts; ++i) {
    cerr << ((double *)trip1->x)[i] << endl;
    ((double *)trip1->x)[i] *= timestep;
    ((double *)trip1->x)[i] += vert_areas[i]/4;
  }
  for (int i = nverts; i < nverts + 6 * ntets; ++i) {
    cerr << ((double *)trip1->x)[i] << endl;
    ((double *)trip1->x)[i] *= timestep;
  }

  mat1 = cholmod_triplet_to_sparse(trip1, 0, workspace);

  L1 = cholmod_analyze(mat1, workspace);
  cholmod_print_factor(L1, "fac1", workspace);
  cholmod_print_factor(L2, "fac2", workspace);
}

vector<int> vertset;
void setupHeat() {
  // Sets up initial heat vector based on vertset (a Kronecker delta on
  // vertset, essentially)
  if (vec1 != NULL)
    cholmod_free_dense(&vec1, workspace);
  vec1 = cholmod_zeros(nverts, 1, CHOLMOD_REAL, workspace);
  for (int v : vertset)
    ((double *)vec1->x)[v] = 1;
}

void setverts(const vector<int> &set) {
  vertset = set;
  setupHeat();
}

void step1(double *dists) {
  // Assumes setverts has already been called with the appropriate set.
  // Outputs into dists, which is assumed to be an array of size nverts.

  // Obviously we need to calculate the gradient and divergence here.
  // Of course they are more or less analogous to the 2-D ones (see
  // Keenan Crane's SIGGRAPH 2013 talk) except we don't get the nice
  // cotan form for the Hodge star

  if (vec2) {
    cholmod_free_dense(&vec2, workspace);
  }
  vec2 = cholmod_solve(CHOLMOD_A, L1, vec1, workspace);
  double *heats = ((double *)vec2->x);
  // This gives us vec2 (u)
  vec3 = cholmod_zeros(nverts, 1, CHOLMOD_REAL, workspace);
  double *vdiv = (double *)vec3->x;
  double tetgrad[3]; //gradient of a tet as a 3-vector

  // The gradient over a tet is just 1/3 the sum over vertices of the
  // value of u at that vertex times the area of the opposing face
  // times the normal of the opposing face divided by the volume

  // I leave out some multiplicative factors (like the volume) in this
  // calculation because they are more or less irrelevant
  // (we were going to normalize it anyway)
  for (int i = 0; i < ntets; ++i) {
    Tet *t = tets[i];

    // Zero the tet gradient
    tetgrad[0] = 0;
    tetgrad[1] = 0;
    tetgrad[2] = 0;

    // first vertex
    // Notice that AN is just half the cross-product of two remaing edges
    // (in the appropriate direction; we can check this)
    {
      double dx1, dx2, dx3;
      double dy1, dy2, dy3;
      double dz1, dz2, dz3;
      double nx, ny, nz;
      // These are 3-0, 3-1, and 3-2 respectively
      dx1 = pos[3*t->verts[3]] - pos[3*t->verts[0]];
      dy1 = pos[3*t->verts[3]+1] - pos[3*t->verts[0]+1];
      dz1 = pos[3*t->verts[3]+2] - pos[3*t->verts[0]+2];
      
      dx2 = pos[3*t->verts[3]] - pos[3*t->verts[1]];
      dy2 = pos[3*t->verts[3]+1] - pos[3*t->verts[1]+1];
      dz2 = pos[3*t->verts[3]+2] - pos[3*t->verts[1]+2];
      
      dx3 = pos[3*t->verts[3]] - pos[3*t->verts[2]];
      dy3 = pos[3*t->verts[3]+1] - pos[3*t->verts[2]+1];
      dz3 = pos[3*t->verts[3]+2] - pos[3*t->verts[2]+2];
      
      // We will use the cross product of the latter two, using the first
      // to check orientation
      nx = dy2 * dz3 - dz2 * dy3;
      ny = dz2 * dx3 - dx2 * dz3;
      nz = dx2 * dy3 - dy2 * dx3;
      if (nx * dx1 + ny * dy1 + nz * dz1 < 0) {
        // Negate them
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      // Add this to the tet gradient
      tetgrad[0] += nx;
      tetgrad[1] += ny;
      tetgrad[2] += nz;
    }

    // Second vertex
    {
      double dx1, dx2, dx3;
      double dy1, dy2, dy3;
      double dz1, dz2, dz3;
      double nx, ny, nz;
      // These are 0-1, 0-2, and 0-3 respectively
      dx1 = pos[3*t->verts[0]] - pos[3*t->verts[1]];
      dy1 = pos[3*t->verts[0]+1] - pos[3*t->verts[1]+1];
      dz1 = pos[3*t->verts[0]+2] - pos[3*t->verts[1]+2];
      
      dx2 = pos[3*t->verts[0]] - pos[3*t->verts[2]];
      dy2 = pos[3*t->verts[0]+1] - pos[3*t->verts[2]+1];
      dz2 = pos[3*t->verts[0]+2] - pos[3*t->verts[2]+2];
      
      dx3 = pos[3*t->verts[0]] - pos[3*t->verts[3]];
      dy3 = pos[3*t->verts[0]+1] - pos[3*t->verts[3]+1];
      dz3 = pos[3*t->verts[0]+2] - pos[3*t->verts[3]+2];
      
      // We will use the cross product of the latter two, using the first
      // to check orientation
      nx = dy2 * dz3 - dz2 * dy3;
      ny = dz2 * dx3 - dx2 * dz3;
      nz = dx2 * dy3 - dy2 * dx3;
      if (nx * dx1 + ny * dy1 + nz * dz1 < 0) {
        // Negate them
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      // Add this to the tet gradient
      tetgrad[0] += nx;
      tetgrad[1] += ny;
      tetgrad[2] += nz;
    }
    // Third vertex
    {
      double dx1, dx2, dx3;
      double dy1, dy2, dy3;
      double dz1, dz2, dz3;
      double nx, ny, nz;
      // These are 1-2, 1-3, and 1-0 respectively
      dx1 = pos[3*t->verts[1]] - pos[3*t->verts[2]];
      dy1 = pos[3*t->verts[1]+1] - pos[3*t->verts[2]+1];
      dz1 = pos[3*t->verts[1]+2] - pos[3*t->verts[2]+2];
      
      dx2 = pos[3*t->verts[1]] - pos[3*t->verts[3]];
      dy2 = pos[3*t->verts[1]+1] - pos[3*t->verts[3]+1];
      dz2 = pos[3*t->verts[1]+2] - pos[3*t->verts[3]+2];
      
      dx3 = pos[3*t->verts[1]] - pos[3*t->verts[0]];
      dy3 = pos[3*t->verts[1]+1] - pos[3*t->verts[0]+1];
      dz3 = pos[3*t->verts[1]+2] - pos[3*t->verts[0]+2];
      
      // We will use the cross product of the latter two, using the first
      // to check orientation
      nx = dy2 * dz3 - dz2 * dy3;
      ny = dz2 * dx3 - dx2 * dz3;
      nz = dx2 * dy3 - dy2 * dx3;
      if (nx * dx1 + ny * dy1 + nz * dz1 < 0) {
        // Negate them
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      // Add this to the tet gradient
      tetgrad[0] += nx;
      tetgrad[1] += ny;
      tetgrad[2] += nz;
    }
    // fourth vertex
    {
      double dx1, dx2, dx3;
      double dy1, dy2, dy3;
      double dz1, dz2, dz3;
      double nx, ny, nz;
      // These are 2-3, 2-0, and 2-1 respectively
      dx1 = pos[3*t->verts[2]] - pos[3*t->verts[3]];
      dy1 = pos[3*t->verts[2]+1] - pos[3*t->verts[3]+1];
      dz1 = pos[3*t->verts[2]+2] - pos[3*t->verts[3]+2];
      
      dx2 = pos[3*t->verts[2]] - pos[3*t->verts[0]];
      dy2 = pos[3*t->verts[2]+1] - pos[3*t->verts[0]+1];
      dz2 = pos[3*t->verts[2]+2] - pos[3*t->verts[0]+2];
      
      dx3 = pos[3*t->verts[2]] - pos[3*t->verts[1]];
      dy3 = pos[3*t->verts[2]+1] - pos[3*t->verts[1]+1];
      dz3 = pos[3*t->verts[2]+2] - pos[3*t->verts[1]+2];
      
      // We will use the cross product of the latter two, using the first
      // to check orientation
      nx = dy2 * dz3 - dz2 * dy3;
      ny = dz2 * dx3 - dx2 * dz3;
      nz = dx2 * dy3 - dy2 * dx3;
      if (nx * dx1 + ny * dy1 + nz * dz1 < 0) {
        // Negate them
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      // Add this to the tet gradient
      tetgrad[0] += nx;
      tetgrad[1] += ny;
      tetgrad[2] += nz;
    }
    double grad_len = len(tetgrad[0], tetgrad[1], tetgrad[2]);
    if (grad_len == 0)
      continue;
    tetgrad[0] /= grad_len;
    tetgrad[1] /= grad_len;
    tetgrad[2] /= grad_len;
    // Now the gradient is normalized

    // Use it to update vertex divergences
    // The vertex divergence is the sum over incident tets of
    // (1/3|A_t|) sum ( (|e*|/|e|) (e dot X) )
    // where X is the normalized gradient over the tet, and the
    // sum is taken over edges e of the tet originating at the vertex.
    
    double dx, dy, dz;
    // Edges 0-1 and 1-0
    int e = map_vert_to_edge(t->verts[0], t->verts[1]);
    dx = pos[3 * t->verts[0]] - pos[3 * t->verts[1]];
    dy = pos[3 * t->verts[0] + 1] - pos[3 * t->verts[1] + 1];
    dz = pos[3 * t->verts[0] + 2] - pos[3 * t->verts[1] + 2];
    double e_len = len(dx, dy, dz);
    double e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[0]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len)) * e_dot_x;
    vdiv[t->verts[1]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len)) * e_dot_x;

    // Edges 0-2 and 2-0
    dx = pos[3 * t->verts[0]] - pos[3 * t->verts[2]];
    dy = pos[3 * t->verts[0] + 1] - pos[3 * t->verts[2] + 1];
    dz = pos[3 * t->verts[0] + 2] - pos[3 * t->verts[2] + 2];
    e_len = len(dx, dy, dz);
    e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[0]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    vdiv[t->verts[2]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    
    // Edges 0-3 and 3-0
    dx = pos[3 * t->verts[0]] - pos[3 * t->verts[3]];
    dy = pos[3 * t->verts[0] + 1] - pos[3 * t->verts[3] + 1];
    dz = pos[3 * t->verts[0] + 2] - pos[3 * t->verts[3] + 2];
    e_len = len(dx, dy, dz);
    e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[0]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    vdiv[t->verts[3]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));

    // Edges 1-2 and 2-1
    dx = pos[3 * t->verts[1]] - pos[3 * t->verts[2]];
    dy = pos[3 * t->verts[1] + 1] - pos[3 * t->verts[2] + 1];
    dz = pos[3 * t->verts[1] + 2] - pos[3 * t->verts[2] + 2];
    e_len = len(dx, dy, dz);
    e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[1]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    vdiv[t->verts[2]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));

    // Edges 1-3 and 3-1
    dx = pos[3 * t->verts[1]] - pos[3 * t->verts[3]];
    dy = pos[3 * t->verts[1] + 1] - pos[3 * t->verts[3] + 1];
    dz = pos[3 * t->verts[1] + 2] - pos[3 * t->verts[3] + 2];
    e_len = len(dx, dy, dz);
    e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[1]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    vdiv[t->verts[3]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));

    // Edges 2-3 and 3-2
    dx = pos[3 * t->verts[2]] - pos[3 * t->verts[3]];
    dy = pos[3 * t->verts[2] + 1] - pos[3 * t->verts[3] + 1];
    dz = pos[3 * t->verts[2] + 2] - pos[3 * t->verts[3] + 2];
    e_len = len(dx, dy, dz);
    e_dot_x = dx * tetgrad[0] + dy * tetgrad[1] + dz * tetgrad[2];
    vdiv[t->verts[2]] -= (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
    vdiv[t->verts[3]] += (1 / (12 * t->volume)) * (0.75 * edge_areas[e] / (e_len * e_len));
  }

  // Now that we're done with divergences, solve the Poisson problem
  cholmod_free_dense(&vec2, workspace); // don't need it anymore

  vec2 = cholmod_solve(CHOLMOD_A, L2, vec3, workspace);
  memcpy(dists, vec2->x, nverts * sizeof(double));

  double min = dists[0];
  for (int i = 1; i < nverts; ++i) {
    if (dists[i] < min)
      min = dists[i];
  }
  for (int i = 0; i < nverts; ++i) {
    dists[i] -= min;
  }

  // And free everything
  cholmod_free_dense(&vec1, workspace);
  cholmod_free_dense(&vec2, workspace);
  cholmod_free_dense(&vec3, workspace);
  return;
}
