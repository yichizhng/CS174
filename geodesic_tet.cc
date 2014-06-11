// The implementation of the geodesic distance finder on a
// 3-simplicial manifold

// The steps are roughly analogous to those of the 2-simplicial manifold:

// Integrate the heat flow for a fixed time t (we can use the same one, really)

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

// Data structures:

// Constants (well, per mesh)
extern int ntets;
extern int nedges;
extern int nverts;
extern Tet **tets;
extern double *pos;

#define SMALL_FLOAT 1000000*DBL_EPSILON

static cholmod_common *workspace;
cholmod_sparse *mat1, *mat2, *mat3;

static inline double len(double x, double y, double z) {
  return sqrt( (x*x) + (y*y) + (z*z) );
}

using namespace std;

// An unordered map, allowing us to retrieve an edge index from vertex pairs,
// although it should not be used directly
unordered_map<Edge, int> vert_edge_map;
// Note that nedges == vert_edge_map.size()

// Returns the edge corresponding to the pair (i,j) should it exist; -1
// otherwise
// Given that we store the mesh as tets, there should be no reason to
// make a call to this that returns -1, but we should try to be safe.
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

// Sum of incident tet volumes on an edge
double *edge_areas;

// Sum of incident tet volumes on a vertex
double *vert_areas;

void init() {
  // Set up things needed for calculations in this file; assumes all extern
  // variables have already been set
  vert_areas = (double *) calloc(nverts, sizeof(double));
  edge_areas = (double *) calloc(nedges, sizeof(double));
  workspace = (cholmod_common *)malloc(sizeof(cholmod_common));
  cholmod_start(workspace);
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
    y2 = pos[3 * t->verts[2]] - pos[3 * t->verts[0]];
    y3 = pos[3 * t->verts[2]] - pos[3 * t->verts[0]];
    z1 = pos[3 * t->verts[3]] - pos[3 * t->verts[0]];
    z2 = pos[3 * t->verts[3]] - pos[3 * t->verts[0]];
    z3 = pos[3 * t->verts[3]] - pos[3 * t->verts[0]];

    t->volume = ( (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2) -
                     (x1 * y3 * z2 + x2 * y1 * z3 + x3 * y2 * z1) ) / 6;
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

  // (i != j)
  // L_(i,j) = (sum over tetrahedra t incident on both i and j)
  // |e*_ij| / |e_ij|
  // |e_ij| is obvious (it's just the length of edge ij) but
  // |e*_ij| isn't; we'll approximate it as edge_areas[ij]
  // Meanwhile L_(i,i) is (sum over j adjacent to i) -L_(i,j)

  cholmod_triplet *trip1;
  // trip1 will be used to build L_c
  
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
        dy = pos[3*v2+1] - pos[3*v2+1];
        dz = pos[3*v2+2] - pos[3*v2+2];
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
}
