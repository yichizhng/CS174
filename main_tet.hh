#ifndef MAIN_TET_H
#define MAIN_TET_H

#include <cstddef>
#include <functional>

class Edge {
public:
  int v1;
  int v2;
  Edge(int a, int b) : v1(a), v2(b) {};
  
  bool operator==(const Edge &other) const {
    return (v1 == other.v1) && (v2 == other.v2);
  }
};

namespace std {
  template <>
  struct hash<Edge> {
    size_t operator()(const Edge& e) const {
      return (hash<int>()(e.v1) ^ hash<int>()(-e.v2));
    }
  };
}

class Tet {
public:
  Tet(int v1, int v2, int v3, int v4) {
    verts[0] = v1;
    verts[1] = v2;
    verts[2] = v3;
    verts[3] = v4;
  };
  int verts[4];
  double volume;
};

#endif /* MAIN_TET_H */
