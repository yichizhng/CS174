#ifndef MAIN_H_
#define MAIN_H_

#include <cstddef> // For size_t
#include <functional>

// Class definition

void boundset();
void edgesetup();
void findboundary();

class Edge {
public:
  // We will guarantee that v1 < v2
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
      return (hash<int>()(e.v1) ^ hash<int>()(e.v2));
    }
  };
}

class Triangle {
public:
  // Constructor
  Triangle(int v1, int v2, int v3 ){
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
  };
  
  int vertices[3];
  // All other data can be accessed from those variables
  
  // Auxiliary data (reusable calculations)
  double area; // Triangle area
  double cotans[3]; // Cotangents of angles
};

#endif /* MAIN_H_ */
