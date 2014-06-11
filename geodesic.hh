#ifndef GEODESIC_H_
#define GEODESIC_H_

#include <vector>

void step0();

void timestep(double *u);

void setverts(std::vector<int>&);
void step1(double *dists);

#endif /* GEODESIC_H_ */
