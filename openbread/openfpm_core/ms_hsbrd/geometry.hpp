#include "math.h"
#include <cstdlib>


#ifndef GEOMETRY_H
#define GEOMETRY_H

// random_unit_vector
// Generates a random 3D unit vector (direction) with a uniform spherical distribution
// Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
Point<3,double> random_unit_vector()
{
	Point<3,double> vect;

	double phi = 2.0*M_PI*(double)rand()/RAND_MAX;
	double theta = M_PI*(double)rand()/RAND_MAX;
	//double theta = acos(costheta);

	vect[0] = sin(theta)*cos(phi);
	vect[1] = sin(theta)*sin(phi);
	vect[2] = cos(theta);

	return vect;
}


#endif //GEOMETRY_H
