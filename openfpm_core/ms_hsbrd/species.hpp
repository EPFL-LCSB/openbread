#include <cstdlib>
#include <map>
#include <vector>
#include <string>

#ifndef SPECIES_H
#define SPECIES_H

// Species object
struct Species
{
	int id;
	double radius;
	double diff;
	double mass;

	Species(int _id, double _radius, double _diff, double _mass)
		{id = _id; radius = _radius; diff = _diff; mass = _mass;}
	~Species(){}
};


#endif //SPECIES_H
