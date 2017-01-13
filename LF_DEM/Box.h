//
//  Box.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Box
 \brief Box object holding identities of the particles in a subset of the whole suspension, to be used in a BoxSet object.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Box__
#define __LF_DEM__Box__
#define __vector_container__
#ifndef __vector_container__
#define __flist_container__
#endif
#include "vec3d.h"
#include <set>
#include <vector>

class Box{
private:
	std::vector <Box*> static_neighbors;
	std::vector <Box*> moving_neighbors;
	std::set <int> container;
	std::vector <int> neighborhood_container;
	vec3d position;

public:
	Box(){};
	~Box();

	void addStaticNeighbor(Box* neigh_box);
	void addMovingNeighbor(Box* neigh_box);
	void resetMovingNeighbors();
	const std::vector <Box*> & getNeighborBox(){return static_neighbors;}

	vec3d getPosition() const {return position;}
	void setPosition(vec3d pos) {position = pos;}

	void add(int);
	void remove(int);

	const std::set <int> & getContainer() const {return container;}
	const std::vector <int> & getNeighborhoodContainer() const {return neighborhood_container;};
	void buildNeighborhoodContainer();
};

#endif /* defined(__LF_DEM__Box__) */
