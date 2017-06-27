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
// #define __vector_container__
// #ifndef __vector_container__
// #define __flist_container__
// #endif
#include "vec3d.h"
#include <set>
#include <vector>

class Box{
private:
	std::vector <Box*> neighbors;
	std::set <unsigned> container;
	std::vector <unsigned>	neighborhood_container;

public:
	void addNeighborBox(Box* neigh_box);
	void add(unsigned);
	void addFromNeighbor(unsigned);

	void remove(unsigned);
	void removeFromNeighbor(unsigned);
	void buildNeighborhoodContainer();

	const std::vector <unsigned> & getNeighborhoodContainer() const {return neighborhood_container;}
	const std::set <unsigned> & getContainer() const {return container;}
};

#endif /* defined(__LF_DEM__Box__) */
