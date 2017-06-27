//
//  BoxSet.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class BoxSet
 \brief Set of Box objects making a partition of the simulation box.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__BoxSet__
#define __LF_DEM__BoxSet__
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "PeriodicBoundary.h"
#include "vec3d.h"
#include "Box.h"

#define DELETE(x) if(x){delete [] x; x = NULL;}
class System;


class BoxSet{
private:
	std::array<std::size_t, 3> box_nb;
	PBC::RhomboidLattice rhomb;
	std::array<vec3d, 3> rhomb_unit_normals;
	std::array<double, 3> rhomb_depths;

	bool _is_boxed;
	double _box_min_size;
	std::vector<Box> boxes;
	std::vector<Box*> boxMap;
	/*****
	 WhichBox(vec3d pos)
	 returns a pointer on the box containg position pos
	 *****/
	Box* whichBox(const vec3d&);

	// init methods
	void assignNeighbors();
	std::array<std::size_t, 3> computeBoxNumber() const;
	void setupBoxes();
public:
	void init(double box_min_size,
	          const PBC::RhomboidLattice &rhomb_,
	          std::size_t np);
	bool setPeriodicBoundaryConditions(const PBC::RhomboidLattice &rhomb_);

	/*****
	 box(unsigned i, vec3d position_i)
	 boxes particles i
	 should be called after moving particle i
	 *****/
	void box(unsigned i, vec3d position_i);
	void update();
	/*****
	 neighborhood_begin(int i) and neighborhood_end(int i)
	 gives iterators to beginning and ending point of the container including
	 all particles in the box containing particle i and in the adjacent boxes.
	 *****/
	const std::vector <unsigned>& neighborhood(unsigned i) const;
	// void printBoxNetwork();
	// void printBoxContainers();
	// void printNeighborhoodContainers();
	// void printBoxMap();
	// void yaplotBox(std::ofstream &fout_boxing);
};
#endif /* defined(__LF_DEM__BoxSet__) */
