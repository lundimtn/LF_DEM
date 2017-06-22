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
	Rhomboid rhomb;
	std::array<vec3d, 3> rhomb_unit_normals;
	bool _is_boxed;
	double _box_min_size;
	std::vector <Box> boxes;
	std::vector<unsigned> boxMap;
	/*****
	 WhichBox(vec3d pos)
	 returns a pointer on the box containg position pos
	 *****/
	Box* whichBox(const vec3d&) const;
	Box* whichBox(unsigned int box_label) const;

	// init methods
	void assignNeighbors(const PBC::PeriodicBoundary &pb);
	bool setPeriodicBoundaryConditions(const PBC::PeriodicBoundary &pb);
	std::array<std::size_t, 3> computeBoxNumber() const;
	void setupBoxes(const PBC::PeriodicBoundary &pb);
public:
	BoxSet(double box_min_size,
	       const PBC::PeriodicBoundary &pb,
	       std::size_t np)
	/*****
	 is_boxed()

	 Can be called before calling an other method of BoxSet.
	 If false, than calls to other method may usually be avoided.
	 They can be performed anyway though, and they are normally safe (and useless)

	 is_boxed() tells if the boxing is effective.
	 If the system size is small, the neighborhood of a box may contain
	 the whole system. If this happens, the boxing consists of only one box (as it
	 is useless, if not ill defined, to do something else), and is_boxed() returns false.

	 *****/
	bool is_boxed() const;

	/*****
	 box(int i)
	 boxes particles i
	 should be called after moving particle i
	 *****/
	void box(unsigned i);
	void update();
	/*****
	 neighborhood_begin(int i) and neighborhood_end(int i)
	 gives iterators to beginning and ending point of the container including
	 all particles in the box containing particle i and in the adjacent boxes.
	 *****/
	const std::vector <unsigned>& neighborhood(unsigned i);
	void printBoxNetwork();
	void printBoxContainers();
	void printNeighborhoodContainers();
	void printBoxMap();
	void yaplotBox(std::ofstream &fout_boxing);
};
#endif /* defined(__LF_DEM__BoxSet__) */
