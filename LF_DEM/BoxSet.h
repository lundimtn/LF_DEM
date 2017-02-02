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
#include "vec3d.h"
#include "LeesEdwards.h"
#include "Box.h"

#define DELETE(x) if(x){delete [] x; x = NULL;}
class System;

class BoxSet{
private:
	double box_xsize;
	double box_ysize;
	double box_zsize;
	std::size_t x_box_nb;
	std::size_t y_box_nb;
	std::size_t z_box_nb;
	std::size_t box_nb;
	bool _is_boxed;
	std::set <Box*> Boxes;
	std::set <Box*> BulkBoxes;
	std::set <Box*> TopBoxes;
	std::set <Box*> BottomBoxes;
	std::set <Box*> TopBottomBoxes;
	std::vector <Box*> box_labels;
	Box* whichBox(const vec3d&);
	void updateNeighbors(const LeesEdwards &);
	// init methods
	void allocateBoxes();
	void positionBoxes();
	void buildProbingPositions();
	void assignNeighbors(const LeesEdwards &);
	void assignNeighborsBulk(const LeesEdwards &);
	void assignNeighborsTop(const LeesEdwards &);
	void assignNeighborsBottom(const LeesEdwards &);
	void assignNeighborsTopBottom(const LeesEdwards &);
	std::vector<Box*> boxMap;
	std::vector<vec3d> top_probing_positions;
	std::vector<vec3d> bottom_probing_positions;
public:
	void init(double interaction_dist,
	          LeesEdwards pbc,
	          unsigned int np);
	// void init(double interaction_dist, System *sys_);
	/*****
	 update()

	 To be called at each time step.
	 It updates the neighborhood relations betwenn boxes.
	 Those relations change at each time step for boxes on top or bottom
	 *****/
	void update(const LeesEdwards &);

	/*****
	 box(int i)
	 boxes particles i
	 should be called after moving particle i
	 *****/
	void box(int i, const vec3d &pos);
	/*****
	 neighborhood_begin(int i) and neighborhood_end(int i)
	 gives iterators to beginning and ending point of the container including
	 all particles in the box containing particle i and in the adjacent boxes.
	 *****/
	const std::vector <int>& neighborhood(int i);
	void inflateZ(double inflation_ratio);
	void printBoxNetwork();
	void printBoxContainers();
	void printNeighborhoodContainers();
	void printBoxMap();
};
#endif /* defined(__LF_DEM__BoxSet__) */
