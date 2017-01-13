#include <stdexcept>
#include <sstream>
#include "BoxSet.h"
using namespace std;

void BoxSet::init(double interaction_dist,
                  LeesEdwards pbc,
                  unsigned int np)
{
	string indent = "  BoxSet::\t";
	cout << indent << "Setting up Cell List System ... ";
	boxMap.resize(np);
	for (auto &bx: boxMap) {
		bx = NULL;
	}
	auto system_dimensions = pbc.dimensions();
	double xratio = system_dimensions.x/interaction_dist;
	double yratio = system_dimensions.y/interaction_dist;
	double zratio = system_dimensions.z/interaction_dist;
	x_box_nb = (unsigned int)xratio;
	y_box_nb = (unsigned int)yratio;
	z_box_nb = (unsigned int)zratio;
	if (x_box_nb == 0) {
		x_box_nb = 1;
	}
	if (y_box_nb == 0) {
		y_box_nb = 1;
	}
	if (z_box_nb == 0) {
		z_box_nb = 1;
	}
	if (x_box_nb < 4 && y_box_nb < 4 && z_box_nb < 4) { // boxing useless: a neighborhood is the whole system
		_is_boxed = false;
		box_xsize = system_dimensions.x;
		box_ysize = system_dimensions.y;
		box_zsize = system_dimensions.z;
		box_nb = 1;

		auto it = Boxes.insert(new Box());
		Box* const b = (*it.first);
		b->position = {0, 0, 0};
		TopBottomBoxes.insert(b);
		box_labels.push_back(b);
	} else {
		_is_boxed = true;
		box_xsize = system_dimensions.x/x_box_nb;
		box_ysize = system_dimensions.y/y_box_nb;
		box_zsize = system_dimensions.z/z_box_nb;
		int m1p1[] = {-1, 1};
		for (int a : m1p1) {
			for (int b : m1p1) {
				auto far_corner = 1.4999999*vec3d(a*box_xsize, b*box_ysize, box_zsize);
				top_probing_positions.push_back(far_corner);
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize, 0, 0));
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize, b*box_ysize, 0));
				top_probing_positions.push_back(far_corner-vec3d(0, b*box_ysize, 0));

				far_corner = 1.4999999*vec3d(a*box_xsize, b*box_ysize, -box_zsize);
				bottom_probing_positions.push_back(far_corner);
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize, 0, 0));
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize, b*box_ysize, 0));
				bottom_probing_positions.push_back(far_corner-vec3d(0, b*box_ysize, 0));
			}
		}
		box_nb = x_box_nb*y_box_nb*z_box_nb;
		allocateBoxes();
		// give them their position
		positionBoxes();
		// tell them their neighbors
		assignNeighbors(pbc);
	}
	cout << " [ok]" << endl;
}

void BoxSet::inflateZ(double inflation_ratio)
{
	for (auto &pos: top_probing_positions) {
		pos.z *= inflation_ratio;
	}
	for (auto &pos: bottom_probing_positions) {
		pos.z *= inflation_ratio;
	}
	box_zsize *= inflation_ratio;
	for (auto &bx: Boxes) {
		bx->position.z *= inflation_ratio;
	}
}

void BoxSet::allocateBoxes()
{
	for (unsigned int i=0; i<box_nb; i++) {
		Boxes.insert(new Box());
	}
	box_labels.resize(box_nb);
}

void BoxSet::positionBoxes()
{
	// position boxes
	auto it = Boxes.begin();

	for (unsigned int ix=0; ix<x_box_nb; ix++) {
		for (unsigned int iy=0; iy<y_box_nb; iy++) {
			for (unsigned int iz=0; iz<z_box_nb; iz++) {
				Box* const bx = (*it);
				bx->position = {box_xsize*(ix+0.5), box_ysize*(iy+0.5), box_zsize*(iz+0.5)}; // the center of the box
				int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;
				box_labels[label] = bx;
				if (iz == 0 && iz < z_box_nb-1) {// bottom box
					BottomBoxes.insert(bx);
				}
				if (iz == z_box_nb-1 && iz > 0) {//top box
					TopBoxes.insert(bx);
				}
				if (iz == 0 && iz == z_box_nb-1) {// bottom box
					TopBottomBoxes.insert(bx);
				}
				if (iz > 0 && iz < z_box_nb-1) {// bulk box
					BulkBoxes.insert(bx);
				}
				it++;
			}
		}
	}
}


void BoxSet::assignNeighborsBulk(const LeesEdwards &pbc)
{
	for (auto& bx : BulkBoxes) {
		vec3d delta;
		int m10p1[] = {-1, 0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10p1) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(pbc.periodized(bx->position+delta)));
				}
			}
		}
	}
}

void BoxSet::assignNeighborsBottom(const LeesEdwards &pbc)
{
	for (auto& bx : BottomBoxes) {
		vec3d delta;
		// boxes  at same level and above first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		int p10[] = {0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : p10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(pbc.periodized(bx->position+delta)));
				}
			}
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}
}

void BoxSet::assignNeighborsTop(const LeesEdwards &pbc)
{
	for (auto& bx : TopBoxes) {
		vec3d delta;
		// boxes  at same level and bottom first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		int m10[] = {-1, 0};
		for (const auto & a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(pbc.periodized(bx->position+delta)));
				}
			}
		}
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}
}

void BoxSet::assignNeighborsTopBottom(const LeesEdwards &pbc)
{
	for (auto& bx : TopBottomBoxes) {
		vec3d delta;

		// boxes at same level first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				delta.z = 0;
				bx->addStaticNeighbor(whichBox(pbc.periodized(bx->position+delta)));
			}
		}

		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}
}

void BoxSet::assignNeighbors(const LeesEdwards &pbc)
{
	assignNeighborsBulk(pbc);
	assignNeighborsBottom(pbc);
	assignNeighborsTop(pbc);
	assignNeighborsTopBottom(pbc);
}

/*****
 UpdateNeighbors()

 At each time step, we need to check if neighborhood on top and bottom boxes have changed.
 Bulk boxes do not need to be updated.
 *****/
void BoxSet::updateNeighbors(const LeesEdwards &pbc)
{
	/**
	 \brief Update the neighbors of top and bottom boxes have changed.

		To be called when the boundary conditions have changed.
	 **/

	for (auto& bx : TopBoxes) {
		bx->resetMovingNeighbors();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}

	for (auto& bx : BottomBoxes) {
		bx->resetMovingNeighbors();
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}

	for (auto& bx : TopBottomBoxes) {
		bx->resetMovingNeighbors();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(pbc.periodized(bx->position+delta_prob)));
		}
	}
}

//public methods
void BoxSet::update(const LeesEdwards &pbc)
{
	if (_is_boxed) {
		updateNeighbors(pbc);
	}
	for (const auto& bx : Boxes) {
		bx->buildNeighborhoodContainer();
	}
}

Box* BoxSet::whichBox(const vec3d &pos)
{
	unsigned int ix = (unsigned int)(pos.x/box_xsize);
	unsigned int iy;
	if (box_ysize > 0) {
		iy = 0;
	} else {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	unsigned int iz = (unsigned int)(pos.z/box_zsize);
	unsigned int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;
	if (label > box_labels.size()-1) {
		ostringstream error_str;
		error_str  << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	return box_labels[label];
}

void BoxSet::box(int i, const vec3d &pos)
{
	Box* b = whichBox(pos);
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

const vector<int>& BoxSet::neighborhood(int i){
	return (boxMap[i])->getNeighborhoodContainer();
}

void BoxSet::printBoxNetwork()
{
	for (const auto& bx : Boxes) {
		const auto& neighbors = bx->getNeighborBox();
		for (const auto& neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->position << " ";
			cerr << neighbor_box->position << endl;
		}
	}
}

void BoxSet::printBoxContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getContainer()) {
			cerr << bx->position << " " << j << endl;
		}
	}
}

void BoxSet::printNeighborhoodContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getNeighborhoodContainer()) {
			cerr << bx->position << " " << j << endl;
		}
	}
}

void BoxSet::printBoxMap()
{
	for (unsigned int i=0; i<boxMap.size(); i++) {
		cerr << i << " " << boxMap[i]->position << endl;
	}
}
