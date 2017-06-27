#ifndef __LF_DEM__LeesEdwards__
#define __LF_DEM__LeesEdwards__
#include <assert.h>
#include <numeric>
#include <algorithm>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "Matrix.h"

namespace PBC {

class RhomboidLattice {
private:
	std::array<vec3d, 3> edges;                  // edges of the rhomboid box
	std::array<vec3d, 3> side_normals;           // unit vectors normal to the box sides
	std::array<double, 3> depths;                // the dimensions of the unslanted box
	std::array<unsigned, 3> ordering;            // ordering of the projection on GS basis for nearest image calc

public:
	void setEdges(const std::array<vec3d, 3> &_edges);
	void applyStrain(const matrix &strain_tensor);
	std::array<double, 3> placeInPrimaryCell(vec3d& pos) const;
	vec3d primaryCellImage(const vec3d &pos) const;
	std::array<double, 3> getNearestImage(vec3d& separation) const;
	std::array<double, 3> getNearestImageLargeSystem(vec3d& separation) const;
	const std::array<vec3d, 3> & getSideNormals() const {return side_normals;}
	const std::array<double, 3> & getDepths() const {return depths;}
	const std::array<vec3d, 3> & getEdges() const {return edges;}
};


void RhomboidLattice::applyStrain(const matrix &strain_tensor) {
	std::array<vec3d, 3> new_edges;
	for (unsigned i=0; i<3; i++) {
		new_edges[i] = strain_tensor*edges[i];
	}
	setEdges(new_edges);
}


void RhomboidLattice::setEdges(const std::array<vec3d, 3> &_edges)
{
	edges = _edges;
	side_normals = {unitvector(cross(edges[1], edges[2])),
					unitvector(cross(edges[2], edges[0])),
					unitvector(cross(edges[0], edges[1]))};
	depths = {dot(edges[0], unitvector(cross(edges[1], edges[2]))),
			  dot(edges[1], unitvector(cross(edges[2], edges[0]))),
			  dot(edges[2], unitvector(cross(edges[0], edges[1])))};
	for (auto d: depths) {
		assert(d>0);
	}

	// finding an ordering to estimate nearest image
	std::array<double, 3> min_aspect_ratio;
	min_aspect_ratio[0] = std::min(depths[0]/edges[1].norm(), depths[0]/edges[2].norm());
	min_aspect_ratio[1] = std::min(depths[1]/edges[2].norm(), depths[1]/edges[0].norm());
	min_aspect_ratio[2] = std::min(depths[2]/edges[1].norm(), depths[2]/edges[0].norm());
	std::iota(ordering.begin(), ordering.end(), 0);
	std::sort(ordering.begin(), ordering.end(),
			   [&min_aspect_ratio](size_t i1, size_t i2) {return min_aspect_ratio[i1] > min_aspect_ratio[i2];});

}


inline std::array<double, 3> RhomboidLattice::placeInPrimaryCell(vec3d& pos) const
{
	std::array<double, 3> box_shifts;
	for (unsigned i=0; i<edges.size(); i++) {
		box_shifts[i] = std::floor(dot(side_normals[i], pos)/depths[i]);
		pos -= box_shifts[i]*edges[i];
	}
	return box_shifts;
}


inline vec3d RhomboidLattice::primaryCellImage(const vec3d &pos) const
{
	vec3d periodized_pos = pos;
	placeInPrimaryCell(periodized_pos);
	return periodized_pos;
}

inline std::array<double, 3> RhomboidLattice::getNearestImageLargeSystem(vec3d& separation) const
{
	// Closest vector problem in rhomboid lattice.
	// The method below does not always find the correct solutions for separations of order system size.
	//
	// Note 1: there is no simple exact algorithm
	// Note 2: if there is a boxing, this limitation is not a problem, as you ask for nearest images across the whole system.
	std::array<double, 3> box_shifts;
	for (auto u: ordering) {
		box_shifts[u] = std::round(dot(side_normals[u], separation)/depths[u]);
		separation -= box_shifts[u]*edges[u];
	}
	return box_shifts;
}


class LeesEdwards { // slanted boxes, not sliding boxes, gradient along z
public:
	void init(vec3d system_dimensions,
	          vec3d shear_displacement,
	          bool keep_init_strain);
	void applyStrain(const matrix &strain_tensor);
	vec3d dimensions() const {return L;}
	double getCumulatedStrain() const {return cumulated_strain_;};
	matrix getStrain() const {return strain_;};
private:
	vec3d L;
	vec3d Lhalf;
	double cumulated_strain_;
	matrix strain_;
	RhomboidLattice box_rhomboid;
	vec3d shear_disp;
};



inline void LeesEdwards::init(vec3d system_dimensions,
                              vec3d shear_displacement,
                              bool keep_init_strain)
{
	L = system_dimensions;
	Lhalf = L/2;
	assert(shear_displacement.z == 0);

	std::array<vec3d, 3> edges = {{ {L.x, 0, 0},
	{0, L.y, 0},
	shear_displacement }};
	box_rhomboid.setEdges(edges);

	if (keep_init_strain) {
		strain_.elm[2] = shear_displacement.x/L.z;
		strain_.elm[5] = shear_displacement.y/L.z;
	}
	cumulated_strain_ = 0;
}

inline void LeesEdwards::applyStrain(const matrix &strain_tensor)
{
	 // make sure there is no shear with gradient along x or y
	assert(strain_tensor.elm[1]==0);
	assert(strain_tensor.elm[3]==0);
	assert(strain_tensor.elm[6]==0);
	assert(strain_tensor.elm[7]==0);

	box_rhomboid.applyStrain(strain_tensor);
	strain_ += strain_tensor;
	cumulated_strain_ += strain_tensor.norm();

	shear_disp = box_rhomboid.getEdges()[2];
	int m = (int)(shear_disp.x/L.x);
	if (shear_disp.x < 0) {
		m--;
	}
	shear_disp.x -= m*L.x;

	m = (int)(shear_disp.y/L.y);
	if (shear_disp.y < 0) {
		m--;
	}
	shear_disp.y -= m*L.y;
}


}

#endif /* defined(__LF_DEM__LeesEdwards__) */
