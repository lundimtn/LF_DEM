#ifndef __LF_DEM__LeesEdwards__
#define __LF_DEM__LeesEdwards__
#include <assert.h>
#include "vec3d.h"
#include "Sym2Tensor.h"

namespace PBC {

class RhomboidLattice {
private:
	std::array<vec3d, 3> edges;                  // edges of the rhomboid box
	std::array<vec3d, 3> side_normals;           // unit vectors normal to the box sides
	std::array<vec3d, 3> GS_basis;               // normal of the unslanted box, for nearest image separation (Graham-Schmidt)
	std::array<double, 3> depths;                // the dimensions of the unslanted box
	std::array<unsigned, 3> ordering;            // ordering of the projection on GS basis for nearest image calc
	void setEdges(const std::array<vec3d, 3> &_edges)

public:
	void applyStrain(const matrix &strain_tensor);
	inline std::array<double, 3> placeInPrimaryCell(vec3d& pos) const;
	inline vec3d primaryCellImage(const vec3d &pos) const;
	inline std::array<double, 3> getNearestImage(vec3d& separation) const;
	inline std::array<double, 3> getNearestImageLargeSystem(vec3d& separation) const;
}


void RhomboidLattice::applyStrain(const matrix &strain_tensor) {
	vec3d new_edges [3];
	for (unsigned i=0; i<3; i++) {
		new_edges = strain_tensor*edges;
	}
	setEdges(new_edges);
}


void RhomboidLattice::setEdges(const std::array<vec3d, 3> &_edges)
{
	edges = _edges;
	side_normals = {unitvector(cross(edges[1], edges[2])),
					unitvector(cross(edges[2], edges[0])),
					unitvector(cross(edges[0], edges[1]))};
	depths = {dot(edges[0], cross(edges[1], edges[2]).unitvector()),
			  dot(edges[1], cross(edges[2], edges[0]).unitvector()),
			  dot(edges[2], cross(edges[0], edges[1]).unitvector())};
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

	// constructing the Graham-Schmidt basis
	GS_basis[ordering[0]] = side_normals[ordering[0]];
	GS_basis[ordering[1]] = side_normals[ordering[1]]
							- dot(side_normals[ordering[1]], side_normals[ordering[0]])*side_normals[ordering[0]];
	GS_basis[ordering[1]].unitvector();
	GS_basis[ordering[2]] = side_normals[ordering[2]];
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
	periodize(periodized_pos);
	return periodized_pos;
}

inline std::array<double, 3> RhomboidLattice::getNearestImageLargeSystem(vec3d& separation) const
{
	// Faster but not correct for boxes smaller than twice the interaction range.
	std::array<double, 3> box_shifts;
	for (auto u: ordering) {
		box_shifts[u] = std::round(dot(side_normals[u], pos)/depths[u]);
		pos -= box_shifts[u]*edges[u];
	}
	return box_shifts;
}


class LeesEdwards { // slanted boxes, not sliding boxes, gradient along z
public:
	void init(vec3d system_dimensions,
	          vec3d shear_displacement,
	          bool keep_init_strain);
	void applyStrain(vec3d strain_increment_zortho);
	vec3d dimensions() const {return L;}
	double getCumulatedStrain() const {return cumulated_shear_strain_;};
	matrix getStrain() const {return shear_strain_;};
	Rhomboid getRhomboid() const {return box_rhomboid;};
private:
	vec3d L;
	vec3d Lhalf;
	double cumulated_strain_;
	matrix strain_;
	RhomboidLattice box_rhomboid;                       // for particle postions
};



inline void LeesEdwards::init(vec3d system_dimensions,
                              vec3d shear_displacement,
                              bool keep_init_strain)
{
	L = system_dimensions;
	Lhalf = L/2;
	assert(shear_displacement.z == 0);
	box_rhomboid.edges[0] = {L.x, 0, 0};
	box_rhomboid.edges[1] = {0, L.y, 0};
	box_rhomboid.edges[2] = shear_displacement;

	if (keep_init_strain) {
		strain_[2] = shear_displacement.x/L.z;
		strain_[5] = shear_displacement.y/L.z;
	}
	cumulated_strain_ = 0;
}

inline void LeesEdwards::apply_strain(matrix strain_tensor)
{
	 // make sure there is no shear with gradient along x or y
	assert(strain_tensor.elm[1]==0);
	assert(strain_tensor.elm[3]==0);
	assert(strain_tensor.elm[6]==0);
	assert(strain_tensor.elm[7]==0);

	box_rhomboid.apply_strain(strain_tensor);
	strain_ += strain_tensor;
	cumulated_strain_ += strain_tensor.norm();

	auto &shear_disp = box_rhomboid.edges[2];
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
