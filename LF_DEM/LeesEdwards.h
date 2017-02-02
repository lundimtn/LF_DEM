#ifndef __LF_DEM__LeesEdwards__
#define __LF_DEM__LeesEdwards__
#include <assert.h>
#include "vec3d.h"

class LeesEdwards {
public:
	void init(vec3d system_dimensions,
	          vec3d shear_displacement,
	          bool keep_init_strain);
	void apply_strain(vec3d strain_increment_zortho);
	int periodize(vec3d&) const;
	int periodizeDiff(vec3d&) const;
	vec3d periodized(const vec3d&) const;
	vec3d dimensions() const {return L;}
	vec3d shear_disp() const {return shear_disp_;};
	double cumulated_strain() const {return cumulated_shear_strain_;};
	vec3d shear_strain() const {return shear_strain_;};
private:
	vec3d L;
	vec3d Lhalf;
	vec3d shear_disp_;
	double cumulated_shear_strain_;
	vec3d shear_strain_;
};

inline void LeesEdwards::init(vec3d system_dimensions,
                              vec3d shear_displacement,
                              bool keep_init_strain)
{
	L = system_dimensions;
	Lhalf = L/2;
	assert(shear_displacement.z == 0);
	shear_disp_ = shear_displacement;
	if (keep_init_strain) {
		shear_strain_ = shear_displacement/L.z;
	} else {
		shear_strain_ = {0, 0, 0};
	}
	cumulated_shear_strain_ = 0;
}

inline void LeesEdwards::apply_strain(vec3d strain_increment_zortho)
{
	// homothety along z
	// we do it before shear because it affects it through L.z
	L.z += L.z*strain_increment_zortho.z;
	Lhalf = L/2;
	// shear: flow in x/y, gradient z
	vec3d shear_increment ({strain_increment_zortho.x,
													strain_increment_zortho.y,
													0});
	shear_disp_ += shear_increment*L.z;
	int m = (int)(shear_disp_.x/L.x);
	if (shear_disp_.x < 0) {
		m--;
	}
	shear_disp_.x = shear_disp_.x-m*L.x;
	m = (int)(shear_disp_.y/L.y);
	if (shear_disp_.y < 0) {
		m--;
	}
	shear_disp_.y = shear_disp_.y-m*L.y;
	cumulated_shear_strain_ += shear_increment.norm();
	shear_strain_ += shear_increment;
}

// [0,l]
inline int LeesEdwards::periodize(vec3d& pos) const
{
	int z_shift = 0;
	if (pos.z >= L.z) {
		pos.z -= L.z;
		pos -= shear_disp_;
		z_shift = -1;
	} else if (pos.z < 0) {
		pos.z += L.z;
		pos += shear_disp_;
		z_shift = 1;
	}
	while (pos.x >= L.x) {
		pos.x -= L.x;
	}
	while (pos.x < 0) {
		pos.x += L.x;
	}
	if (pos.y >= L.y) {
		pos.y -= L.y;
	} else if (pos.y < 0) {
		pos.y += L.y;
	}
	return z_shift;
}

// [0,l]
inline vec3d LeesEdwards::periodized(const vec3d &pos) const
{
	vec3d periodized_pos = pos;
	periodize(periodized_pos);
	return periodized_pos;
}

inline int LeesEdwards::periodizeDiff(vec3d& pos_diff) const
{
	/** Periodize a separation vector with Lees-Edwards boundary condition

		On return pos_diff is the separation vector corresponding to the closest copies,
		and zshift is the shift applied in the z direction (-1, 0 or +1)
	 */
	int zshift = 0;
	if (pos_diff.z > Lhalf.z) {
		pos_diff.z -= L.z;
		pos_diff -= shear_disp_;
		zshift = -1;
	} else if (pos_diff.z < -Lhalf.z) {
		pos_diff.z += L.z;
		pos_diff += shear_disp_;
		zshift = 1;
	}
	while (pos_diff.x > Lhalf.x) {
		pos_diff.x -= L.x;
	}
	while (pos_diff.x < -Lhalf.x) {
		pos_diff.x += L.x;
	}
	if (pos_diff.y > Lhalf.y) {
		pos_diff.y -= L.y;
	} else if (pos_diff.y < -Lhalf.y) {
		pos_diff.y += L.y;
	}
	return zshift;
}
#endif /* defined(__LF_DEM__LeesEdwards__) */
