#ifndef __LF_DEM__LeesEdwards__
#define __LF_DEM__LeesEdwards__
#include "vec3d.h"

class LeesEdwards {
public:
	void set(vec3d system_dimensions,
	         vec3d shear_displacement) {
		L = system_dimensions;
		Lhalf = L/2;
		shear_disp_ = shear_displacement;
	}
	int periodize(vec3d&) const;
	int periodizeDiff(vec3d&) const;
	vec3d periodized(const vec3d&) const;
	vec3d dimensions() const {return L;}
	vec3d shear_disp() const {return shear_disp_;};
private:
	vec3d L;
	vec3d Lhalf;
	vec3d shear_disp_;
};

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
