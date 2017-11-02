#include "TimeActivatedAdhesion.h"

namespace TActAdhesion {


void TimeActivatedAdhesion::update(double time_now, double gap, vec3d &nvec)
{
	if (gap > params.adhesion_range) {
		if (state.activity != TAAActivity::inactive) {
			state.activity = TAAActivity::inactive;
			force_amplitude = 0;
		}
	} else {
		switch (state.activity) {
			case TAAActivity::inactive:
				state.activity = TAAActivity::dormant;
				state.contacting_time = time_now;
				// no break, as we want to allow for immediate activation
			case TAAActivity::dormant:
				if ((time_now - state.contacting_time) >= params.activation_time) {
					state.activity = TAAActivity::active;
				}
			default:
				break;

		}
		if (state.activity == TAAActivity::active) {
			if (gap > 0) {// we also know gap < range 
				force_amplitude = params.adhesion_max_force*gap/params.adhesion_range; 
			} else {
				force_amplitude = 0;
			}
			force_on_p0 = -force_amplitude*nvec;
		}
	}
}

void TimeActivatedAdhesion::deactivate() 
{
	state.activity = TAAActivity::inactive;
	force_on_p0 = 0;
	force_amplitude = 0;
}

void TimeActivatedAdhesion::addUpForce(vec3d &force_p0, vec3d &force_p1) const 
{
	if (state.activity == TAAActivity::active) {
		force_p0 += force_on_p0;
		force_p1 -= force_on_p0;
	}
}

void TimeActivatedAdhesion::addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, const vec3d &rvec) const
{
	auto sc = outer_sym(rvec, force_on_p0);
	stress_p0 += stress_split_p0*sc;
	stress_p1 += (1-stress_split_p0)*sc;
}

} // namespace TActAdhesion