#ifndef __LF_DEM__GradU__
#define __LF_DEM__GradU__
#include "Matrix.h"
#include "AntiSym2Tensor.h"
#include "Sym2Tensor.h"

namespace Algebra
{
struct GradU {
	matrix grad_u;
	Sym2Tensor E;
	Sym2Tensor E_shape;
	AntiSym2Tensor O;
	AntiSym2Tensor O_shape;
};

struct GradU buildGradU(matrix grad_u_) {
	struct GradU gu;
	gu.grad_u = grad_u_;
	gu.O = antisymmetric(grad_u_);
	gu.E = symmetric(grad_u_);

    // we define the shape of the flow as E and O normalized by the largest (in amplitude) eigenvalue of E.
	auto evals = [](double det)->std::array<double,3> // eigenvalues of a 3-by-3 real symmetric traceless matrix
					{auto theta = acos(0.5*det)/3.;
					 return {2*cos(theta), 2*cos(theta+2*M_PI/3), 2*cos(theta+4*M_PI/3)};};
	auto eigvals = evals(det(gu.E));
	std::sort(eigvals.begin(), eigvals.end(), [](double a, double b){return abs(b) < abs(a);});
	auto deformation_rate = abs(eigvals[0]);
    gu.O_shape = gu.O/deformation_rate;
    gu.E_shape = gu.E/deformation_rate;

	return gu;
};
}
#endif /* defined(__LF_DEM__GradU__) */
