#ifndef __LF_DEM__GradU__
#define __LF_DEM__GradU__
#include "Matrix.h"
#include "AntiSym2Tensor.h"
#include "Sym2Tensor.h"

namespace Algebra
{
struct GradU {
    matrix full;
    Sym2Tensor sym;
    AntiSym2Tensor antisym;
};

struct GradU buildGradU(matrix grad_u_) {
    struct GradU gu;
    gu.full = grad_u_;
    gu.antisym = antisymmetric(grad_u_);
    gu.sym = symmetric(grad_u_);
    return gu;
};
}
#endif /* defined(__LF_DEM__GradU__) */
