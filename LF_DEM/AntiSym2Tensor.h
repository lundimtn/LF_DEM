//
//  AntiSym2Tensor.h
//  LF_DEM
//
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class AntiSym2Tensor
 \brief Anti-Symmetric rank 2 tensor object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_AntiSym2Tensor_h
#define LF_DEM_AntiSym2Tensor_h
#include "vec3d.h"
#include <iostream>
#include <iomanip>
#include "TensorBase.h"
#include "Matrix.h"


class matrix;

namespace Algebra
{
class AntiSym2Tensor : public TensorBase<3, AntiSym2Tensor> {// elm = (12, 13, 23)
public:
	void transpose() {
		for (auto &e: elm) {
			e = -e;
		}
	}
};

// Helper functions
inline vec3d dot(const AntiSym2Tensor& s,
                 const vec3d& v)
{
	return {s[0]*v.y + s[1]*v.z,
	        -s[0]*v.x + s[2]*v.z,
	        -s[1]*v.x - s[2]*v.y};
}

inline vec3d dot(const vec3d& v,
                 const AntiSym2Tensor& s)
{
	return -dot(s, v); // s anti-symmetric
}

inline AntiSym2Tensor antisymmetric(const matrix &M) {
	AntiSym2Tensor a;
	a.elm = {{0.5*(M.elm[1]-M.elm[3]),
	          0.5*(M.elm[2]-M.elm[6]),
	          0.5*(M.elm[5]-M.elm[7])}};
	return a;
}

inline AntiSym2Tensor transpose(const AntiSym2Tensor &s) {
	return -s;
}

}
#endif // #ifndef LF_DEM_AntiSym2Tensor_h
