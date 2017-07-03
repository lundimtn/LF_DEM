//
//  Sym2Tensor.h
//  LF_DEM
//
//  Created by Ryohei Seto on 4/16/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

/**
 \class Sym2Tensor
 \brief Stress tensor object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_Sym2Tensor_h
#define LF_DEM_Sym2Tensor_h
#include "vec3d.h"
#include <iostream>
#include <iomanip>
#include <array>
#include "Matrix.h"

class matrix;

class Sym2Tensor {
public:
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 */
	double elm [6];

	Sym2Tensor()
	{
		reset();
	}

	Sym2Tensor(double a)
	{
		for (unsigned int i=0; i<6; i++) {
			elm[i] = a;
		}
	}

	Sym2Tensor(std::initializer_list<double> il)
	{
		std::copy(il.begin(), il.end(), elm);
	}

	Sym2Tensor& operator += (const Sym2Tensor& s)
	{
		for (int i=0; i<6; i++) {
			elm[i] += s.elm[i];
		}
		return *this;
	}

	Sym2Tensor& operator -= (const Sym2Tensor& s)
	{
		for (int i=0; i<6; i++) {
			elm[i] -= s.elm[i];
		}
		return *this;
	}

	template <typename T>
	Sym2Tensor& operator *= (const T& d)
	{
		for (int i=0; i<6; i++) {
			elm[i] *= d;
		}
		return *this;
	}

	template <typename T>
	Sym2Tensor& operator /= (const T& d)
	{
		double d_inv = 1.0/d;
		for (int i=0; i<6; i++) {
			elm[i] *= d_inv;
		}
		return *this;
	}

	void reset() {
		for (int i=0; i<6; i++) {
			elm[i] = 0;
		}
	}

	vec3d diag()
	{
		return {elm[0], elm[4], elm[5]};
	}

	double trace()
	{
		return elm[0]+elm[4]+elm[5];
	}

	double selfdoubledot()
	{
		return elm[0]*elm[0]
		       + 2*elm[1]*elm[1]
		       + 2*elm[2]*elm[2]
		       + 2*elm[3]*elm[3]
		       + elm[4]*elm[4]
		       + elm[5]*elm[5];
	}

	matrix getMatrix()
	{
		matrix m(elm[0], elm[1], elm[2],
		         elm[1], elm[4], elm[3],
		         elm[2], elm[3], elm[5]);
		return m;
	}
};

inline Sym2Tensor operator + (const Sym2Tensor& s)
{
	return s;
}

inline Sym2Tensor operator + (const Sym2Tensor& a1,
                              const Sym2Tensor& a2)
{
	return {a1.elm[0]+a2.elm[0],
	        a1.elm[1]+a2.elm[1],
	        a1.elm[2]+a2.elm[2],
	        a1.elm[3]+a2.elm[3],
	        a1.elm[4]+a2.elm[4],
	        a1.elm[5]+a2.elm[5]};
}

/* subtraction */
inline Sym2Tensor operator - (const Sym2Tensor& s)
{
	return {-s.elm[0],
	        -s.elm[1],
	        -s.elm[2],
	        -s.elm[3],
	        -s.elm[4],
	        -s.elm[5]};
}

inline Sym2Tensor operator - (const Sym2Tensor& a1,
                              const Sym2Tensor& a2)
{
	return {a1.elm[0]-a2.elm[0],
	        a1.elm[1]-a2.elm[1],
	        a1.elm[2]-a2.elm[2],
	        a1.elm[3]-a2.elm[3],
	        a1.elm[4]-a2.elm[4],
	        a1.elm[5]-a2.elm[5]};
}

// output stream operator
inline std::ostream& operator << (std::ostream& out,
                                  const Sym2Tensor& st)
{
	out << st.elm[0] << ' '\
	    << st.elm[1] << ' '\
	    << st.elm[2] << ' '\
	    << st.elm[3] << ' '\
	    << st.elm[4] << ' '\
	    << st.elm[5] << ' ';
	return out;
}

template <typename T>
inline Sym2Tensor operator * (const T& a,
                              const Sym2Tensor& s)
{
	return {a*s.elm[0],
	        a*s.elm[1],
	        a*s.elm[2],
	        a*s.elm[3],
	        a*s.elm[4],
	        a*s.elm[5]};
}

template <typename T>
inline  Sym2Tensor operator * (const Sym2Tensor& s,
                               const T& a)
{
	return a*s;
}

template <typename T>
inline Sym2Tensor operator / (const Sym2Tensor& s,
                              const T& a)
{
	return {s.elm[0]/a,
	        s.elm[1]/a,
	        s.elm[2]/a,
	        s.elm[3]/a,
	        s.elm[4]/a,
	        s.elm[5]/a};
}

// Helper functions

inline vec3d dot(const Sym2Tensor& s,
                 const vec3d& v)
{
	//  s = (xx, xy, xz, yz, yy, zz)
	return {s.elm[0]*v.x + s.elm[1]*v.y + s.elm[2]*v.z,
	        s.elm[1]*v.x + s.elm[4]*v.y + s.elm[3]*v.z,
	        s.elm[2]*v.x + s.elm[3]*v.y + s.elm[5]*v.z};
}

inline vec3d dot(const vec3d& v,
                 const Sym2Tensor& s)
{
	return dot(s, v);   // s is symmetric
}

inline double doubledot(const Sym2Tensor& s1, const Sym2Tensor& s2)
{
	return s1.elm[0]*s2.elm[0]
	       + 2*s1.elm[1]*s2.elm[1]
	       + 2*s1.elm[2]*s2.elm[2]
	       + 2*s1.elm[3]*s2.elm[3]
	       + s1.elm[4]*s2.elm[4]
	       + s1.elm[5]*s2.elm[5];
}

inline Sym2Tensor outer_sym(const vec3d& v1, const vec3d& v2)
{
	return {v1.x*v2.x, //xx
	        0.5*(v1.x*v2.y+v1.y*v2.x), //xy
	        0.5*(v1.x*v2.z+v1.z*v2.x), //xz
	        0.5*(v1.y*v2.z+v1.z*v2.y), //yz
	        v1.y*v2.y, //yy
	        v1.z*v2.z}; //zz
}

inline Sym2Tensor outer(const vec3d& v)
{
	return {v.x*v.x,
	        v.x*v.y,
	        v.x*v.z,
	        v.y*v.z,
	        v.y*v.y,
	        v.z*v.z};
}

namespace Algebra
{
	inline Sym2Tensor symmetric(const matrix &M) {
		return {M.elm[0], 0.5*(M.elm[1]+M.elm[3]), 0.5*(M.elm[2]+M.elm[6]),
		                  M.elm[4],                0.5*(M.elm[5]+M.elm[7]),
     	                                           M.elm[8]};
	}
	inline double det(const Sym2Tensor &T) {
		return T.elm[0]*(T.elm[4]*T.elm[5]-T.elm[3]*T.elm[3])
				- T.elm[1]*(T.elm[1]*T.elm[5]-T.elm[3]*T.elm[2])
				+ T.elm[2]*(T.elm[1]*T.elm[3]-T.elm[4]*T.elm[2]);
	}
}

#endif
