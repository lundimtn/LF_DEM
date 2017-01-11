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
#include "Matrix.h"
#include "vec3d.h"
#include <iostream>
#include <iomanip>
#include <array>

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
    for (unsigned int i=0; i<6; i++)
      elm[i] = a;
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

  inline void set(const matrix &m)
  {
    elm[0] = m.elm[0]; // xx
    elm[1] = m.elm[1]; // xy
    elm[2] = m.elm[2]; // xz
    elm[3] = m.elm[5]; // yz
    elm[4] = m.elm[4]; // yy
    elm[5] = m.elm[8]; // zz
  }

  void reset() {
    for (int i=0; i<6; i++) {
      elm[i] = 0;
    }
  }

  /*  N1 = Sxx-Szz;
   */
  double getNormalStress1()
  {
    return elm[0]-elm[5];
  }

  /*  N2 = Szz-Syy;
   */
  double getNormalStress2()
  {
    return elm[5]-elm[4];
  }

  double getParticlePressure()
  {
    return -(1./3)*(elm[0]+elm[4]+elm[5]);
  }

};

// Helper functions
inline double shearStressComponent(const Sym2Tensor& s, double theta_shear)
{
  return cos(theta_shear)*s.elm[2]+sin(theta_shear)*s.elm[3];
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
  out << st.elm[0] << ' ';
  out << st.elm[1] << ' ';
  out << st.elm[2] << ' ';
  out << st.elm[3] << ' ';
  out << st.elm[4] << ' ';
  out << st.elm[5] << ' ';
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
friend Sym2Tensor operator / (const Sym2Tensor& s,
                              const T& a)
{
  return {s.elm[0]/a,
          s.elm[1]/a,
          s.elm[2]/a,
          s.elm[3]/a,
          s.elm[4]/a,
          s.elm[5]/a};
}

#endif
