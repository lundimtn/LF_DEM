//
//  DimensionalValue.h
//  LF_DEM
//
//  Copyright (c) 2017 Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Dimensional__
#define __LF_DEM__Dimensional__
#include <stdexcept>
#include <map>
#include <assert.h>

namespace Dimensional {
struct DimensionalValue{
  std::string type;
  double *value; // a pointer to the actual Param::ParameterSet member
  std::string unit;
};

enum Dimension {
  Force,
  Time,
  Viscosity,
  Stress,
  Rate,
  Velocity,
  none
};

namespace Unit {
enum Unit {
  hydro,
  repulsion,
  brownian,
  cohesion,
  critical_load,
  ft_max,
  kn,
  kt,
  kr,
  stress,
  sigma_zz,
  none
};
} // namespace Unit

namespace Param {
enum Parameter {
  hydro,
  repulsion,
  brownian,
  cohesion,
  critical_load,
  ft_max,
  kn,
  kt,
  kr,
  stress,
  sigma_zz
};
}// namespace Param

inline Unit::Unit getUnit(std::string s) {
  if (s=="h") {
    return Unit::hydro;
  }
  if (s=="r" || s=="repulsion") {
    return Unit::repulsion;
  }
  if (s=="b" || s=="brownian") {
    return Unit::brownian;
  }
  if (s=="c" || s=="cohesion") {
    return Unit::cohesion;
  }
  if (s=="cl" || s=="critical_load") {
    return Unit::critical_load;
  }
  if (s=="ft" || s=="ft_max") {
    return Unit::ft_max;
  }
  if (s=="kn") {
    return Unit::kn;
  }
  if (s=="kt") {
    return Unit::kt;
  }
  if (s=="kr") {
    return Unit::kr;
  }
  if (s=="s") {
    return Unit::stress;
  }
  if (s=="sz" || s=="sigma_zz") {
    return Unit::sigma_zz;
  }
  return Unit::none;
}

struct NewDimensionalValue{
  Dimension dimension;
  double value;
  Unit::Unit unit;
};

class UnitSystem {
public:
  // void add(Param::Parameter param, NewDimensionalValue value);
  void add(Unit::Unit unit, NewDimensionalValue value);
  void setInternalUnit(Unit::Unit unit);
  void expressInUnit(NewDimensionalValue &value);
  void convertToInternalUnit(NewDimensionalValue &value);
  void convertFromInternalUnit(NewDimensionalValue &value, Unit::Unit unit);
private:
  std::map<Unit::Unit, NewDimensionalValue> unit_nodes;
  void convertToParentUnit(NewDimensionalValue &node);
  void flipDependency(Unit::Unit node_name);
  void convertUnit(NewDimensionalValue& node, Unit::Unit unit);
  void convertUnits(NewDimensionalValue &value, const NewDimensionalValue &new_unit); // for arbitrary Dimension
};



} // namespace Dimensional
#endif // #ifndef __LF_DEM__Dimensional__
