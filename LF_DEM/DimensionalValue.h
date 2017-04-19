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
#include <iostream>

namespace Dimensional {

enum Dimension {
  Force,
  Time,
  Viscosity,
  Stress,
  Rate,
  Velocity,
  TimeOrStrain,
  Strain,
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

inline Unit getUnit(std::string s) {
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

inline std::string unit2suffix(Unit unit) {
  if (unit==hydro) {
    return "h";
  }
  if (unit==repulsion) {
    return "r";
  }
  if (unit==brownian) {
    return "b";
  }
  if (unit==cohesion) {
    return "c";
  }
  if (unit==critical_load) {
    return "cl";
  }
  if (unit==ft_max) {
    return "ft";
  }
  if (unit==kn) {
    return "kn";
  }
  if (unit==kt) {
    return "kt";
  }
  if (unit==kr) {
    return "kr";
  }
  if (unit==stress) {
    return "s";
  }
  if (unit==sigma_zz) {
    return "sz";
  }
  if (unit==none) {
    return "";
  }
  return "";
}
} // namespace Unit

template<typename T>
struct DimensionalValue {
  Dimension dimension;
  T value;
  Unit::Unit unit;
};


inline bool getSuffix(const std::string& str, std::string& value, std::string& suffix)
{
	size_t suffix_pos = str.find_first_of("abcdfghijklmnopqrstuvwxyz"); // omission of "e" is intended, to allow for scientific notation like "1e5h"
	value = str.substr(0, suffix_pos);
	if (suffix_pos != str.npos) {
		suffix = str.substr(suffix_pos, str.length());
		return true;
	} else {
		return false;
	}
}

inline void errorNoSuffix(std::string quantity)
{
	std::cerr << "Error : no unit scale (suffix) provided for " << quantity << std::endl; exit(1);
}


inline DimensionalValue<double> str2DimensionalValue(Dimension dimension,
                                                     std::string value_str,
                                                     std::string name)
{
	DimensionalValue<double> inv;
	inv.dimension = dimension;

	std::string numeral, suffix;
	bool caught_suffix = true;
	caught_suffix = getSuffix(value_str, numeral, suffix);
	if (!caught_suffix) {
		errorNoSuffix(name);
	}
	inv.value = stod(numeral);
	inv.unit = Unit::getUnit(suffix);

  if (inv.dimension == TimeOrStrain) {
    if (inv.unit == Unit::hydro) {
      inv.dimension = Strain;
      inv.unit = Unit::none;
    } else {
      inv.dimension = Time;
    }
  }
	return inv;
}

class UnitSystem {
public:
  // void add(Param::Parameter param, DimensionalValue value);
  void add(Unit::Unit unit, DimensionalValue<double> value);
  void setInternalUnit(Unit::Unit unit);
  template<typename T> void convertToInternalUnit(DimensionalValue<T> &value);
  template<typename T> void convertFromInternalUnit(DimensionalValue<T> &value, Unit::Unit unit);
  const std::map<Unit::Unit, DimensionalValue<double>> getForceTree() {return unit_nodes;};
private:
  std::map<Unit::Unit, DimensionalValue<double>> unit_nodes;
  void convertToParentUnit(DimensionalValue<double> &node);
  void flipDependency(Unit::Unit node_name);
  void convertNodeUnit(DimensionalValue<double> &node, Unit::Unit unit);
  template<typename T> void convertUnits(DimensionalValue<T> &value,
                                         const DimensionalValue<double> &new_unit); // for arbitrary Dimension
};

template<typename T>
void UnitSystem::convertUnits(DimensionalValue<T> &value, const DimensionalValue<double> &new_unit)
{
  switch (value.dimension) {
    case none:
      break;
    case Force:
      value.value *= new_unit.value;
      break;
    case Time:
      value.value /= new_unit.value;
      break;
    case Rate:
      value.value *= new_unit.value;
      break;
    case Viscosity:
      value.value *= 6*M_PI;  // viscosity in solvent viscosity units irrespective the unit system chosen
      break;
    case Stress:
      value.value *= 6*M_PI*new_unit.value;
      break;
    case Velocity:
      value.value *= new_unit.value;
      break;
    case Strain:
      break;
    case TimeOrStrain:
      throw std::runtime_error("UnitSystem::convertUnits : cannot convert units of a TimeOrStrain quantity.");
      break;
  }
}

template<typename T>
void UnitSystem::convertToInternalUnit(DimensionalValue<T> &value)
{
  if (value.unit != Unit::none && unit_nodes.count(value.unit) == 0) {
    throw std::runtime_error(" UnitSystem::convertToInternalUnit : unknown unit "+Unit::unit2suffix(value.unit));
  }
  auto &unit_node = unit_nodes[value.unit];
  convertUnits(value, unit_node);
}

template<typename T>
void UnitSystem::convertFromInternalUnit(DimensionalValue<T> &value, Unit::Unit unit)
{
  if (value.unit != Unit::none && unit_nodes.count(value.unit) == 0) {
    throw std::runtime_error(" UnitSystem::convertFromInternalUnit : unknown unit "+Unit::unit2suffix(value.unit));
  }
  DimensionalValue<double> internal_force = {Force, 1/unit_nodes[unit].value, unit};
  convertUnits(value, internal_force);
}

} // namespace Dimensional
#endif // #ifndef __LF_DEM__Dimensional__
