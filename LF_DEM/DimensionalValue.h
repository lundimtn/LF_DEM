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
  template<typename T> void convertToInternalUnit(DimensionalValue<T> &value) const;
  template<typename T> void convertFromInternalUnit(DimensionalValue<T> &value, Unit::Unit unit) const;
  const std::map<Unit::Unit, DimensionalValue<double>> getForceTree() const {return unit_nodes;};
  Unit::Unit getLargestUnit() const;
private:
  std::map<Unit::Unit, DimensionalValue<double>> unit_nodes;
  void convertToParentUnit(DimensionalValue<double> &node);
  void flipDependency(Unit::Unit node_name);
  void convertNodeUnit(DimensionalValue<double> &node, Unit::Unit unit);
  template<typename T> void convertUnit(DimensionalValue<T> &value,
                                        Unit::Unit new_unit) const; // for arbitrary Dimension
};

template<typename T>
void UnitSystem::convertUnit(DimensionalValue<T> &quantity, Unit::Unit new_unit) const
{
  /**
    \brief Convert quantity to new_unit.
    */

  // special cases
  if (quantity.dimension == none || quantity.dimension == Strain) {
    return;
  }
  if (quantity.dimension == TimeOrStrain) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert units of a TimeOrStrain quantity.");
  }

  if (!unit_nodes.count(quantity.unit)) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert from "+unit2suffix(quantity.unit)+" unit.");
  }
  if (!unit_nodes.count(new_unit)) {
    throw std::runtime_error("UnitSystem::convertUnits : cannot convert to "+unit2suffix(new_unit)+" unit.");
  }

  // value of the old force unit in new_unit
  double force_unit_ratio = unit_nodes.at(quantity.unit).value/unit_nodes.at(new_unit).value;
  switch (quantity.dimension) {
    case Force:
      quantity.value *= force_unit_ratio;
      break;
    case Time:
      quantity.value /= force_unit_ratio;
      break;
    case Rate:
      quantity.value *= force_unit_ratio;
      break;
    case Viscosity:
      quantity.value *= 6*M_PI;  // viscosity in solvent viscosity units irrespective the unit system chosen
      break;
    case Stress:
      quantity.value *= 6*M_PI*force_unit_ratio;
      break;
    case Velocity:
      quantity.value *= force_unit_ratio;
      break;
    case none: case Strain: case TimeOrStrain:
      // Keep these cases here even if unreachable code (see above if statements)
      // in order to avoid gcc -Wswitch warnings.
      // While default: case would achieve a similar results,
      // it would also switch off ALL -Wswitch warnings, which is dangerous.
      break;
  }
  quantity.unit = new_unit;
}

template<typename T>
void UnitSystem::convertToInternalUnit(DimensionalValue<T> &quantity) const
{
  auto internal_unit = unit_nodes.cbegin()->second.unit;
  convertUnit(quantity, internal_unit);
}

template<typename T>
void UnitSystem::convertFromInternalUnit(DimensionalValue<T> &quantity, Unit::Unit unit) const
{
  auto internal_unit = unit_nodes.cbegin()->second.unit;
  assert(quantity.unit == internal_unit);

  convertUnit(quantity, unit);
}

} // namespace Dimensional
#endif // #ifndef __LF_DEM__Dimensional__
