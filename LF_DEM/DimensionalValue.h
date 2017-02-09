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
  Time
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
  if (s=="r" || "repulsion") {
    return Unit::repulsion;
  }
  if (s=="b" || "brownian") {
    return Unit::brownian;
  }
  if (s=="c" || "cohesion") {
    return Unit::cohesion;
  }
  if (s=="cl" || "critical_load") {
    return Unit::critical_load;
  }
  if (s=="ft" || "ft_max") {
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
  if (s=="sz" || "sigma_zz") {
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
  void add(Param::Parameter param, NewDimensionalValue value);
  void add(Unit::Unit unit, NewDimensionalValue value);
private:
  std::map<Param::Parameter, NewDimensionalValue> param_nodes;
  std::map<Unit::Unit, NewDimensionalValue> unit_nodes;
  void convertToParentUnit(NewDimensionalValue&);
  void convertUnit(NewDimensionalValue&, Unit::Unit);
  void resolveUnitSystem(Unit::Unit);
};



} // namespace Dimensional
