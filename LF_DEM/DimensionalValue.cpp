#include "DimensionalValue.h"

namespace Dimensional {

void UnitSystem::add(Param::Parameter param, NewDimensionalValue value)
{
  if (param_nodes.count(param) == 0) {
    param_nodes[param] = value;
  } else {
    throw std::runtime_error(" Multiple definition of "+param);
  }
}

void UnitSystem::add(Unit::Unit unit, NewDimensionalValue value)
{
  if (unit_nodes.count(unit) == 0) {
    unit_nodes[unit] = value;
  } else {
    throw std::runtime_error(" Multiple definition of "+unit);
  }
}

void UnitSystem::convertToParentUnit(NewDimensionalValue &value)
{
  auto &unit_value = unit_nodes[value.unit];
  switch (value.dimension) {
    case Force:
      switch (unit_value.dimension) {
        case Force:
          value.value *= unit_value.value;
          break;
        case Time:
          value.value /= unit_value.value;
          break;
      }
      break;
    case Time:
      switch (unit_value.dimension) {
        case Force:
          value.value /= unit_value.value;
          break;
        case Time:
          value.value *= unit_value.value;
          break;
      }
      break;
  }
}

void UnitSystem::convertUnit(NewDimensionalValue &value, Unit::Unit unit_force)
{
  if (value.unit != unit_force) {
    auto &unit_value = unit_nodes[value.unit];
    convertUnit(unit_value, unit_force);
    convertToParentUnit(value);
  }
}

void UnitSystem::resolveUnitSystem(Unit::Unit unit_force)
{
	/**
		\brief Check force units consistency, expresses all input forces in the unit "unit_force".
	 */

	// the unit_force has a value of 1*unit_force (says captain obvious)
	unit_nodes[unit_force].unit = unit_force;
  unit_nodes[unit_force].value = 1;

  for (auto &node: unit_nodes) {
    convertUnit(node.second, unit_force);
  }
  for (auto &node: param_nodes) {
    convertUnit(node.second, unit_force);
  }
}

} // namespace Dimensional
