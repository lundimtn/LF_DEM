//
//  DimensionalValue.cpp
//  LF_DEM
//
//  Copyright (c) 2017 Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include <cmath> // for M_PI
#include "DimensionalValue.h"

namespace Dimensional {

void UnitSystem::add(Unit::Unit unit, NewDimensionalValue value)
{
  assert(value.dimension == Force || value.dimension == Stress);
  if (value.dimension==Stress) {
    value.value /= 6*M_PI; // at some point we have to get rid of this weird unit choice
  }
  unit_nodes[unit] = value;
  auto parent_node_name = value.unit;
  // no orphans!
  if (unit_nodes.count(parent_node_name) == 0) {
    unit_nodes[parent_node_name] = {Force, 1, parent_node_name};
  }
}

void UnitSystem::convertToParentUnit(NewDimensionalValue &node)
{
  auto &parent_node = unit_nodes[node.unit];
  node.value *= parent_node.value;
  node.unit = parent_node.unit;
}

void UnitSystem::convertUnit(NewDimensionalValue &node, Unit::Unit unit)
{
  if (node.unit != unit) {
    auto &parent_node = unit_nodes[node.unit];
    if (parent_node.unit != node.unit) {
      convertUnit(parent_node, unit);
      convertToParentUnit(node);
    } else { // parent_node is a root, but is not unit: the unit system is not closed
      throw std::runtime_error(" UnitSystem:: unable to solve dependency problem. You did not provide a consistent set of dimensionless numbers.");
    }
  }
}

void UnitSystem::flipDependency(Unit::Unit node_name)
{
  auto &node = unit_nodes[node_name];
  const auto &parent_node_name =  node.unit;
  auto &parent_node = unit_nodes[parent_node_name];

  flipDependency(parent_node_name);
  parent_node.value = 1/node.value;
  parent_node.unit = node_name;
}

void UnitSystem::setInternalUnit(Unit::Unit unit)
{
	/**
		\brief Check force units consistency, expresses all input forces in the unit "unit".
	 */

	// the unit has a value of 1*unit (says captain obvious)
  if (unit_nodes[unit].unit != unit) { // if the unit force is expressed in other units than itself
    flipDependency(unit);
  }
  unit_nodes[unit] = {Force, 1, unit};

  for (auto &node: unit_nodes) {
    convertUnit(node.second, unit);
  }
}

void UnitSystem::convertUnits(NewDimensionalValue &value, const NewDimensionalValue &new_unit)
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
  }
}

void UnitSystem::convertToInternalUnit(NewDimensionalValue &value)
{
  if (unit_nodes.count(value.unit) == 0) {
    throw std::runtime_error(" UnitSystem::expressInUnit : unknown unit for one of the parameters.");
  }
  auto &unit_node = unit_nodes[value.unit];
  convertUnits(value, unit_node);
}

void UnitSystem::convertFromInternalUnit(NewDimensionalValue &value, Unit::Unit unit)
{
  if (unit_nodes.count(value.unit) == 0) {
    throw std::runtime_error(" UnitSystem::expressInUnit : unknown unit for one of the parameters.");
  }
  NewDimensionalValue internal_force = {Force, 1/unit_nodes[unit].value, unit};
  convertUnits(value, internal_force);
}

} // namespace Dimensional
