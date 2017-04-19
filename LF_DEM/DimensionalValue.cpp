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

void UnitSystem::add(Unit::Unit unit, DimensionalValue<double> value)
{
  assert(value.dimension == Force || value.dimension == Stress);
  if (value.value == 0) {
    return;
  }
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

Unit::Unit UnitSystem::getLargestUnit()
{
  auto largest_unit = unit_nodes.cbegin()->first;
  double largest_value = 1;

  setInternalUnit(largest_unit);
  for (auto &node: unit_nodes) {
    if (node.second.value > largest_value) {
      largest_value = node.second.value;
      largest_unit = node.first;
    }
  }
  return largest_unit;
}

void UnitSystem::convertToParentUnit(DimensionalValue<double> &node)
{
  auto &parent_node = unit_nodes[node.unit];
  node.value *= parent_node.value;
  node.unit = parent_node.unit;
}

void UnitSystem::convertNodeUnit(DimensionalValue<double> &node, Unit::Unit unit)
{
  if (node.unit != unit) {
    auto &parent_node = unit_nodes[node.unit];
    if (parent_node.unit != node.unit) {
      convertNodeUnit(parent_node, unit);
      convertToParentUnit(node);
    } else { // parent_node is a root, but is not unit: the unit system is not closed
      throw std::runtime_error(" UnitSystem:: cannot express "
                               +Unit::unit2suffix(node.unit)+" in "+Unit::unit2suffix(unit)+" units ");
    }
  }
}

void UnitSystem::flipDependency(Unit::Unit node_name)
{
  auto &node = unit_nodes[node_name];
  const auto &parent_node_name =  node.unit;
  auto &parent_node = unit_nodes[parent_node_name];
  if (node_name==parent_node_name) {
    return;
  }
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
  if (unit_nodes.at(unit).unit != unit) { // if the unit force is expressed in other units than itself
    flipDependency(unit);
  }
  unit_nodes[unit] = {Force, 1, unit};

  for (auto &node: unit_nodes) {
    convertNodeUnit(node.second, unit);
  }
}

} // namespace Dimensional
