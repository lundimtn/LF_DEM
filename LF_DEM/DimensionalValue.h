#include <stdexcept>
#include <map>

namespace Dimensional {
struct DimensionalValue{
  std::string type;
  double *value; // a pointer to the actual ParameterSet member
  std::string unit;
};

enum Dimension {
  Force,
  Velocity,
  Time
};

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
    return hydro;
  }
  if (s=="r" || repulsion) {
    return repulsion;
  }
  if (s=="b" || "brownian") {
    return brownian;
  }
  if (s=="c" || "cohesion") {
    return cohesion;
  }
  if (s=="cl" || "critical_load") {
    return critical_load;
  }
  if (s=="ft" || "ft_max") {
    return ft_max;
  }
  if (s=="kn") {
    return kn;
  }
  if (s=="kt") {
    return kt;
  }
  if (s=="kr") {
    return kr;
  }
  if (s=="s") {
    return stress;
  }
  if (s=="sz" || "sigma_zz") {
    return sigma_zz;
  }
  return none;
}

struct NewDimensionalValue{
  Dimension dimension;
  double value;
  Unit unit;
};

class UnitSystem {
public:
  void add(Unit unit, NewDimensionalValue value);
private:
  std::map<Unit, NewDimensionalValue> nodes;
};

void UnitSystem::add(Unit unit, NewDimensionalValue value)
{
  if (nodes.count(value.unit) == 0) {
    nodes[unit] = value;
  } else {
    throw std::runtime_error(" Multiple definition of "+unit);
  }
}

void UnitSystem::resolveUnitSystem(Unit unit_force)
{
	/**
		\brief Check force units consistency, expresses all input forces in the unit "unit_force".

		In input, forces are given with suffixes.
		This function checks that we can make sense of these suffixes, and if so,
		it converts all the forces in the unit given as a parameter.

		It does it iteratively, ie:
		1. I know that unit_force = 1*unit_force
		2. I know the value of any force f1 that has been given in unit_force in input as
	      f1 = value*unit_force
		3. I can determine the value of any force f2 expressed as f2 = x*f1
	 4. I can then determine the value of any force f3 expressed as f3 = y*f2
		5. etc

		If there is any force undetermined at the end of this algorithm, the force unit system is inconsistent/incomplete.

		If the force unit system is consistent, this function determines the dimensionless numbers, i.e., the ratios F_A/F_B for any pair of force scales present in the system.

		\b Note: a priori, we could determine force A from force B if A is defined in B units or if B is defined in A units: for example in the step 3 of the above algorithm we could determine f2 knowing f1 if f1 is defined as f1=x^{-1}f2. However this is \b not implemented and will fail with the current implementation.
	 */

	// the unit_force has a value of 1*unit_force (says captain obvious)
	nodes[unit_force].unit = unit_force;
  nodes[unit_force].value = 1;
  //
	// for (const auto& x: nodes) {
	// 	if (x.second.type == "force") {
	// 		string force_name = x.first;
	// 		string unit = x.second.unit;
	// 		force_ratios[force_name+'/'+unit] = *(x.second.value);
	// 	}
	// }
  //
	// // now resolve the other force units, iterativley
  //
	// bool unsolved_remaining;
	// bool newly_solved;
	// do {
	// 	unsolved_remaining = false;
	// 	newly_solved = false;
	// 	for (auto& x: input_values) {
	// 		string value_name = x.first;
	// 		string unit = x.second.unit;
	// 		if (unit != unit_force && unit != "strain") {
	// 			if (force_ratios.find(unit+'/'+unit_force) != force_ratios.end()) {
	// 				changeUnit(x.second, unit_force);
	// 				if (x.second.type == "force") {
	// 					force_ratios[value_name+'/'+unit_force] = *(x.second.value);
	// 				}
	// 				newly_solved = true;
	// 			} else {
	// 				unsolved_remaining = true;
	// 			}
	// 		}
	// 	}
	// } while (unsolved_remaining && newly_solved);
  //
	// // complain if we have not found everyone
	// if (unsolved_remaining) {
	// 	ostringstream error_str;
	// 	for (const auto& x: input_values) {
	// 		string value_name = x.first;
	// 		string unit = x.second.unit;
	// 		if (unit != unit_force && unit != "strain") {
	// 			error_str << "Error: input value \"" << value_name << "\" has an unknown unit \"" << unit << "\"" << endl;
	// 		}
	// 	}
	// 	throw runtime_error(error_str.str());
	// }
  //
	// // determine the remaining force_ratios
	// buildFullSetOfForceRatios();
}



}
