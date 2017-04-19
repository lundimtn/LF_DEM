//
//  SimulationInit.cpp
//  LF_DEM
//
//  Created by Romain Mari on 08/10/15.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Simulation.h"
#include <string>
#include <sstream>
#include <stdexcept>
#include "Configuration.h"


using namespace std;

void Simulation::importPreSimulationData(string filename)
{
	// @@@ DEPRECATED?
	ifstream fin_PreSimulationData;
	fin_PreSimulationData.open(filename.c_str());
	if (!fin_PreSimulationData) {
		ostringstream error_str;
		error_str  << " Pre-simulation data file '" << filename << "' not found." << endl;
		throw runtime_error(error_str.str());
	}
	double stress_, shear_rate_;
	while (fin_PreSimulationData >> stress_ >> shear_rate_) {
		if (stress_ == target_stress_input) {
			break;
		}
	}
	shear_rate_expectation = shear_rate_;
}

void Simulation::echoInputFiles(string in_args,
                                vector<string>& input_files)
{
	/**
	 \brief Print the entire information needed to reproduce the simulation in Simulation::fout_input
	 */
	fout_input << "# LF_DEM version " << GIT_VERSION << ", called with:" << endl;
	fout_input << in_args << endl << endl;
	for (const string& in_file : input_files) {
		ifstream in_f;
		string line;
		in_f.open(in_file.c_str());
		if (in_f.is_open()) {
			fout_input << "********** File " << in_file << " ************" << endl << endl;
			while (in_f.good()) {
				getline(in_f, line);
				fout_input << line << endl;
			}
			fout_input << endl << endl;
		}
		in_f.close();
	}
	fout_input.close();
}

void Simulation::setLowPeclet()
{
	sys.lowPeclet = true;
	double scale_factor_SmallPe = p.Pe_switch/force_ratios["hydro/thermal"];
	p.dt *= p.Pe_switch; // to make things continuous at Pe_switch
}

Dimensional::Unit::Unit Simulation::pickInternalUnitsRateControl()
{
	/**
	 \brief Determine the best internal force scale to run the simulation (rate controlled case).

		If the system is non-Brownian, the hydrodynamic force unit is taken (\b note: this will change in the future). If the system is Brownian, the Brownian force unit is selected at low Peclet (i.e., Peclet numbers smaller that ParameterSet::Pe_switch) and the hydrodynamic force unit is selected at high Peclet.
	 */
	units.setInternalUnit(Dimensional::Unit::hydro);
	auto unit_tree = units.getForceTree();
	bool is_brownian = unit_tree.count(Dimensional::Unit::brownian) > 0;
	if (is_brownian) {
		double inverse_Peclet = unit_tree[Dimensional::Unit::brownian].value;
		if (inverse_Peclet < 1/p.Pe_switch && !sys.zero_shear) {
			return Dimensional::Unit::hydro;
		} else { // low Peclet mode
			setLowPeclet();
			return Dimensional::Unit::brownian;
		}
	} else {
		return Dimensional::Unit::hydro;
	}
}

void Simulation::exportForceAmplitudes()
{
	/**
	 \brief Copy the input force alues in the ForceAmplitude struct of the System class
	 */
	string indent = "  Simulation::\t";
	cout << indent+"Forces used:" << endl;
	indent += "\t";

	using namespace Dimensional::Unit;
	auto forces = units.getForceTree();
	sys.repulsiveforce = forces.count(repulsion) > 0;
	if (sys.repulsiveforce) {
		sys.p.repulsion = forces[repulsion].value;
		cout << indent+"Repulsive force (in \"" << Dimensional::Unit::unit2suffix(forces[repulsion].unit) << "\" units): " << sys.p.repulsion << endl;
	}

	sys.critical_load = forces.count(critical_load) > 0;
	if (sys.critical_load) {
		sys.p.critical_load = forces[critical_load].value;
		cout << indent+"Critical Load (in \"" << Dimensional::Unit::unit2suffix(forces[critical_load].unit) << "\" units): " << sys.p.critical_load << endl;
	}

	sys.cohesion = forces.count(cohesion) > 0;
	if (sys.cohesion) {
		sys.p.cohesion = forces[cohesion].value;
		cout << indent+"Cohesion (in \"" << Dimensional::Unit::unit2suffix(forces[cohesion].unit) << "\" units): " << sys.p.cohesion << endl;
	}

	bool is_ft_max = forces.count(ft_max) > 0;
	if (is_ft_max) {
		sys.p.ft_max = forces[ft_max].value;
		cout << indent+"Max Tangential Load (in \"" << Dimensional::Unit::unit2suffix(forces[ft_max].unit) << "\" units): " << sys.p.ft_max << endl;
	}

	sys.brownian = forces.count(brownian) > 0;
	if (sys.brownian) {
		sys.p.brownian = forces[brownian].value;
		cout << indent+"Brownian force (in \"" << Dimensional::Unit::unit2suffix(forces[brownian].unit) << "\" units): " << sys.p.brownian << endl;
	}
	sys.p.kn = forces[kn].value;
	cout << indent+"Normal contact stiffness (in \"" << Dimensional::Unit::unit2suffix(forces[kn].unit) << "\" units): " << sys.p.kn << endl;
	sys.p.kt = forces[kt].value;
	cout << indent+"Sliding contact stiffness (in \"" << Dimensional::Unit::unit2suffix(forces[kt].unit) << "\" units): " << sys.p.kt << endl;
	sys.p.kr = forces[kr].value;
	cout << indent+"Rolling contact stiffness (in \"" << Dimensional::Unit::unit2suffix(forces[kr].unit) << "\" units): " << sys.p.kr << endl;

	if (forces.count(hydro) > 0) { // == if rate controlled
		sys.set_shear_rate(forces[hydro].value);
	}
}

void Simulation::setupNonDimensionalization(Dimensional::DimensionalValue<double> control_value){
	/**
	 \brief Non-dimensionalize the simulation.

		This function determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled), and converts all the input values in these units.
	 */
	if (control_var == ControlVariable::rate) {
		input_rate = control_value.value; // @@@ Renaming is required?
	}
	if (control_var == ControlVariable::rate || control_var == ControlVariable::viscnb) {
		units.add(Dimensional::Unit::hydro, control_value);
		internal_units = pickInternalUnitsRateControl();
	} else if (control_var == ControlVariable::stress) {
		units.add(Dimensional::Unit::stress, control_value);
		internal_units = control_value.unit;
	}
	output_units = control_value.unit;
	string indent = "  Simulation::\t";
	cout << indent << "internal units = " << Dimensional::Unit::unit2suffix(internal_units) << endl;
	cout << indent << "output units = " << Dimensional::Unit::unit2suffix(output_units) << endl;
	units.setInternalUnit(internal_units);
	exportForceAmplitudes();
	for (auto dimval: dimensional_input_params) {
		units.convertToInternalUnit(dimval.second);
	}
}

void Simulation::assertParameterCompatibility()
{
	// test for incompatibilities
	if (sys.brownian == true) {
		if (sys.pairwise_resistance && p.integration_method != 1) {
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (control_var == ControlVariable::stress) {
		if (p.integration_method != 0) {
			cerr << "Warning : use of the Predictor-Corrector method for the stress controlled simulation is experimental." << endl;
		}
		if (sys.brownian == true) {
			throw runtime_error("No Brownian stress-controlled simulations.");
		}
	}
	if (control_var == ControlVariable::viscnb) {
		if (sys.mobile_fixed) {
			throw runtime_error("Cannot run viscous number controlled simulations with fixed particles.");
		}
	}
	if (sys.critical_load) {
		p.friction_model = 2;
		cerr << "Warning : critical load simulation -> switched to friction_model=2" << endl;
	}
}

void Simulation::resolveTimeOrStrainParameters(const map<string, Dimensional::DimensionalValue<double>> &dim_params)
{
	/**
		\brief Interpret time units.

		We have to treat times as a special case, because we sometimes want to
		use times expressed in units of inverse shear rate (i.e. time=strain).
		Because the shear rate might change in the simulation, it is not possible
		to perform a simple change of unit at any time in the simulation.

		Ex 1: I am running a stress controlled simulation in repulsive time units.
		I set kn = 100r and time_end = 1kn;
		Then I know that I need to stop the simulation when sys.time()==1/100
		(i.e. when the time expressed in repulsive units reaches 1/100)

		Ex 2: I am running a stress controlled simulation in repulsive time units.
		I set kn = 100r but this time time_end = 1h;
		There is no way to know what 1h corresponds to in repulsive units,
		as the shear rate is evolving with time (stopping at sys.time()==???).
		The only way is to stop when sys.shear_strain()==1.

		Because these 2 cases lead to 2 different tests, we have to inform
		the System class which test to perform. If time_end > 0, System tests with
		System::time()==time_end. If time_end==-1, the test is done with
		System::shear_strain()==strain_end.

		We have to do this not only for time_end, but also for every time defined
		in the parameters.
	 */
	if (dim_params.at("time_end").dimension == Dimensional::Strain) {
		time_end = -1;
		strain_end = dim_params.at("time_end").value;
	} else {
		time_end = dim_params.at("time_end").value;
	}
	if (p.log_time_interval) {
		if (dim_params.at("time_end").dimension != dim_params.at("initial_log_time").dimension &&
			(dim_params.at("time_end").dimension == Dimensional::Strain || dim_params.at("initial_log_time").dimension == Dimensional::Strain)) {
			throw runtime_error(" If one of time_end or initial_log_time is a strain (\"h\" unit), than both must be.\n");
		}
	}
}

void Simulation::setConfiguration(bool binary_conf,
                                  string filename_import_positions)
{
	if (binary_conf) {
		auto format = getBinaryConfigurationFileFormat(filename_import_positions);
		ifstream file_import;
		file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
		if (!file_import) {
			ostringstream error_str;
			error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
			throw runtime_error(error_str.str());
		}

		switch(format) {
			case bin_format_base_new:
				{
					auto conf = readBinaryBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case bin_format_fixed_vel:
				{
					auto conf = readBinaryFixedVeloConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			default:
				{
					throw std::runtime_error("Unrecognized binary conf file format.");
				}
		}
	} else {
		auto format = getTxtConfigurationFileFormat(filename_import_positions);
		ifstream file_import;
		file_import.open(filename_import_positions.c_str());
		if (!file_import) {
			ostringstream error_str;
			error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
			throw runtime_error(error_str.str());
		}

		switch(format) {
			case txt_format_base_old:
				{
					auto conf = readTxtBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case txt_format_base_new:
				{
					auto conf = readTxtBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case txt_format_fixed_vel:
				{
					auto conf = readTxtFixedVeloConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case txt_format_circular_couette:
				{
					auto conf = readTxtCircularCouetteConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			default:
				{
					throw std::runtime_error("Unrecognized text conf file format.");
				}
		}
	}
}

void Simulation::setupSimulation(string in_args,
                                 vector<string>& input_files,
                                 bool binary_conf,
                                 Dimensional::DimensionalValue<double> control_value,
                                 string simu_identifier)
{
	/**
	 \brief Set up the simulation.

		This function is intended to be generically used to set up the simulation. It processes the input parameters, non-dimensionalizes the system and starts up a System class with the relevant parameters.
	 */
	string indent = "  Simulation::\t";
	cout << indent << "Simulation setup starting... " << endl;
	string filename_import_positions = input_files[0];
	string filename_parameters = input_files[1];

	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cout << "init_relax" << endl;
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	setDefaultParameters(control_value);
	readParameterFile(filename_parameters);
	if (p.impose_sigma_zz) {
		if (control_var == ControlVariable::stress) {
			throw runtime_error("controlling shear and normal stresses at once not implemented yet");
		} else {
			control_var = ControlVariable::viscnb;
		}
	}
	setupOptionalSimulation(indent);

	setupNonDimensionalization(control_value);
	assertParameterCompatibility();

	if (input_files[3] != "not_given") {
		throw runtime_error("pre-simulation data deprecated?");
	}
	resolveTimeOrStrainParameters(dimensional_input_params);
	setFromMap(p, dimensional_input_params);

	setConfiguration(binary_conf, filename_import_positions);

	p_initial = p;
	simu_name = prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
	                                  simu_identifier, control_value);
	openOutputFiles(simu_name);
	echoInputFiles(in_args, input_files);
	cout << indent << "Simulation setup [ok]" << endl;
}

void Simulation::autoSetParameters(const string &keyword, const string &value)
{
	/**
	 \brief Parse an input parameter
	 */
	string numeral, suffix;
	if (keyword == "simulation_mode") {
		p.simulation_mode = atoi(value.c_str());
	} else if (keyword == "lubrication_model") {
		p.lubrication_model = value;
	} else if (keyword == "friction_model") {
		p.friction_model = atoi(value.c_str());
	} else if (keyword == "repulsion") {
		units.add(Dimensional::Unit::repulsion, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "cohesion") {
		units.add(Dimensional::Unit::cohesion, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "brownian") {
		units.add(Dimensional::Unit::brownian, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "critical_load") {
		units.add(Dimensional::Unit::critical_load, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "repulsive_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "repulsive_max_length") {
		p.repulsive_max_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Time, value, keyword);
	} else if (keyword == "contact_relaxation_time_tan"){
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Time, value, keyword);
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "time_end") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::TimeOrStrain, value, keyword);
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max_gap") {
		p.lub_max_gap = atof(value.c_str());
	} else if (keyword == "interaction_range") {
		p.interaction_range = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		units.add(Dimensional::Unit::kn, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "kt") {
		units.add(Dimensional::Unit::kt, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "kr") {
		units.add(Dimensional::Unit::kr, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "dt") {
		p.dt = atof(value.c_str());
	} else if (keyword == "Pe_switch") {
		p.Pe_switch = atof(value.c_str());
	} else if (keyword == "mu_static") {
		p.mu_static = atof(value.c_str());
	} else if (keyword == "mu_dynamic") {
		p.mu_dynamic = atof(value.c_str());
	} else if (keyword == "mu_rolling") {
		p.mu_rolling = atof(value.c_str());
	} else if (keyword == "time_interval_output_config") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::TimeOrStrain, value, keyword);
	} else if (keyword == "time_interval_output_data") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::TimeOrStrain, value, keyword);
	} else if (keyword == "log_time_interval") {
		p.log_time_interval = str2bool(value);
	} else if (keyword == "initial_log_time") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::TimeOrStrain, value, keyword);
	} else if (keyword == "nb_output_data_log_time") {
		p.nb_output_data_log_time = atoi(value.c_str());
	} else if (keyword == "nb_output_config_log_time") {
		p.nb_output_config_log_time = atoi(value.c_str());
	} else if (keyword == "out_data_particle") {
		p.out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		p.out_data_interaction = str2bool(value);
	} else if (keyword == "out_data_vel_components") {
		p.out_data_vel_components = str2bool(value);
	} else if (keyword == "origin_zero_flow") {
		p.origin_zero_flow = str2bool(value);
	} else if (keyword == "auto_determine_knkt") {
		p.auto_determine_knkt = str2bool(value.c_str());
	} else if (keyword == "overlap_target") {
		p.overlap_target = atof(value.c_str());
	} else if (keyword == "disp_tan_target") {
		p.disp_tan_target = atof(value.c_str());
	} else if (keyword == "memory_strain_avg") {
		p.memory_strain_avg = atof(value.c_str());
	} else if (keyword == "memory_strain_k") {
		p.memory_strain_k = atof(value.c_str());
	} else if (keyword == "start_adjust") {
		p.start_adjust = atof(value.c_str());
	} else if (keyword == "min_kn") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword);
	} else if (keyword == "max_kn") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword);
	} else if (keyword == "min_kt") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword);
	} else if (keyword == "max_kt") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword);
	} else if (keyword == "min_dt") {
		p.min_dt = atof(value.c_str());
	} else if (keyword == "max_dt") {
		p.max_dt = atof(value.c_str());
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		units.add(Dimensional::Unit::ft_max, Dimensional::str2DimensionalValue(Dimensional::Force, value, keyword));
	} else if (keyword == "fixed_dt") {
		p.fixed_dt = str2bool(value);
	} else if (keyword == "theta_shear") {
		p.theta_shear = atof(value.c_str());
		p.theta_shear *= M_PI/180;  // convert in radians
	} else if (keyword == "strain_reversal") {
		p.strain_reversal = atof(value.c_str());
	} else if (keyword == "event_handler") {
		p.event_handler = value;
		p.event_handler.erase(remove(p.event_handler.begin(), p.event_handler.end(), '\"' ), p.event_handler.end());
	} else if (keyword == "out_particle_stress") {
		p.out_particle_stress = value;
		p.out_particle_stress.erase(remove(p.out_particle_stress.begin(), p.out_particle_stress.end(), '\"' ), p.out_particle_stress.end());
	} else if (keyword == "out_binary_conf") {
		p.out_binary_conf = str2bool(value);
	} else if (keyword == "np_fixed") {
		p.np_fixed = atoi(value.c_str());
	} else if (keyword == "keep_input_strain") {
		p.keep_input_strain = str2bool(value);
	} else if (keyword == "sigma_zz") {
		units.add(Dimensional::Unit::sigma_zz, Dimensional::str2DimensionalValue(Dimensional::Stress, value, keyword));
	} else if (keyword == "impose_sigma_zz") {
		p.impose_sigma_zz = str2bool(value);
	} else {
		ostringstream error_str;
		error_str  << "keyword " << keyword << " is not associated with an parameter" << endl;
		throw runtime_error(error_str.str());
	}
}

void Simulation::readParameterFile(const string& filename_parameters)
{
	/**
	 \brief Read and parse the parameter file
	 */
	ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		ostringstream error_str;
		error_str  << " Parameter file '" << filename_parameters << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	string keyword, value;
	while (!fin.eof()) {
		string line;
		if (!getline(fin, line, ';')) {
			break;
		}
		if (fin.eof()) {
			break;
		}
		string str_parameter;
		removeBlank(line);
		str_parameter = line;
		string::size_type begin_comment;
		string::size_type end_comment;
		do {
			begin_comment = str_parameter.find("/*");
			end_comment = str_parameter.find("*/");
			if (begin_comment > 10000) {
				break;
			}
			str_parameter = str_parameter.substr(end_comment+2);
		} while (true);
		if (begin_comment > end_comment) {
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			throw runtime_error("syntax error in the parameter file.");
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if (pos_slashslash != string::npos) {
			throw runtime_error(" // is not the syntax to comment out. Use /* comment */");
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void Simulation::setDefaultParameters(Dimensional::DimensionalValue<double> control_value)
{

	/**
	 \brief Set default values for ParameterSet parameters.
	 */
	auto input_scale = Dimensional::Unit::unit2suffix(control_value.unit);
	autoSetParameters("Pe_switch", "5");
	autoSetParameters("dt", "1e-4");
	autoSetParameters("disp_max", "1e-3");
	autoSetParameters("monolayer", "false");
	autoSetParameters("rest_threshold", "1e-4");
	autoSetParameters("integration_method", "1");
	autoSetParameters("np_fixed", "0");
	autoSetParameters("sd_coeff", "1");
	autoSetParameters("lubrication_model", "tangential");
	autoSetParameters("friction_model", "1");
	autoSetParameters("time_end", "10h");
	autoSetParameters("lub_max_gap", "0.5");
	autoSetParameters("interaction_range", "-1");
	autoSetParameters("lub_reduce_parameter", "1e-3");
	autoSetParameters("contact_relaxation_time", "1e-3"+input_scale);
	autoSetParameters("contact_relaxation_time_tan", "-1"+input_scale);
	if (input_scale != "kn") {
		autoSetParameters("kn", "2000"+input_scale);
		autoSetParameters("min_kn", "1000"+input_scale);
		autoSetParameters("max_kn", "1000000"+input_scale);
	}
	if (input_scale != "kt") {
		autoSetParameters("kt", "0.5kn");
		autoSetParameters("min_kt", "1000"+input_scale);
		autoSetParameters("max_kt", "1000000"+input_scale);
	}
	if (input_scale != "kr") {
		autoSetParameters("kr", "0kn");
	}
	autoSetParameters("auto_determine_knkt", "false");
	autoSetParameters("overlap_target", "0.05");
	autoSetParameters("disp_tan_target", "0.05");
	autoSetParameters("memory_strain_avg", "0.01");
	autoSetParameters("memory_strain_k", "0.02");
	autoSetParameters("start_adjust", "0.2");
	autoSetParameters("repulsive_length", "0.05");
	autoSetParameters("repulsive_max_length", "-1");
	autoSetParameters("mu_static", "1");
	autoSetParameters("mu_dynamic", "-1");
	autoSetParameters("mu_rolling", "0");
	autoSetParameters("time_interval_output_data", "1e-2h");
	autoSetParameters("time_interval_output_config", "1e-1h");
	autoSetParameters("log_time_interval", "false");
	autoSetParameters("initial_log_time", "1e-4h");
	autoSetParameters("nb_output_data_log_time", "100");
	autoSetParameters("nb_output_config_log_time", "100");
	autoSetParameters("origin_zero_flow", "true");
	autoSetParameters("out_data_particle", "true");
	autoSetParameters("out_data_interaction", "true");
	autoSetParameters("out_particle_stress", "");
	autoSetParameters("out_binary_conf", "false");
	autoSetParameters("out_data_vel_components", "false");
	autoSetParameters("fixed_dt", "false");
	autoSetParameters("theta_shear", "0");
	autoSetParameters("event_handler", "");
	autoSetParameters("simulation_mode", "0");
	autoSetParameters("keep_input_strain", "false");
	autoSetParameters("impose_sigma_zz", "false");
}


inline string columnDefinition(int &cnb, const string &type, const string &name)
{
	stringstream defs;
	if (type == "vec3d") {
		array<string, 3> xyz = {"x", "y", "z"};
		for (auto &u : xyz) {
			stringstream col_def_complement;
			defs << "#" << cnb << ": "<< name << " " << u << "\n";
			cnb ++;
		}
	} else if (type=="scalar") {
		defs << "#" << cnb << ": "<< name << "\n";
	} else {
		throw runtime_error(" unknown type for column def\n");
	}
	return defs.str();
}

void Simulation::openOutputFiles(string simu_name)
{
	/**
	 \brief Set up the output files

		This function determines a simulation name from the parameters, opens the output files with the corresponding name and prints their header.
	 */

	stringstream data_header;
	createDataHeader(data_header);
	outdata.setFile("data_"+simu_name+".dat", data_header.str(), force_to_run);
	outdata_st.setFile("st_"+simu_name+".dat", data_header.str(), force_to_run);
	if (!p.out_particle_stress.empty()) {
		outdata_pst.setFile("pst_"+simu_name+".dat", data_header.str(), force_to_run);
	}
	string time_filename = "t_"+simu_name+".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_"+simu_name+".dat";
	fout_input.open(input_filename.c_str());
	if (p.out_data_particle) {
		outdata_par.setFile("par_"+simu_name+".dat", data_header.str(), force_to_run);
	}
	if (p.out_data_interaction) {
		outdata_int.setFile("int_"+simu_name+".dat", data_header.str(), force_to_run);
	}
}

string Simulation::prepareSimulationName(bool binary_conf,
                                         const string& filename_import_positions,
                                         const string& filename_parameters,
                                         const string& simu_identifier,
                                         Dimensional::DimensionalValue<double> control_value)
{
	/**
	 \brief Determine simulation name.
	 */
	ostringstream ss_simu_name;
	string::size_type pos_name_end = filename_import_positions.find_last_of(".");
	string::size_type param_name_end = filename_parameters.find_last_of(".");
	string::size_type pos_name_start;
	if (binary_conf) { // TO DO: improve name generation for binary input
		pos_name_start = filename_import_positions.find_last_of("/");
	} else {
		pos_name_start = filename_import_positions.find_last_of("/");
	}
	string::size_type param_name_start = filename_parameters.find_last_of("/");
	if (pos_name_start == std::string::npos) {
		pos_name_start = -1;
	}
	if (param_name_start == std::string::npos) {
		param_name_start = -1;
	}
	pos_name_start += 1;
	param_name_start += 1;
	ss_simu_name << filename_import_positions.substr(pos_name_start, pos_name_end-pos_name_start);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(param_name_start, param_name_end-param_name_start);
	ostringstream string_control_parameters;
	if (long_file_name) {
		for (const auto& f: dimensional_input_params) {
			if (f.first.find("stiffness") == std::string::npos) {
				string_control_parameters << "_" << f.first << f.second.value << Dimensional::Unit::unit2suffix(f.second.unit);
			}
		}
	}
	if (control_var==ControlVariable::rate || control_var==ControlVariable::viscnb) {
		string_control_parameters << "_" << "rate";
	}
	if (control_var==ControlVariable::stress) {
		string_control_parameters << "_" << "stress";
	}
	string_control_parameters << control_value.value << Dimensional::Unit::unit2suffix(control_value.unit);
	ss_simu_name << string_control_parameters.str();
	if (simu_identifier != "") {
		ss_simu_name << "_";
		ss_simu_name << simu_identifier;
	}
	string indent = "  Simulation::\t";
	cout << indent << "filename: " << ss_simu_name.str() << endl;

	return ss_simu_name.str();
}

TimeKeeper Simulation::initTimeKeeper() {
	TimeKeeper tk;
	if (p.log_time_interval) {
		tk.addClock("data", LogClock(p.initial_log_time,
		            p.time_end,
		            p.nb_output_data_log_time,
		            dimensional_input_params["time_end"].dimension == Dimensional::Strain));
	} else {
		tk.addClock("data", LinearClock(p.time_interval_output_data,
		            dimensional_input_params["time_interval_output_data"].dimension == Dimensional::Strain));
	}
	if (p.log_time_interval) {
		tk.addClock("config", LogClock(p.initial_log_time,
		            p.time_end,
		            p.nb_output_config_log_time,
		            dimensional_input_params["time_end"].dimension == Dimensional::Strain));
	} else {
		tk.addClock("config", LinearClock(p.time_interval_output_config,
		            dimensional_input_params["time_interval_output_config"].dimension == Dimensional::Strain));
	}
	return tk;
}
