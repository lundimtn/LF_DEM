//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include "Simulation.h"
#ifndef GIT_VERSION
#include "VersionInfo.h"
#endif
/*
 * VersionInfo.h is automatically generated
 * before compiling source codes.
 * In Xcode, the following script is run in the Pre-Action of Build.
 * -----------------------------------
 * git=/usr/bin/git
 * cd ${PROJECT_DIR}/LF_DEM
 * version=`$git describe --dirty`
 * echo "#define GIT_VERSION \"$version\"" > VersionInfo.h
 * -----------------------------------
 */
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <cctype>

Simulation::Simulation():
shear_rate_expectation(-1),
unit_scales("hydro")
{};

Simulation::~Simulation()
{
	if (fout_rheo.is_open()) {
		fout_rheo.close();
	}
	if (fout_particle.is_open()) {
		fout_particle.close();
	}
	if (fout_interaction.is_open()) {
		fout_interaction.close();
	}
};

void Simulation::contactForceParameter(string filename)
{
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." << endl;
		exit(1);
	}
	
	// temporal variables to keep imported values.
	double phi_, kn_, kt_, dt_;
	// To find parameters for considered volume fraction phi.
	bool found = false;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	
	if (found) {
		// Set the parameter object
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		cerr << " Input for kn, kt, dt = " << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_ << endl;
	} else {
		cerr << " Error: file " << filename.c_str() << " contains no data for vf = " << phi_ << endl;
		exit(1);
	}
}

void Simulation::contactForceParameterBrownian(string filename)
{
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." <<endl;
		exit(1);
	}
	
	// temporal variables to keep imported values.
	double phi_, peclet_, kn_, kt_, dt_;
	bool found = false;
	while (fin_knktdt >> phi_ >> peclet_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction && peclet_ == dimensionless_numbers["b"]) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	
	if (found) {
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		cout << "Input for vf = " << phi_ << " and Pe = " << peclet_ << " : kn = " << kn_ << ", kt = " << kt_ << " and dt = " << dt_ << endl;
	} else {
		cerr << " Error: file " << filename.c_str() << " contains no data for vf = " << volume_or_area_fraction << " and Pe = " << dimensionless_numbers["b"] << endl;
		exit(1);
	}
}

void Simulation::importPreSimulationData(string filename)
{
	ifstream fin_PreSimulationData;
	fin_PreSimulationData.open(filename.c_str());
	if (!fin_PreSimulationData) {
		cerr << " Pre-simulation data file '" << filename << "' not found." << endl;
		exit(1);
	}
	
	double stress_, shear_rate_;
	while (fin_PreSimulationData >> stress_ >> shear_rate_) {
		if (stress_ == sys.target_stress_input) {
			break;
		}
	}
	shear_rate_expectation = shear_rate_;
}

void Simulation::echoInputFiles(string in_args, vector<string> &input_files)
{
	fout_input << "# LF_DEM version " << GIT_VERSION << ", called with:" << endl;
	fout_input << in_args << endl << endl;
	for (const string &in_file : input_files) {
		ifstream in_f;
		string line;
		in_f.open(in_file.c_str());
		if(in_f.is_open()){
			fout_input << "********** File " << in_file << " ************" << endl << endl;
			while(in_f.good()){
				getline(in_f, line);
				fout_input << line << endl;
			}
			fout_input << endl << endl;
		}
		in_f.close();
	}
	fout_input.close();
}

// void Simulation::setUnitScalesBrownian(double dimensionless_number)
// {
// 	sys.dimensionless_number = dimensionless_number; // Peclet number
// 	if (sys.dimensionless_number > p.Pe_switch && !sys.zero_shear) {
// 		unit_scales = "hydro";
// 		sys.amplitudes.sqrt_temperature = 1/sqrt(sys.dimensionless_number);
// 		sys.set_shear_rate(1);
// 		if (p.repulsiveforce == false
// 			&& p.cohesion == false
// 			&& p.critical_load == false
// 			&& p.magnetic == false) {
// 			cerr << "Only Brownian" << endl;
// 			string_control_parameters << "_p" << sys.dimensionless_number;
// 		} else if (p.repulsiveforce == true
// 				   && p.critical_load == false
// 				   && p.cohesion == false
// 				   && p.magnetic == false) {
// 			cerr << "Repulsive force, ratio to kT/a^3 : " << p.ratio_repulsion << endl;
// 			/* When both Brownian and repulsive forces exist
			 
// 			 * `ratio_repulsion' = F_rep(0)/(kT/a)
// 			 * Filename includes "rXXXX_pXXXX".
// 			 */
// 			sys.amplitudes.repulsion = p.ratio_repulsion/sys.dimensionless_number;
// 			string_control_parameters << "_r" << p.ratio_repulsion  << "_p" << sys.dimensionless_number;
// 		} else if (p.critical_load == true
// 				   && p.repulsiveforce == false
// 				   && p.cohesion == false
// 				   && p.magnetic == false) {
// 			cerr << "Critical load, ratio to kT/a^3 : " << p.ratio_critical_load << endl;
// 			sys.critical_normal_force = p.ratio_critical_load/sys.dimensionless_number;
// 			p.friction_model = 2;
// 			string_control_parameters << "_c" << p.ratio_critical_load << "_p" << sys.dimensionless_number;
// 		} else if (p.cohesion == true
// 				   && p.repulsiveforce == false
// 				   && p.critical_load == false
// 				   && p.magnetic == false) {
// 			cerr << "Cohesion, ratio to kT/a^3 : " << p.ratio_cohesion << endl;
// 			sys.dimensionless_cohesive_force = p.ratio_cohesion/sys.dimensionless_number;
// 			string_control_parameters << "_a" << p.ratio_cohesion << "_p" << sys.dimensionless_number;
// 		} else if (p.magnetic == true
// 				   && p.repulsiveforce == false
// 				   && p.cohesion == false
// 				   && p.critical_load == false) {
// 			cerr << "Brownian Magnetic" << endl;
// 			string_control_parameters << "_mag_p" << sys.dimensionless_number;
// 		} else {
// 			cerr << "repulsiveforce: " << p.repulsiveforce << endl;
// 			cerr << "critical_load: " << p.critical_load << endl;
// 			cerr << "cohesion: " << p.cohesion << endl;
// 			cerr << "not yet implemented" << endl;
// 			exit(1);
// 		}
// 	} else {
// 		unit_scales = "thermal";
// 		sys.amplitudes.sqrt_temperature = 1;
// 		sys.set_shear_rate(sys.dimensionless_number);
// 		if (p.repulsiveforce == false
// 			&& p.cohesion == false
// 			&& p.critical_load == false
// 			&& p.magnetic == false) {
// 			cerr << "Only Brownian" << endl;
// 			string_control_parameters << "_p" << sys.dimensionless_number;
// 		} else if (p.repulsiveforce == true
// 				   && p.critical_load == false
// 				   && p.cohesion == false
// 				   && p.magnetic == false) {
// 			cerr << "Repulsive force, ratio to kT/a^3 : " << p.ratio_repulsion << endl;
// 			/* When both Brownian and repulsive forces exist
// 			 *
// 			 * `ratio_repulsion' = F_rep(0)/(kT/a)
// 			 * Filename includes "rXXXX_pXXXX".
// 			 */
// 			sys.amplitudes.repulsion = p.ratio_repulsion;
// 			string_control_parameters << "_r" << p.ratio_repulsion << "_p" << sys.dimensionless_number;
// 		} else if (p.critical_load == true
// 				   && p.repulsiveforce == false
// 				   && p.cohesion == false
// 				   && p.magnetic == false) {
// 			cerr << "Critical load, ratio to kT/a^3 : " << p.ratio_critical_load << endl;
// 			sys.critical_normal_force = p.ratio_critical_load;
// 			p.friction_model = 2;
// 			string_control_parameters << "_c" << p.ratio_critical_load << "_p" << sys.dimensionless_number;
// 		} else if (p.cohesion == true
// 				   && p.repulsiveforce == false
// 				   && p.critical_load == false
// 				   && p.magnetic == false) {
// 			cerr << "Cohesion, ratio to kT/a^3 : " << p.ratio_cohesion << endl;
// 			sys.cohesive_force = p.ratio_cohesion;
// 			sys.dimensionless_cohesive_force = sys.cohesive_force;
// 			string_control_parameters << "_a" << p.ratio_cohesion << "_p" << sys.dimensionless_number;
// 			cerr << "dimensionless_cohesive_force : " << sys.dimensionless_cohesive_force << endl;
// 		} else if (p.magnetic == true
// 				   && p.repulsiveforce == false
// 				   && p.cohesion == false
// 				   && p.critical_load == false) {
// 			cerr << "Brownian Magnetic" << endl;
// 			string_control_parameters << "_mag_p" << sys.dimensionless_number;
// 		} else {
// 			cerr << "repulsiveforce: " << p.repulsiveforce << endl;
// 			cerr << "critical_load: " << p.critical_load << endl;
// 			cerr << "cohesion: " << p.cohesion << endl;
// 			cerr << "not yet implemented" << endl;
// 			exit(1);
// 		}

// 	}

// }


// void Simulation::setUnitScalesNonBrownianRate(double dimensionlessnumber)
// {
// 	if (p.repulsiveforce == true
// 		&& p.critical_load == false
// 		&& p.cohesion == false
// 		&& p.magnetic == false) {
// 		cerr << "Repulsive force, shear rate (in units of F_R(0)/(6 pi eta_0 a^2)): " << dimensionlessnumber << endl; //@???
// 		sys.dimensionless_number = dimensionlessnumber;
// 		sys.amplitudes.repulsion = 1/sys.dimensionless_number;
// 		string_control_parameters << "_r" << sys.dimensionless_number;
// 	} else if (p.repulsiveforce == false
// 			   && p.critical_load == false
// 			   && p.cohesion == false
// 			   && p.magnetic == false) {
// 		cerr << "Infinite shear rate (quasi-Newtonian)" << endl;
// 		sys.dimensionless_number = -1;
// 		sys.amplitudes.repulsion = 0;
// 		string_control_parameters << "_quasi_Newtonian";
// 	} else if (p.critical_load == true
// 			   && p.repulsiveforce == false
// 			   && p.cohesion == false
// 			   && p.magnetic == false) {
// 		cerr << "Critical load, shear rate (in units of F_R(0)/(6 pi eta_0 a^2)): " << dimensionlessnumber << endl;
// 		sys.dimensionless_number = dimensionlessnumber;
// 		p.friction_model = 2;
// 		sys.critical_normal_force = 1/sys.dimensionless_number;
// 		string_control_parameters << "_c" << sys.dimensionless_number;
// 	} else if (p.cohesion == true
// 			   && p.repulsiveforce == false
// 			   && p.critical_load == false
// 			   && p.magnetic == false) {
// 		cerr << "Cohesive force, sheafr rate (in units of F_R(0)/6 pi eta_0 a^2)): " << dimensionlessnumber << endl;
// 		sys.dimensionless_number = dimensionlessnumber;
// 		/* In the rate control simulation,
// 		 * dimensionless_cohesive_force can be given.
// 		 */
// 		string_control_parameters << "_a" << dimensionlessnumber;
// 	} else if (p.repulsiveforce == true
// 			   && p.cohesion == true
// 			   && p.critical_load == false
// 			   && p.magnetic == false) {
// 		cerr << "Repulsive force + Cohesive force" << endl;
// 		sys.dimensionless_number = dimensionlessnumber;
// 		sys.amplitudes.repulsion = 1/sys.dimensionless_number;
// 		sys.cohesive_force = p.ratio_cohesion;
// 		string_control_parameters << "_a" <<  p.ratio_cohesion << "_r" << sys.dimensionless_number;
// 	} else if (p.magnetic == true
// 			   && p.repulsiveforce == false
// 			   && p.cohesion == false
// 			   && p.critical_load == false) {
// 		cerr << "Magnetic interaction" << endl;
// 		sys.dimensionless_number = dimensionlessnumber;
// 		sys.amplitudes.magnetic = 1/sys.dimensionless_number;
// 		string_control_parameters << "_mag" << sys.dimensionless_number << "_r" << sys.dimensionless_number;

// 	} else {
// 		cerr << "strain -> non-Brownian -> ???" << endl;
// 		exit(1);
// 	}
// }

// void Simulation::setUnitScalesStressControlled(double dimensionlessnumber){
//   if (p.repulsiveforce == false && p.cohesion == false) {
// 	cerr << " Stress controlled simulations need a repulsive force ! " << endl;
// 	cerr << " ===> This is not correct. We can make stress controlled simulation without any additional force." << endl;
// 	exit(1);
//   } else {
// 	if (p.repulsiveforce == true
// 		&& p.critical_load == false
// 		&& p.cohesion == false) {
// 	  cerr << "Repulsive force" << endl;
//  			sys.amplitudes.repulsion = 1;
//  			sys.target_stress_input = dimensionlessnumber;
//  			sys.target_stress = sys.target_stress_input/6/M_PI;
//  			string_control_parameters << "_s" << sys.target_stress_input;
// 	} else if (p.cohesion == true
// 				   && p.repulsiveforce == false
// 				   && p.critical_load == false) {
// 			cerr << "Cohesive force" << endl;
// 			p.unscaled_contactmodel = false;
// 			sys.cohesive_force = 1;
// 			sys.target_stress_input = dimensionlessnumber;
// 			sys.target_stress = sys.target_stress_input/6/M_PI;
// 			/* Initial relaxation for stress control simulation.
// 			 * (To avoid breaking bonds due to startup flows.)
// 			 */
// 			sys.init_strain_shear_rate_limit = -9999;
// 			sys.init_shear_rate_limit = 9999;
// 			string_control_parameters << "_b" << sys.target_stress_input;
// 		} else if (p.cohesion == true
// 				   && p.repulsiveforce == true
// 				   && p.critical_load == false) {
// 			string_control_parameters << "_b" <<  p.ratio_cohesion << "_r" << sys.dimensionless_number;
// 			cerr << "not yet implemented" << endl;
// 			exit(1);
// 		} else {
// 			cerr << "stress -> non-Brownian -> ???" << endl;
// 			exit(1);
// 		}
// 		sys.dimensionless_number = 1; // needed for 1st time step
// 		/* The target stress (``ratio_repulsion'') is given trough the command argument
// 		 * with an unit stres: eta_0*gammmadot_0.
// 		 * However, in the code, sys.target_stress is computed as an unit F_rep/a^2.
// 		 */
// 	}
// }


void Simulation::determineDimensionlessNumbers(double dimensionlessnumber, string rate_unit)
{
	// determine the dimensionless numbers
	
	string force_type = rate_unit; // our force defining the shear rate
	dimensionless_numbers[force_type] = dimensionlessnumber;
	if(values[force_type]>0){
		cerr << "Error: redefinition of the rate (given both in the command line and in the parameter file with \"" << force_type << "\" force)" << endl; exit(1);
	}
	// switch this force in hydro units
	values[force_type] = 1/dimensionless_numbers[force_type];
	suffixes[force_type] = "h";

	
	// now resolve the other force units
	set <string> resolved_units;
	resolved_units.clear();

	// already done the one which gives the shear rate
	resolved_units.insert(force_type);

	// some of them are already in hydro
	for(auto&& f: suffixes){
	  force_type = f.first;
	  string suffix = f.second;
	  if ( suffix == "h" ){
		dimensionless_numbers[force_type] = 1./values[force_type];
		resolved_units.insert(force_type);
	  }
	}

	// now the "non-trivial" ones, we solve iteratively
	unsigned int resolved = resolved_units.size();
	unsigned int previous_resolved;
	do{
	  previous_resolved = resolved;
	  for(auto&& f: suffixes){
		force_type = f.first;
		string suffix = f.second;
		if (resolved_units.find(suffix) != resolved_units.end()){  // then we know how to convert to hydro
		  values[force_type] /= dimensionless_numbers[suffix];
		  suffixes[force_type] = "h";
		  dimensionless_numbers[force_type] = 1./values[force_type];
		  resolved_units.insert(force_type);
		}
	  }
	  resolved = resolved_units.size();
	}while(previous_resolved < resolved);

	// check we found everyone
	if( resolved < suffixes.size() ){
	  for(auto&& f: suffixes){
		force_type = f.first;
		string suffix = f.second;
		if (resolved_units.find(suffix) == resolved_units.end()){
		  cerr << "Error: force type \"" << force_type << "\" has an unknown scale \"" << suffix << "\"" << endl;
		}
	  }
	  exit(1);
	}
	  
}

void Simulation::setLowPeclet()
{
  sys.lowPeclet = true;
  double scale_factor_SmallPe = p.Pe_switch/dimensionless_numbers["b"];
  p.memory_strain_k /= scale_factor_SmallPe;
  p.memory_strain_avg /= scale_factor_SmallPe;
  p.start_adjust /= scale_factor_SmallPe;
  p.dt *= p.Pe_switch; // to make things continuous at Pe_switch
}


void Simulation::convertForceValues()
{
	double converter = 1;
	if (unit_scales == "thermal") {
		converter = dimensionless_numbers["b"];
	}
	if (unit_scales == "repulsive") {
		converter = dimensionless_numbers["r"];
	}
	for(auto&& f: suffixes){
		string force_type = f.first;
		string suffix = f.second;
		values[force_type] *= converter;
	}	
}

void Simulation::setUnitScale()
{
  bool is_brownian = dimensionless_numbers.find("b") != dimensionless_numbers.end();
  if(is_brownian){
	  //	  sys.dimensionless_number = dimensionless_numbers["b"]; // Peclet number
	  if (dimensionless_numbers["b"] > p.Pe_switch && !sys.zero_shear) {
		unit_scales = "hydro";
		sys.amplitudes.sqrt_temperature = 1/sqrt(dimensionless_numbers["b"]);
		sys.set_shear_rate(1);
	  }
	  else{ // low Peclet mode
 		unit_scales = "thermal";
 		sys.amplitudes.sqrt_temperature = 1;
 		sys.set_shear_rate(dimensionless_numbers["b"]);
		setLowPeclet();
	  }
  }
  else{
	unit_scales = "hydro";
	sys.set_shear_rate(1);
  }
  
  // convert from hydro scale to chosen scale
  convertForceValues();

  if(is_brownian){
	sys.brownian = true;
	p.brownian_amplitude = values["b"];
	cerr << "Brownian, Peclet number " << dimensionless_numbers["b"] << endl;
  }
  else{
	cerr << "non-Brownian" << endl;
  }
  bool is_repulsive = dimensionless_numbers.find("r") != dimensionless_numbers.end();
  if(is_repulsive){
	sys.repulsiveforce = true;
	sys.amplitudes.repulsion = values["r"];
  }
  bool is_critical_load = dimensionless_numbers.find("cl") != dimensionless_numbers.end();
  if(is_critical_load){
	sys.critical_load = true;
	sys.amplitudes.critical_normal_force = values["cl"];
  }
  bool is_cohesive = dimensionless_numbers.find("c") != dimensionless_numbers.end();
  if(is_cohesive){
	sys.cohesion = true;
	sys.amplitudes.cohesion = values["c"];
  }
  bool is_magnetic = dimensionless_numbers.find("m") != dimensionless_numbers.end();
  if(is_magnetic){
	sys.magnetic = true;
	p.magnetic_amplitude = values["m"];
  }
  bool is_ft_max = dimensionless_numbers.find("ft") != dimensionless_numbers.end();
  if(is_ft_max){
	sys.amplitudes.ft_max = values["ft"];
  }


}


void Simulation::setupSimulationSteadyShear(string in_args,
											vector<string> &input_files,
											bool binary_conf,
											double dimensionlessnumber,
											string input_scale,
											string control_variable)
{
	control_var = control_variable;
	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];
	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cerr << "init_relax" << endl;
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	
	setDefaultParameters();
	readParameterFile();

	if (control_var == "rate") {
	  determineDimensionlessNumbers(dimensionlessnumber, input_scale);	
	  setUnitScale();
	}
	else if (control_var == "stress") {
		p.unscaled_contactmodel = true;
		unit_scales = "repulsion";
		sys.amplitudes.repulsion = 1;
		sys.set_shear_rate(1);
		bool is_critical_load = values.find("cl") != values.end();
		if (is_critical_load) {
			cerr << " Stress controlled simulations for CLM not implemented ! " << endl;
			exit(1);
		}
		cerr << " stress controlled temporarily disabled " << endl; exit(1);
		//		setUnitScalesNonBrownianStress(dimensionlessnumber);
	}

	// test for incompatibilities
	if (sys.brownian == true) {
		if (p.integration_method != 1) {
			cerr << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			cerr << "Modify the parameter file: " << filename_parameters << endl;
			exit(1);
		}
	}
	if (control_var == "stress") {
		if (p.integration_method != 0) {
			cerr << "Must be Euler method for stress controlled simulation" << endl;
		}
		p.integration_method = 0;
	}
	if(sys.critical_load){
		p.friction_model = 2;
	}
	
	if (binary_conf) {
		importConfigurationBinary();
	} else {
		importInitialPositionFile();
	}
	if (initial_lees_edwards_disp > 0) {
		sys.shear_disp = initial_lees_edwards_disp;
	} else {
		sys.shear_disp = 0;
	}
	
	if (input_files[2] != "not_given") {
		if (sys.brownian && !p.auto_determine_knkt) {
			contactForceParameterBrownian(input_files[2]);
		} else {
			contactForceParameter(input_files[2]);
		}
	}
	
	if (input_files[3] != "not_given") {
		importPreSimulationData(input_files[3]);
		time_interval_output_data = p.time_interval_output_data/shear_rate_expectation;
		time_interval_output_config = p.time_interval_output_config/shear_rate_expectation;
	} else {
		time_interval_output_data = p.time_interval_output_data;
		time_interval_output_config = p.time_interval_output_config;
	}
	
	cerr << "  time_interval_output_data = " << time_interval_output_data << endl;
	cerr << "  time_interval_output_config = " << time_interval_output_config << endl;
	
	//exportParameterSet();
	sys.importParameterSet(p);

	if (sys.brownian) {
		sys.setupBrownian();
	}
	sys.setupSystem(control_var);

	openOutputFiles(binary_conf);
	echoInputFiles(in_args, input_files);
}

/*
 * Main simulation
 */
void Simulation::simulationSteadyShear(string in_args,
									   vector<string> &input_files,
									   bool binary_conf,
									   double dimensionless_number,
									   string input_scale,
									   string control_variable)
{
	user_sequence = false;
	control_var = control_variable;
	setupSimulationSteadyShear(in_args, input_files, binary_conf, dimensionless_number, input_scale, control_var);
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	//	double strain_output_data = 0;
	double strain_output_config = 0;
	double time_output_data = 0;
	double time_output_config = 0;
	if (sys.cohesion) {
		sys.new_contact_gap = 0.02;
	} else {
		sys.new_contact_gap = 0;
	}
	int jammed = 0;
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;

	/******************** OUTPUT INITIAL DATA ********************/
	evaluateData();
	outputRheologyData();
	outputStressTensorData();
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	while (sys.get_time() < p.time_end-1e-8) {
		time_output_data = cnt_simu_loop*time_interval_output_data;
		time_output_config = cnt_config_out*time_interval_output_config;
		sys.timeEvolution(time_output_data);
		cnt_simu_loop ++;

		/******************** OUTPUT DATA ********************/
		evaluateData();
		outputRheologyData();
		outputStressTensorData();
		outputConfigurationBinary();
		if (time_interval_output_data == -1) {
			if (sys.get_shear_strain() >= strain_output_config-1e-8) {
				outputConfigurationData();
				cnt_config_out ++;
			}
		} else {
			if (sys.get_time() >= time_output_config-1e-8) {
				outputConfigurationData();
				cnt_config_out ++;
			}
		}
		/*****************************************************/
		
		cerr << "time: " << sys.get_time() << " / " << p.time_end << endl;
		if (!sys.zero_shear
			&& abs(sys.get_shear_rate()) < p.rest_threshold){
			cerr << "shear jamming " << jammed << endl;
			jammed ++;
			if (jammed > 10) {
				cerr << "shear jamming";
				break;
			}
		} else {
			jammed = 0;
		}
		sys.new_contact_gap = 0;
		if (time_strain_1 == 0 && sys.get_shear_strain() > 1) {
			now = time(NULL);
			time_strain_1 = now;
			timestep_1 = sys.get_total_num_timesteps();
		}
	}
	now = time(NULL);
	time_strain_end = now;
	timestep_end = sys.get_total_num_timesteps();
	outputComputationTime();
	if (filename_parameters.find("init_relax", 0)) {
		/* To prepare relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to export the position data.
		 */
		outputFinalConfiguration();
	}
}

void Simulation::outputComputationTime()
{
	int time_from_1 = time_strain_end-time_strain_1;
	int time_from_0 = time_strain_end-time_strain_0;
	int timestep_from_1 = timestep_end-timestep_1;
	fout_time << "# np time_from_0 time_from_1 timestep_end timestep_from_1" << endl;
	fout_time << sys.get_np() << ' ';
	fout_time << time_from_0 << ' ';
	fout_time << time_from_1 << ' ';
	fout_time << timestep_end << ' ';
	fout_time << timestep_from_1 << endl;
}

/*
 * Main simulation
 */
// void Simulation::simulationUserDefinedSequence(string seq_type,
// 											   string in_args,
// 											   vector<string> &input_files,
// 											   bool binary_conf,
// 											   string control_variable)
// {
// 	user_sequence = true;
// 	control_var = control_variable;
// 	filename_import_positions = input_files[0];
// 	filename_parameters = input_files[1];
// 	filename_sequence = input_files[4];
// 	string::size_type pos_ext_sequence = filename_sequence.find(".dat");
// 	sys.brownian = false;
// 	sys.target_stress_input = 0;
// 	sys.target_stress = 0;
// 	cerr << seq_type << endl;
// 	if (seq_type == "S") {
// 		p.unscaled_contactmodel = true;
// 		cerr << "Repulsive force" << endl;
// 		sys.repulsiveforce = true;
// 		sys.amplitudes.repulsion = 1;
// 		string_control_parameters << "_S" << filename_sequence.substr(0, pos_ext_sequence);
// 	} else if (seq_type == "R") {
// 		//p.unscaled_contactmodel
// 		sys.repulsiveforce = true;
// 		cerr << "Repulsive force" << endl;
// 		cerr << " User Defined Sequence only implemented for ....\n";
// 		exit(1);
// 	} else if (seq_type == "B") {
// 		cerr << "Cohesive force" << endl;
// 		sys.set_shear_rate(1);
// 		sys.repulsiveforce = false;
// 		sys.cohesion = true;
// 		sys.cohesive_force = 1;
// 		string_control_parameters << "_B" << filename_sequence.substr(0, pos_ext_sequence);
// 	} else {
// 		cerr << " User Defined Sequence only implemented for ....\n";
// 		exit(1);
// 	}
// 	setDefaultParameters();
// 	readParameterFile();
// 	if (binary_conf) {
// 		importConfigurationBinary();
// 	} else {
// 		importInitialPositionFile();
// 	}
// 	if (input_files[3] != "not_given") {
// 		importPreSimulationData(input_files[3]);
// 		// time_interval_out
// 		time_interval_output_data = p.time_interval_output_data/shear_rate_expectation;
// 		time_interval_output_config = p.time_interval_output_config/shear_rate_expectation;
// 	} else {
// 		time_interval_output_data = p.time_interval_output_data;
// 		time_interval_output_config = p.time_interval_output_config;
// 	}
// 	if (control_var == "stress") {
// 		if (p.integration_method != 0) {
// 			cerr << "Must be Euler method for stress controlled simulation" << endl;
// 		}
// 		p.integration_method = 0;
// 	}
// 	sys.importParameterSet(p);
// 	sys.setupSystem(control_var);
// 	openOutputFiles(binary_conf);
// 	echoInputFiles(in_args, input_files);
// 	outputConfigurationData();
// 	vector <double> strain_sequence;
// 	vector <double> rsequence;
// 	ifstream fin_seq;
// 	fin_seq.open(filename_sequence.c_str());
// 	if (!fin_seq) {
// 		cerr << " Sequence file '" << filename_sequence << "' not found." <<endl;
// 		exit(1);
// 	}
// 	double strain, targ_st;
// 	while (fin_seq >> strain >> targ_st) {
// 		strain_sequence.push_back(strain);
// 		rsequence.push_back(targ_st);
// 	}
// 	int cnt_simu_loop = 1;
// 	int cnt_config_out = 1;
// 	double next_strain = 0;
// 	double strain_output_config = 0;
// 	double time_output_data = 0;
// 	double time_output_config = 0;
// 	int jammed = 0;
// 	/******************** OUTPUT INITIAL DATA ********************/
// 	evaluateData();
// 	outputRheologyData();
// 	outputStressTensorData();
// 	outputConfigurationBinary();
// 	outputConfigurationData();
// 	/*************************************************************/
// 	for (unsigned int step = 0; step<strain_sequence.size(); step++) {
// 		/* The target stress (``rsequence'') is given trough the command argument
// 		 * with an unit stres: eta_0*gammmadot_0.
// 		 * However, in the code, sys.target_stress is computed as an unit F_rep/a^2.
// 		 */
// 		sys.target_stress_input = rsequence[step];
// 		sys.target_stress = rsequence[step]/6/M_PI;
// 		cerr << "Target stress " << sys.target_stress_input << endl;
// 		sys.updateUnscaledContactmodel();
// 		sys.amplitudes.repulsion = 1; // needed for 1st time step
// 		sys.dimensionless_number = 1;
// 		next_strain = sys.get_shear_strain()+strain_sequence[step];
// 		while (sys.get_shear_strain() < next_strain-1e-8) {
// 			time_output_data = cnt_simu_loop*time_interval_output_data;
// 			time_output_config = cnt_config_out*time_interval_output_config;
// 			sys.timeEvolution(time_output_data);
// 			cnt_simu_loop ++;
			
// 			/******************** OUTPUT DATA ********************/
// 			evaluateData();
// 			outputRheologyData();
// 			outputStressTensorData();
// 			outputConfigurationBinary();
// 			if (time_interval_output_data == -1) {
// 				if (sys.get_shear_strain() >= strain_output_config-1e-8) {
// 					outputConfigurationData();
// 					cnt_config_out ++;
// 				}
// 			} else {
// 				if (sys.get_time() >= time_output_config-1e-8) {
// 					outputConfigurationData();
// 					cnt_config_out ++;
// 				}
// 			}
// 			/******************************************************/
			
// 			if (abs(sys.get_shear_rate()) < p.rest_threshold) {
// 				cerr << "shear jamming " << jammed << endl;
// 				jammed ++;
// 				if (jammed > 10) {
// 					jammed = 0;
// 					cerr << "shear jamming";
// 					break;
// 				}
// 			} else {
// 				jammed = 0;
// 			}
// 			cerr << "strain: " << sys.get_time() << " / " << p.time_end;
// 			cerr << "      stress = " << sys.target_stress_input << endl;
// 		}
// 	}
// }

bool str2bool(const string &value)
{
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		cerr << "The value should be true or false" << endl;
		exit(1);
	}
}

vec3d str2vec3d(const string &value)
{
	string::size_type l1 = value.find("(", 0);
	if (l1 == string::npos) {
		exit(1);
	}
	string::size_type l2 = value.find(",", l1);
	if (l2 == string::npos) {
		exit(1);
	}
	string::size_type l3 = value.find(",", l2+1);
	if (l3 == string::npos) {
		exit(1);
	}
	string::size_type l4 = value.find(")", l3+1);
	if (l4 == string::npos) {
		exit(1);
	}
	double vx = atof(value.substr(l1+1, l2-l1-1).c_str());
	double vy = atof(value.substr(l2+1, l3-l2-1).c_str());
	double vz = atof(value.substr(l3+1, l4-l3-1).c_str());
	return vec3d(vx,vy,vz);
}

void Str2KeyValue(const string &str_parameter,
				  string &keyword,
				  string &value)
{
	string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}


void Simulation::autoSetParameters(const string &keyword, const string &value)
{
	string numeral, suffix;
	if (keyword == "lubrication_model") {
		p.lubrication_model = atoi(value.c_str());
	} else if (keyword == "friction_model") {
		if (p.friction_model == 2) {
			cerr << "!!Neglected friction_model in parameter file!!" << endl;
		} else {
			p.friction_model = atoi(value.c_str());
		}
	} else if (keyword == "rolling_friction") {
		p.rolling_friction = str2bool(value);
	} else if (keyword == "repulsion_amplitude") {
		getSuffix(value, numeral, suffix);
		suffixes["r"] = suffix;
		values["r"] = stof(numeral);
	} else if (keyword == "cohesion_amplitude") {
		getSuffix(value, numeral, suffix);
		suffixes["c"] = suffix;
		values["c"] = stof(numeral);
	} else if (keyword == "brownian_amplitude") {
		getSuffix(value, numeral, suffix);
		suffixes["b"] = suffix;
		values["b"] = stof(numeral);
	} else if (keyword == "critical_load_amplitude") {
		getSuffix(value, numeral, suffix);
		suffixes["cl"] = suffix;
		values["cl"] = stof(numeral);
	} else if (keyword == "magnetic_amplitude") {
		getSuffix(value, numeral, suffix);
		suffixes["m"] = suffix;
		values["m"] = stof(numeral);
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "unscaled_contactmodel") {
		p.unscaled_contactmodel = str2bool(value);
	} else if (keyword == "repulsiveforce_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		p.contact_relaxation_time = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time_tan"){
		p.contact_relaxation_time_tan =  atof(value.c_str());
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "time_end") {
		p.time_end = atof(value.c_str());
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max_gap") {
		p.lub_max_gap = atof(value.c_str());
	} else if (keyword == "interaction_range") {
		p.interaction_range = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		p.kn = atof(value.c_str());
	} else if (keyword == "kt") {
		p.kt = atof(value.c_str());
	} else if (keyword == "kr") {
		p.kr = atof(value.c_str());
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
		p.time_interval_output_config = atof(value.c_str());
	} else if (keyword == "time_interval_output_data") {
		p.time_interval_output_data = atof(value.c_str());
	} else if (keyword == "out_data_particle") {
		p.out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		p.out_data_interaction = str2bool(value);
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
		p.min_kn = atof(value.c_str());
	} else if (keyword == "max_kn") {
		p.max_kn = atof(value.c_str());
	} else if (keyword == "min_kt") {
		p.min_kt = atof(value.c_str());
	} else if (keyword == "max_kt") {
		p.max_kt = atof(value.c_str());
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		getSuffix(value, numeral, suffix);
		suffixes["ft"] = suffix;
		values["ft"] = stof(numeral);
	} else if (keyword == "fixed_dt") {
		p.fixed_dt = str2bool(value);
	} else if (keyword == "ratio_nonmagnetic") {
		p.ratio_nonmagnetic = atof(value.c_str());
	} else if (keyword == "magnetic_dipole_moment") {
		p.magnetic_dipole_moment = atof(value.c_str());
	} else if (keyword == "external_magnetic_field") {
		p.external_magnetic_field = str2vec3d(value);
	} else if (keyword == "dipole_orientation") {
		p.dipole_orientation = atoi(value.c_str());
	} else {
		cerr << "keyword " << keyword << " is not associated with an parameter" << endl;
		exit(1);
	}
}

void Simulation::readParameterFile()
{
	ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		cerr << " Parameter file '" << filename_parameters << "' not found." <<endl;
		exit(1);
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
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if (pos_slashslash != string::npos) {
			cerr << " // is not the syntax to comment out. Use /* comment */" << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void Simulation::openOutputFiles(bool binary_conf)
{
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName(binary_conf);
	string st_filename = "st_" +sys.simu_name + ".dat";
	fout_st.open(st_filename.c_str());
	outputDataHeader(fout_st);
	string rheo_filename = "rheo_" + sys.simu_name + ".dat";
	fout_rheo.open(rheo_filename.c_str());
	string time_filename = "t_" + sys.simu_name + ".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_" + sys.simu_name + ".dat";
	fout_input.open(input_filename.c_str());
	outputDataHeader(fout_rheo);
	//
	string fout_rheo_col_def =
	"#1: shear strain\n"
	"#2: Viscosity\n"
	"#3: N1 viscosity\n"
	"#4: N2 viscosity\n"
	"#5: Viscosity(lub)\n"
	"#6: N1 viscosity(lub)\n"
	"#7: N2 viscosity(lub)\n"
	"#8: Viscosity(xF_contact part)\n"
	"#9: N1 viscosity(xF_contact part)\n"
	"#10: N2 viscosity(xF_contact part)\n"
	"#11: Viscosity(GU_contact part)\n"
	"#12: N1 viscosity(GU_contact part)\n"
	"#13: N2 viscosity(GU_contact part)\n"
	"#14: Viscosity(friction)\n"
	"#15: N1 viscosity(friction)\n"
	"#16: N2 viscosity(friction)\n"
	"#17: Viscosity(repulsive force XF)\n"
	"#18: N1 viscosity(repulsive force XF)\n"
	"#19: N2 viscosity(repulsive force XF)\n"
	"#20: Viscosity(repulsive force GU)\n"
	"#21: N1 viscosity(repulsive force GU)\n"
	"#22: N2 viscosity(repulsive force GU)\n"
	"#23: Viscosity(brownian)\n"
	"#24: N1 viscosity(brownian)\n"
	"#25: N2 viscosity(brownian)\n"
	"#26: particle pressure\n"
	"#27: particle pressure contact\n"
	"#28: min gap (non-dim)\n"
	"#29: max tangential displacement\n"
	"#30: max Fc_normal\n"
	"#31: max Fc_tan\n"
	"#32: max velocity\n"
	"#33: max angular velocity\n"
	"#34: ave contact normal velocity\n"
	"#35: max contact normal velocity\n"
	"#36: ave contact tangential velocity\n"
	"#37: max contact tangential velocity\n"
	"#38: ave sliding velocity\n"
	"#39: max sliding velocity\n"
	"#40: ave contact number per particle\n"
	"#41: num of interaction\n"
	"#42: num of contacts\n"
	"#43: num of frictional contacts\n"
	"#44: kn\n"
	"#45: kt\n"
	"#46: dt\n"
	"#47: time\n"
	"#48: dimensionless_number\n"
	"#49: stress\n"
	"#50: shear_disp\n"
	"#51: max rolling displacement\n"
	"#52: max_contact_gap\n"
	"#53: total_energy\n"
	"#54: magnetic_energy\n";
	
	fout_rheo << fout_rheo_col_def << endl;
	if (p.out_data_particle) {
		string particle_filename = "par_" + sys.simu_name + ".dat";
		fout_particle.open(particle_filename.c_str());
		outputDataHeader(fout_particle);
		//
		string fout_par_col_def =
		"#1: number of the particle\n"
		"#2: radius\n"
		"#3: position x\n"
		"#4: position y\n"
		"#5: position z\n"
		"#6: velocity x\n"
		"#7: velocity y\n"
		"#8: velocity z\n"
		"#9: angular velocity x\n"
		"#10: angular velocity y\n"
		"#11: angular velocity z\n"
		"#12: viscosity contribution of lubrication\n"
		"#13: viscosity contributon of contact GU xz\n"
		"#14: viscosity contributon of brownian xz\n"
		"#15: angle (for 2D simulation only)\n";
		//
		fout_particle << fout_par_col_def << endl;
	}
	if (p.out_data_interaction) {
		string interaction_filename = "int_" + sys.simu_name + ".dat";
		fout_interaction.open(interaction_filename.c_str());
		outputDataHeader(fout_interaction);
		string fout_int_col_def =
		"#1: particle 1 label\n"
		"#2: particle 2 label\n"
		"#3: contact state (0 = no contact, 1 = frictionless contact, 1 = non-sliding frictional, 2 = sliding frictional)\n"
		"#4: normal vector, oriented from particle 1 to particle 2 x\n"
		"#5: normal vector, oriented from particle 1 to particle 2 y\n"
		"#6: normal vector, oriented from particle 1 to particle 2 z\n"
		"#7: dimensionless gap = s-2, s = 2r/(a1+a2)\n"
		"#8: norm of the normal part of the lubrication force\n"
		"#9: tangential part of the lubrication force x\n"
		"#10: tangential part of the lubrication force y\n"
		"#11: tangential part of the lubrication force z\n"
		"#12: norm of the normal part of the contact force\n"
		"#13: tangential part of the contact force, x\n"
		"#14: tangential part of the contact force, y\n"
		"#15: tangential part of the contact force, z\n"
		"#16: norm of the normal repulsive force\n"
		"#17: Viscosity contribution of contact xF\n";
		fout_interaction << fout_int_col_def << endl;
	}
}

void Simulation::setDefaultParameters()
{
	p.brownian_amplitude = 0;
	p.repulsion_amplitude = 0;
	p.cohesion_amplitude = 0;
	p.critical_load_amplitude = 0;
	p.magnetic_amplitude = 0;

	p.Pe_switch = 5;
	p.dt = 1e-4;
	p.disp_max = 2e-3;
	p.monolayer = false;
	p.rest_threshold = 1e-4;
	p.integration_method = 1;
	p.interaction_range = 5;
	/*
	 * Stokes drag coeffient
	 */
	p.sd_coeff = 1;
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 ???
	 */
	p.lubrication_model = 2;
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	p.friction_model = 1;
	p.rolling_friction = false;
	p.time_end = 10;
	p.lub_max_gap = 0.5;
	/*
	 * reduced_gap_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	p.lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxation_factor:
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	p.contact_relaxation_time = 1e-3;
	p.contact_relaxation_time_tan = 0;
	if (control_var == "stress") {
		p.unscaled_contactmodel = true;
		p.kn = 2000;
		p.kt = 1000;
		p.kr = 1000;
	} else {
		p.unscaled_contactmodel = false;
		p.kn = 10000;
		p.kt = 6000;
		p.kr = 6000;
	}
	p.auto_determine_knkt = false;
	p.overlap_target = 0.05;
	p.disp_tan_target = 0.05;
	p.memory_strain_avg = 0.01;
	p.memory_strain_k = 0.02;
	p.start_adjust = 0.2;
	p.min_kn = 1000;
	p.max_kn = 1000000;
	p.min_kt = 1000;
	p.max_kt = 1000000;
	p.repulsive_length = 0.05;
	p.mu_static = 1;
	p.mu_dynamic = -1;
	p.time_interval_output_data = 0.01;
	p.time_interval_output_config = 0.1;
	p.origin_zero_flow = true;
	p.out_data_particle = true;
	p.out_data_interaction = true;
	p.ft_max = 1;
	p.fixed_dt = false;
	p.magnetic_dipole_moment = 1;
	p.ratio_nonmagnetic = 0;
	p.dipole_orientation = 0;
}

void Simulation::importInitialPositionFile()
{
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	bool include_magnetic_moment = false;
	if (filename_import_positions.find("mag", 0) != string::npos) {
		cerr << "The initial configuration file includes magnetic moment." << endl;
		include_magnetic_moment = true;
	}
	char buf;
	int n1, n2;
	double lx, ly, lz, vf1, vf2;
	getline(file_import, import_line[0]);
	getline(file_import, import_line[1]);
	stringstream ss(import_line[1]);
	ss >> buf >> n1 >> n2 >> volume_or_area_fraction >> lx >> ly >> lz >> vf1 >> vf2 >> initial_lees_edwards_disp;
	double x_, y_, z_, a_;
	vector<vec3d> initial_position;
	vector <double> radius;
	if (include_magnetic_moment == false) {
		while (file_import >> x_ >> y_ >> z_ >> a_) {
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
		}
	} else {
		double mx_, my_, mz_;
		while (file_import >> x_ >> y_ >> z_ >> a_ >> mx_ >> my_ >> mz_ ) {
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
			sys.init_magnetic_moment.push_back(vec3d(mx_, my_, mz_));
		}
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void Simulation::outputConfigurationBinary()
{
	string conf_filename;
	//	conf_filename =  "conf_" + sys.simu_name + "_strain" + to_string(sys.get_shear_strain()) + ".dat";
	conf_filename =  "conf_" + sys.simu_name + ".dat";
	outputConfigurationBinary(conf_filename);
}

void Simulation::outputConfigurationBinary(string conf_filename)
{
	vector < vector <double> > pos;
	int np = sys.get_np();
	int dims = 4;
	pos.resize(np);
	for (int i=0; i<np; i++) {
		pos[i].resize(dims);
		pos[i][0] = sys.position[i].x;
		pos[i][1] = sys.position[i].y;
		pos[i][2] = sys.position[i].z;
		pos[i][3] = sys.radius[i];
	}
	ofstream conf_export;
	double lx = sys.get_lx();
	double ly = sys.get_ly();
	double lz = sys.get_lz();
	double shear_disp = sys.shear_disp;
	conf_export.open(conf_filename.c_str(), ios::binary | ios::out);
	conf_export.write((char*)&np, sizeof(int));
	conf_export.write((char*)&volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&lx, sizeof(double));
	conf_export.write((char*)&ly, sizeof(double));
	conf_export.write((char*)&lz, sizeof(double));
	conf_export.write((char*)&shear_disp, sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	conf_export.close();
}

void Simulation::importConfigurationBinary()
{
	ifstream file_import;
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	int np;
	double lx, ly, lz;
	file_import.read((char*)&np, sizeof(int));
	file_import.read((char*)&volume_or_area_fraction, sizeof(double));
	file_import.read((char*)&lx, sizeof(double));
	file_import.read((char*)&ly, sizeof(double));
	file_import.read((char*)&lz, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp, sizeof(double));
	double x_, y_, z_, r_;
	vector <vec3d> initial_position;
	vector <double> radius;
	for (int i=0; i<np; i++) {
		file_import.read((char*)&x_, sizeof(double));
		file_import.read((char*)&y_, sizeof(double));
		file_import.read((char*)&z_, sizeof(double));
		file_import.read((char*)&r_, sizeof(double));
		initial_position.push_back(vec3d(x_,y_,z_));
		radius.push_back(r_);
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void Simulation::prepareSimulationName(bool binary_conf)
{
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
	ss_simu_name << string_control_parameters.str();
	sys.simu_name = ss_simu_name.str();
	cerr << "filename: " << sys.simu_name << endl;
}

void Simulation::evaluateData()
{
	/**
	 \brief Get rheological data from the System class.
	 
	 Data are converted in hydrodynamic units, independently from the actual units used in the System class.
	 */
	
	sys.analyzeState();
	sys.calcStress();
	sys.calcLubricationForce();
	
	double stress_unit_converter = 0;
	if (unit_scales == "hydro") {
		stress_unit_converter = 1;
	}
	if (unit_scales == "thermal") {
		stress_unit_converter = 1/dimensionless_numbers["b"];
	}
	if (unit_scales == "repulsion") {
		stress_unit_converter = 1/dimensionless_numbers["r"];
	}
	
	viscosity = stress_unit_converter*(sys.einstein_stress+sys.total_stress.getStressXZ());
	normalstress_diff_1 = stress_unit_converter*sys.total_stress.getNormalStress1();
	normalstress_diff_2 = stress_unit_converter*sys.total_stress.getNormalStress2();
	particle_pressure = stress_unit_converter*sys.total_stress.getParticlePressure();
	viscosity_hydro = stress_unit_converter*sys.total_hydro_stress.getStressXZ();
	normalstress_diff_1_hydro = stress_unit_converter*sys.total_hydro_stress.getNormalStress1();
	normalstress_diff_2_hydro = stress_unit_converter*sys.total_hydro_stress.getNormalStress2();
	viscosity_cont_XF = stress_unit_converter*sys.total_contact_stressXF.getStressXZ();
	normalstress_diff_1_cont_XF = stress_unit_converter*sys.total_contact_stressXF.getNormalStress1();
	normalstress_diff_2_cont_XF = stress_unit_converter*sys.total_contact_stressXF.getNormalStress2();
	particle_pressure_cont = stress_unit_converter*sys.total_contact_stressXF.getParticlePressure();
	viscosity_friction = stress_unit_converter*sys.total_contact_stressXF_tan.getStressXZ();
	normalstress_diff_1_friction = stress_unit_converter*sys.total_contact_stressXF_tan.getNormalStress1();
	normalstress_diff_2_friction = stress_unit_converter*sys.total_contact_stressXF_tan.getNormalStress2();
	viscosity_cont_GU = stress_unit_converter*sys.total_contact_stressGU.getStressXZ();
	normalstress_diff_1_cont_GU = stress_unit_converter*sys.total_contact_stressGU.getNormalStress1();
	normalstress_diff_2_cont_GU = stress_unit_converter*sys.total_contact_stressGU.getNormalStress2();
	if (sys.repulsiveforce) {
		viscosity_repulsive_XF = stress_unit_converter*sys.total_repulsive_stressXF.getStressXZ();
		normalstress_diff_1_repulsive_XF = stress_unit_converter*sys.total_repulsive_stressXF.getNormalStress1();
		normalstress_diff_2_repulsive_XF = stress_unit_converter*sys.total_repulsive_stressXF.getNormalStress2();
		particle_pressure_repulsive = stress_unit_converter*sys.total_repulsive_stressXF.getParticlePressure();
		viscosity_repulsive_GU = stress_unit_converter*sys.total_repulsive_stressGU.getStressXZ();
		normalstress_diff_1_repulsive_GU = stress_unit_converter*sys.total_repulsive_stressGU.getNormalStress1();
		normalstress_diff_2_repulsive_GU = stress_unit_converter*sys.total_repulsive_stressGU.getNormalStress2();
	} else {
		viscosity_repulsive_XF = 0;
		normalstress_diff_1_repulsive_XF = 0;
		normalstress_diff_2_repulsive_XF = 0;
		particle_pressure_repulsive = 0;
		viscosity_repulsive_GU = 0;
		normalstress_diff_1_repulsive_GU = 0;
		normalstress_diff_2_repulsive_GU = 0;
	}
	if (sys.brownian) {
		viscosity_brownian = stress_unit_converter*sys.total_brownian_stressGU.getStressXZ();
		normalstress_diff_1_brownian = stress_unit_converter*sys.total_brownian_stressGU.getNormalStress1();
		normalstress_diff_2_brownian = stress_unit_converter*sys.total_brownian_stressGU.getNormalStress2();
	} else {
		viscosity_brownian = 0;
		normalstress_diff_1_brownian = 0;
		normalstress_diff_2_brownian = 0;
	}
}

void Simulation::outputStressTensorData()
{
	fout_st << sys.get_shear_strain() << ' ';
	fout_st << 6*M_PI*viscosity << ' ';
	/* total_stress = sys.total_hydro_stress;
	 * + total_contact_stressXF + total_repulsive_stress;
	 */
	// As it is, the output stress lacks a 6pi factor (as the viscosity)
	sys.total_stress.outputStressTensor(fout_st); // (3,4,5,6,7,8)
	sys.total_hydro_stress.outputStressTensor(fout_st); // (9,10,11,12,13,14)
	sys.total_contact_stressXF.outputStressTensor(fout_st); // (15,16,17,18,19,20)
	sys.total_contact_stressGU.outputStressTensor(fout_st); // (21,22,23,24,25,26)
	sys.total_repulsive_stress.outputStressTensor(fout_st); // (27,28,29,30,31,32)
	sys.total_brownian_stressGU.outputStressTensor(fout_st); // (33,34,35,36,37,38)
	//	fout_st << sys.dimensionless_number << ' '; // 39
	fout_st << endl;
}

void Simulation::outputRheologyData()
{
	/**
	 \brief Output rheological data.
	 
	 
	 Stress data are converted in units of
	 \f$\eta_0\dot\gamma\f$. Other data are output in the units used
	 in the System class (these can be hydrodynamic, Brownian or
	 repulsive force units).
	 
	 \b NOTE: this behavior should be changed
	 and made more consistent in the future.
	 */
	
	
	
	/*
	 * Output the sum of the normal forces.
	 *
	 *  Viscosity = S_{xz} / shear_rate
	 *  N1 = S_{xx}-S_{zz}
	 *  N2 = S_{zz}-S_{yy} = S_zz-(-S_xx-S_zz) = S_xx+2*S_zz
	 *
	 * Relative viscosity = Viscosity / viscosity_solvent
	 */
	
	/*
	 * hat(...) indicates dimensionless quantities.
	 * (1) relative viscosity = Sxz/(eta0*shear_rate) = 6*pi*hat(Sxz)
	 * (2) N1/(eta0*shear_rate) = 6*pi*hat(N1)
	 * (3) N2/(eta0*shear_rate) = 6*pi*hat(N2)
	 *
	 * In simulation, we use the force unit where Stokes drag is F = -(U-U^inf)
	 *
	 * [note] In stress controlled simulation,
	 * Averaged viscosity need to be calculated with dimensionless_number_averaged,
	 * i.e. <viscosity> = taget_stress / dimensionless_number_averaged.
	 */
	fout_rheo << sys.get_shear_strain() << ' '; //1
	fout_rheo << 6*M_PI*viscosity << ' '; //2
	fout_rheo << 6*M_PI*normalstress_diff_1 << ' '; //3
	fout_rheo << 6*M_PI*normalstress_diff_2 << ' '; //4
	/*
	 * Hydrodynamic contribution means
	 * stresslet_hydro_GU_i+stresslet_ME_i from vel_hydro
	 * vel_hydro is obtained with GE for the rhs.
	 *
	 * "_hydro" might be bit confusing.
	 * Something indicating "E_inf" would be better.
	 */
	fout_rheo << 6*M_PI*viscosity_hydro << ' '; //5
	fout_rheo << 6*M_PI*normalstress_diff_1_hydro << ' '; //6
	fout_rheo << 6*M_PI*normalstress_diff_2_hydro << ' '; //7
	/*
	 * Contact force contribution seems to be
	 * the sum of viscosity_cont_XF and viscosity_cont_GU.
	 */
	fout_rheo << 6*M_PI*viscosity_cont_XF << ' '; //8
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_XF << ' '; //9
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_XF << ' '; //10
	fout_rheo << 6*M_PI*viscosity_cont_GU << ' ' ; //11
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_GU << ' ' ; //12
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_GU << ' ' ; //13
	/*
	 *
	 */
	fout_rheo << 6*M_PI*viscosity_friction << ' '; //14
	fout_rheo << 6*M_PI*normalstress_diff_1_friction << ' '; //15
	fout_rheo << 6*M_PI*normalstress_diff_2_friction  << ' '; //16
	fout_rheo << 6*M_PI*viscosity_repulsive_XF << ' '; //17
	fout_rheo << 6*M_PI*normalstress_diff_1_repulsive_XF << ' '; //18
	fout_rheo << 6*M_PI*normalstress_diff_2_repulsive_XF << ' '; //19
	fout_rheo << 6*M_PI*viscosity_repulsive_GU << ' '; //20
	fout_rheo << 6*M_PI*normalstress_diff_1_repulsive_GU << ' '; //21
	fout_rheo << 6*M_PI*normalstress_diff_2_repulsive_GU << ' '; //22
	fout_rheo << 6*M_PI*viscosity_brownian << ' ' ; //23
	fout_rheo << 6*M_PI*normalstress_diff_1_brownian << ' ' ; //24
	fout_rheo << 6*M_PI*normalstress_diff_2_brownian << ' ' ; //25
	fout_rheo << 6*M_PI*particle_pressure << ' ';//26
	fout_rheo << 6*M_PI*particle_pressure_cont << ' ';//27
	fout_rheo << sys.min_reduced_gap << ' '; //28
	fout_rheo << sys.max_disp_tan << ' '; //29
	fout_rheo << sys.max_fc_normal << ' '; //30
	fout_rheo << sys.max_fc_tan << ' ';//31
	fout_rheo << sys.max_velocity << ' '; //32
	fout_rheo << sys.max_ang_velocity << ' '; //33
	fout_rheo << sys.ave_contact_velo_normal << ' '; //34
	fout_rheo << sys.max_contact_velo_normal << ' '; //35
	fout_rheo << sys.ave_contact_velo_tan << ' '; //36
	fout_rheo << sys.max_contact_velo_tan << ' '; //37
	fout_rheo << sys.ave_sliding_velocity << ' ' ; //38
	fout_rheo << sys.max_sliding_velocity << ' ' ; //39
	fout_rheo << sys.getParticleContactNumber() << ' ';//40
	fout_rheo << sys.get_nb_of_active_interactions() << ' ';//41
	fout_rheo << sys.contact_nb << ' '; //42
	fout_rheo << sys.fric_contact_nb << ' '; //43
	fout_rheo << sys.kn << ' '; //44
	fout_rheo << sys.kt << ' '; //45
	fout_rheo << sys.dt << ' '; //46
	fout_rheo << sys.get_time() << ' ' ; //47
	/* In stress control simulation,
	 * shear jammed state may cause oscilation of dimensionless_number around 0.
	 * Then, time step also oscilate.
	 * This is why we need to take time average to have correct value of dimensionless_number.
	 */
	//	fout_rheo << sys.dimensionless_number << ' '; // 48
	// if (control_var == "stress") {
	// 	fout_rheo << sys.target_stress_input << ' '; // 49
	// } else {
	// 	fout_rheo << 6*M_PI*viscosity*sys.dimensionless_number << ' '; // 49
	// }
	fout_rheo << sys.shear_disp << ' '; // 50
	fout_rheo << sys.max_disp_rolling << ' '; //51
	fout_rheo << sys.max_contact_gap << ' '; //52
	fout_rheo << sys.get_total_energy() << ' '; // 53;
	fout_rheo << sys.get_magnetic_energy() << ' ';// 54;
	fout_rheo << endl;
}

vec3d Simulation::shiftUpCoordinate(double x, double y, double z)
{
	if (p.origin_zero_flow) {
		z += sys.Lz_half();
		if (z > sys.Lz_half()) {
			x -= sys.shear_disp;
			if (x < -sys.Lx_half()) {
				x += sys.get_lx();
			}
			z -= sys.get_lz();
		}
	}
	return vec3d(x,y,z);
}

void Simulation::outputDataHeader(ofstream &fout)
{
	fout << "# LF_DEM version " << GIT_VERSION << endl;
	fout << "# np " << sys.get_np() << endl;
	fout << "# VF " << sys.volume_fraction << endl;
	fout << "# Lx " << sys.get_lx() << endl;
	fout << "# Ly " << sys.get_ly() << endl;
	fout << "# Lz " << sys.get_lz() << endl;
}

void Simulation::outputConfigurationData()
{
	vector<vec3d> pos;
	vector<vec3d> vel;
	int np = sys.get_np();
	pos.resize(np);
	vel.resize(np);
	for (int i=0; i<np; i++) {
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.Lx_half(),
								   sys.position[i].y-sys.Ly_half(),
								   sys.position[i].z-sys.Lz_half());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	if (p.origin_zero_flow) {
		for (int i=0; i<np; i++) {
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0) {
				vel[i].x -= sys.get_shear_rate()*sys.get_lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (p.out_data_particle) {
		cerr << "   out config: " << sys.get_shear_strain() << endl;
		
		fout_particle << "# " << sys.get_shear_strain() << ' ';
		fout_particle << sys.shear_disp << ' ';
		//		fout_particle << sys.dimensionless_number << ' ';
		fout_particle << sys.target_stress_input << ' ';
		fout_particle << sys.get_time() << endl;
		
		for (int i=0; i<np; i++) {
			const vec3d &r = pos[i];
			const vec3d &v = vel[i];
			const vec3d &o = sys.ang_velocity[i];
			double lub_xzstress = sys.lubstress[i].getStressXZ();
			double contact_xzstressGU = sys.contactstressGU[i].getStressXZ();
			double brownian_xzstressGU = 0;
			if (sys.brownian) {
				brownian_xzstressGU = sys.brownianstressGU[i].getStressXZ();
			}
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << ' ' << r.x << ' ' << r.y << ' ' << r.z; //3, 4, 5: position
			fout_particle << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
			fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions
			fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions
			fout_particle << ' ' << 6*M_PI*brownian_xzstressGU; //14: xz stress contributions
			if (sys.magnetic) {
				fout_particle << ' ' << sys.magnetic_moment[i].x;
				fout_particle << ' ' << sys.magnetic_moment[i].y;
				fout_particle << ' ' << sys.magnetic_moment[i].z;
			} else {
				if (sys.twodimension) {
					fout_particle << ' ' << sys.angle[i]; // 15
				}
			}
			fout_particle << endl;
		}
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			cnt_interaction ++;
		}
	}
	if (p.out_data_interaction) {
		fout_interaction << "# " << sys.get_shear_strain();
		fout_interaction << ' ' << cnt_interaction;
		fout_interaction << ' ' << sys.get_time();
		fout_interaction << endl;
		for (int k=0; k<sys.nb_interaction; k++) {
			if (sys.interaction[k].is_active()) {
				unsigned short i, j;
				sys.interaction[k].get_par_num(i, j);
				vec3d& nr_vec = sys.interaction[k].nvec;
				StressTensor stress_contact = sys.interaction[k].contact.getContactStressXF();
				fout_interaction << i << ' ' << j << ' '; // 1, 2
				/* contact.state:
				 * 0 no contact
				 * 1 Friction is not activated (critical load model)
				 * 2 Static friction
				 * 3 Sliding
				 */
				fout_interaction << sys.interaction[k].contact.state << ' '; //3
				fout_interaction << nr_vec.x << ' '; // 4
				fout_interaction << nr_vec.y << ' '; // 5
				fout_interaction << nr_vec.z << ' '; // 6
				fout_interaction << sys.interaction[k].get_reduced_gap() << ' '; // 7
				/* [NOTE]
				 * Lubrication forces are reference values
				 * in the Brownian case. The force balancing
				 * velocities are recalculated without
				 * including the Brownian forces.
				 * It seems there is no better way to visualize
				 * the lubrication forces.
				 */
				fout_interaction << sys.interaction[k].lubrication.get_lubforce_normal() << ' '; // 8
				fout_interaction << sys.interaction[k].lubrication.get_lubforce_tan() << ' '; // 9, 10, 11
				/*
				 * Contact forces include only spring forces.
				 */
				fout_interaction << sys.interaction[k].contact.get_f_contact_normal_norm() << ' '; // 12
				fout_interaction << sys.interaction[k].contact.get_f_contact_tan() << ' '; // 13, 14, 15
				fout_interaction << sys.interaction[k].repulsion.getForceNorm() << ' '; // 16
				fout_interaction << 6*M_PI*stress_contact.getStressXZ() << ' '; // 17
				sys.interaction[k].contact.addUpContactForceTorque();
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << ' ';
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << ' ';
				fout_interaction << endl;
			}
		}
	}
}

void Simulation::outputFinalConfiguration()
{
	ofstream fout_finalconfig;
	string filename_final_configuration = "./after_relax/"+filename_import_positions;
	fout_finalconfig.open(filename_final_configuration.c_str());
	fout_finalconfig << import_line[0] << endl;
	fout_finalconfig << import_line[1] << endl;
	int np = sys.get_np();
	for (int i=0; i<np; i++) {
		fout_finalconfig << sys.position[i].x << ' ';
		fout_finalconfig << sys.position[i].y << ' ';
		fout_finalconfig << sys.position[i].z << ' ';
		fout_finalconfig << sys.radius[i] << endl;
	}
	string filename_bin = filename_final_configuration;
	string ext=".dat";
	size_t start_pos = filename_bin.find(ext);
	if (start_pos == string::npos) {
		cerr << " WARNING, no binary output generated " << endl;
		return;
	}
	filename_bin.replace(start_pos, ext.length(), ".bin");
	outputConfigurationBinary(filename_bin);
}
