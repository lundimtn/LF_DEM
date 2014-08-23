//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#define VERSION "3.0"
#include "Simulation.h"
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype>

Simulation::Simulation(){};

Simulation::~Simulation(){
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

void
Simulation::contactForceParameter(string filename){
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	double phi_;
	double kn_;
	double kt_;
	double dt_max_;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_max_) {
		if (phi_ == volume_fraction) {
			break;
		}
	}
	fin_knktdt.clear();
	sys.kn = kn_;
	sys.kt = kt_;
	sys.dt_max = dt_max_;
	cerr << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_max_ << endl;
}


void
Simulation::setupSimulationSteadyShear(vector<string> &input_files,
							double peclet_num, double scaled_repulsion,
							double scaled_cohesion,
							double scaled_critical_load,
							string control_variable){
	
	control_var = control_variable;
	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];
	if (control_var == "strain") {
		if (scaled_repulsion > 0 &&
			scaled_critical_load > 0) {
			cerr << " Repulsion AND Critical Load cannot be used at the same time" << endl;
			exit(1);
		}
		if (peclet_num > 0) {
			cerr << "Brownian" << endl;
			sys.brownian = true;
			sys.dimensionless_shear_rate = peclet_num;
			if (scaled_repulsion > 0) {
				cerr << "Repulsive force" << endl;
				sys.repulsiveforce_amplitude = scaled_repulsion/peclet_num;
				sys.repulsiveforce = true;
			}
			if (scaled_critical_load > 0) {
				cerr << "Critical load" << endl;
				sys.critical_normal_force = scaled_critical_load/peclet_num;
				sys.friction_model = 2;
			}
		} else {
			cerr << "non-Brownian" << endl;
			if (scaled_repulsion > 0 && scaled_cohesion == 0) {
				cerr << "Repulsive force" << endl;
				sys.dimensionless_shear_rate = 1/scaled_repulsion;
				sys.repulsiveforce_amplitude = scaled_repulsion;
				sys.repulsiveforce = true;
			}
			if (scaled_critical_load > 0 && scaled_cohesion == 0) {
				sys.dimensionless_shear_rate = 1/scaled_critical_load;
				sys.friction_model = 2;
				cerr << "Critical load" << endl;
				sys.critical_normal_force = scaled_critical_load;
			}
			if (scaled_repulsion == 0 && scaled_cohesion > 0) {
				cerr << "Cohesive force" << endl;
				sys.dimensionless_shear_rate = 1/scaled_cohesion;
				sys.cohesive_force = scaled_cohesion;
			}
			if (scaled_repulsion > 0 && scaled_cohesion > 0) {
				cerr << "Repulsive force + Cohesive force" << endl;
				sys.dimensionless_shear_rate = 1/scaled_repulsion;
				sys.repulsiveforce_amplitude = 1/sys.dimensionless_shear_rate;
				sys.repulsiveforce = true;
				sys.cohesive_force = scaled_cohesion;
			}
		}
	} else if (control_var == "stress") {
		sys.brownian = false;
		if (scaled_critical_load > 0) {
			cerr << " Stress controlled simulations for CLM not implemented ! " << endl;
			exit(1);
		}
		if (scaled_repulsion == 0) {
			cerr << " Stress controlled simulations need a repulsive force ! " << endl;
			exit(1);
		} else {
			sys.repulsiveforce = true;
			sys.repulsiveforce_amplitude = 1;
			sys.target_stress = 1/scaled_repulsion;
			sys.dimensionless_shear_rate = 1; // needed for 1st time step
		}
	}
	setDefaultParameters();
	readParameterFile();
	if (sys.cohesive_force > 0) {
		sys.cohesive_force = sys.cohesive_force/sys.dimensionless_shear_rate; // < Why is that here? It seems it gives sys.cohesive_force = scaled_cohesion^2 for pure cohesion. Is this what is intended?
	}
	importInitialPositionFile();

	int fnb = input_files.size();
	if (fnb == 3) {
		contactForceParameter(input_files[2]);
	}
	openOutputFiles();
	if (sys.brownian) {
		sys.setupBrownian();
	}
	sys.setupSystem(control_var);
	if (filename_parameters == "init_relax.txt") {
		sys.zero_shear = true;
	}
	outputConfigurationData();
	sys.setupShearFlow(true);

	if(control_var == "stress"){
		sys.set_integration_method(0);
	}


}

/*
 * Main simulation
 */
void
Simulation::simulationSteadyShear(vector<string> &input_files,
								  double peclet_num, double scaled_repulsion, double scaled_cohesion,
								  double scaled_critical_load, string control_variable){
	user_sequence = false;
	control_var = control_variable;

	setupSimulationSteadyShear(input_files, peclet_num,
							   scaled_repulsion, scaled_cohesion, scaled_critical_load, control_var);

	int cnt_simu_loop = 1;
	//int cnt_knkt_adjustment = 1;
	int cnt_config_out = 1;
	while (sys.get_shear_strain() < sys.shear_strain_end-1e-8) {
		//double strain_knkt_adjustment = cnt_knkt_adjustment*strain_interval_knkt_adjustment;
		double strain_next_config_out = cnt_config_out*sys.strain_interval_output;
		double strain_next = cnt_simu_loop*sys.strain_interval_output_data;
		sys.timeEvolution(strain_next);
		evaluateData();
		outputRheologyData();
		outputStressTensorData();
		if (sys.get_shear_strain() >= strain_next_config_out-1e-8) {
			outputConfigurationData();
			cnt_config_out ++;
		}
		//if (sys.kn_kt_adjustment) {
		//			if (sys.get_shear_strain() >= strain_knkt_adjustment-1e-8) {
		//				if (sys.adjustContactModelParameters() == 1){
		//					cout << "phi kn kt dt" << endl;
		//					cout << volume_fraction << ' ';
		//					cout << sys.get_kn() << ' ' ;
		//					cout << sys.get_kt() << ' ';
		//					cout << sys.get_dt() << endl;
		//					if (sys.get_kn() > sys.max_kn){
		//						cout << "kn cannot be determined. It can be larger than the upper limit." << endl;
		//					}
		//					return;
		//				}
		//				cnt_knkt_adjustment ++;
		//			}
		//}
		cnt_simu_loop ++;
		cerr << "strain: " << sys.get_shear_strain() << " / " << sys.shear_strain_end << endl;
	}
	if (filename_parameters == "init_relax.txt") {
		/* To prepar relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to expoert the position data.
		 */
		outputFinalConfiguration();
	}
}


/*
 * Main simulation
 */
void
Simulation::simulationUserDefinedSequence(string seq_type, vector<string> &input_files, string control_variable){
	user_sequence = true;
	control_var = control_variable;
	if (seq_type != "r" || control_var != "stress") {
		cerr << " User Defined Sequence only implemented for pure repulsive force under stress control at the moment "; exit(1);
	}
	
	cout << " User Defined Sequence, Repulsive Force " << endl;
	sys.repulsiveforce = true;
	sys.repulsiveforce_amplitude = 1;
	sys.dimensionless_shear_rate = 1; // needed for 1st time step

	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];

	setDefaultParameters();
	readParameterFile();
	
	if (control_var == "stress") {
		sys.set_integration_method(0);
	}

	importInitialPositionFile();

	int fnb = input_files.size();
	if (fnb == 4) {
		contactForceParameter(input_files[2]);
	}
	filename_sequence = input_files[fnb-1];
	
	openOutputFiles();

	sys.setupSystem(control_var);
	outputConfigurationData();

	sys.setupShearFlow(true);

	vector <double> strain_sequence;
	vector <double> rsequence;

	ifstream fin_seq;
	fin_seq.open(filename_sequence.c_str());

	double strain;
	double targ_st;

	while (fin_seq >> strain >> targ_st){
		strain_sequence.push_back(strain);
		rsequence.push_back(targ_st);
	}
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	double next_strain=0;
	for(unsigned int step=0; step<strain_sequence.size();step++){
		sys.target_stress = 1/rsequence[step];
		next_strain += strain_sequence[step];
		while (sys.get_shear_strain() < next_strain-1e-8) {
			double strain_next_config_out = cnt_config_out*sys.strain_interval_output;
			double strain_next = cnt_simu_loop*sys.strain_interval_output_data;
			sys.timeEvolution(strain_next);
			evaluateData();
			outputRheologyData();
			outputStressTensorData();
			if (sys.get_shear_strain() >= strain_next_config_out-1e-8) {
				outputConfigurationData();
				cnt_config_out ++;
			}
			cnt_simu_loop ++;
			cerr << "strain: " << sys.get_shear_strain() << " / " << sys.shear_strain_end << endl;
		}
	}
}

bool
str2bool(string value){
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		cerr << "The value should be true or false" << endl;
		exit(1);
	}
}

void
Str2KeyValue(string &str_parameter,
			 string &keyword,
			 string &value){
	string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

void
removeBlank(string &str){
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

void
Simulation::autoSetParameters(const string &keyword, const string &value){
	if (keyword == "lubrication_model") {
		sys.set_lubrication_model(atoi(value.c_str()));
	} else if (keyword == "friction_model") {
		if (sys.friction_model == 2) {
			cerr << "!!Neglected friction_model in parameter file!!" << endl;
		} else {
			sys.friction_model = atoi(value.c_str());
		}
	} else if (keyword == "rolling_friction") {
		sys.rolling_friction = str2bool(value);
	} else if (keyword == "kn_kt_adjustment") {
		sys.kn_kt_adjustment = str2bool(value);
	} else if (keyword == "strain_interval_knkt_adjustment") {
		strain_interval_knkt_adjustment = atof(value.c_str());
	} else if (keyword == "repulsiveforce_length") {
		if (sys.repulsiveforce) {
			sys.set_repulsiveforce_length(atof(value.c_str()));
		}
	} else if (keyword == "lub_reduce_parameter") {
		sys.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		sys.contact_relaxation_time = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time_tan"){
		sys.contact_relaxation_time_tan =  atof(value.c_str());
	} else if (keyword == "disp_max") {
		sys.set_disp_max(atof(value.c_str()));
	} else if (keyword == "shear_strain_end") {
		sys.shear_strain_end = atof(value.c_str());
	} else if (keyword == "integration_method") {
		sys.set_integration_method(atoi(value.c_str()));
	} else if (keyword == "lub_max") {
		sys.set_lub_max(atof(value.c_str()));
	} else if (keyword == "sd_coeff") {
		sys.set_sd_coeff(atof(value.c_str()));
	} else if (keyword == "kn") {
		sys.kn = atof(value.c_str());
	} else if (keyword == "kt") {
		sys.kt = atof(value.c_str());
	} else if (keyword == "kr") {
		sys.kr = atof(value.c_str());
	} else if (keyword == "dt_max") {
		sys.dt_max = atof(value.c_str());
	} else if (keyword == "kn_lowPeclet") {
		sys.kn_lowPeclet = atof(value.c_str());
	} else if (keyword == "kt_lowPeclet") {
		sys.kt_lowPeclet = atof(value.c_str());
	} else if (keyword == "dt_lowPeclet") {
		sys.dt_lowPeclet = atof(value.c_str());
	} else if (keyword == "Pe_switch") {
		sys.Pe_switch = atof(value.c_str());
	} else if (keyword == "mu_static") {
		sys.set_mu_static(atof(value.c_str()));
	} else if (keyword == "strain_interval_out") {
		sys.strain_interval_output = atof(value.c_str());
	} else if (keyword == "strain_interval_out_data") {
		sys.strain_interval_output_data = atof(value.c_str());
	} else if (keyword == "out_data_particle") {
		out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		out_data_interaction = str2bool(value);
	} else if (keyword == "origin_zero_flow") {
		origin_zero_flow = str2bool(value);
	} else if (keyword == "overlap_target") {
		sys.overlap_target = atof(value.c_str());
	} else if (keyword == "disp_tan_target") {
		sys.disp_tan_target = atof(value.c_str());
	} else {
		cerr << "keyword " << keyword << " is not associated with an parameter" << endl;
		exit(1);
	}
}

void
Simulation::readParameterFile(){
	ifstream fin;
	fin.open(filename_parameters.c_str());
	string keyword;
	string value;
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
		if (begin_comment > end_comment ) {
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if( pos_slashslash != string::npos) {
			cerr << " // is not syntax for comment out." << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void
Simulation::openOutputFiles(){
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName();
	string particle_filename = "par_" + sys.simu_name + ".dat";
	string interaction_filename = "int_" + sys.simu_name + ".dat";
	string vel_filename = "rheo_" + sys.simu_name + ".dat";
	string st_filename = "st_" +sys.simu_name + ".dat";
	fout_particle.open(particle_filename.c_str());
	fout_interaction.open(interaction_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	fout_st.open(st_filename.c_str());
	outputDataHeader(fout_particle);
	outputDataHeader(fout_interaction);
	outputDataHeader(fout_rheo);
	outputDataHeader(fout_st);
	//
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
	//
	string fout_rheo_col_def =
	"#1: shear strain\n"
	"#2: Viscosity\n"
	"#3: N1\n"
	"#4: N2\n"
	"#5: Viscosity(lub)\n"
	"#6: N1(lub)\n"
	"#7: N2(lub)\n"
	"#8: Viscosity(xF_contact part)\n"
	"#9: N1(xF_contact part)\n"
	"#10: N2(xF_contact part)\n"
	"#11: Viscosity(GU_contact part)\n"
	"#12: N1(GU_contact part)\n"
	"#13: N2(GU_contact part)\n"
	"#14: Viscosity(friction)\n"
	"#15: N1(friction)\n"
	"#16: N2(friction)\n"
	"#17: Viscosity(repulsive force XF)\n"
	"#18: N1(repulsive force XF)\n"
	"#19: N2(repulsive force XF)\n"
	"#20: Viscosity(repulsive force GU)\n"
	"#21: N1(repulsive force GU)\n"
	"#22: N2(repulsive force GU)\n"
	"#23: Viscosity(brownian)\n"
	"#24: N1(brownian)\n"
	"#25: N2(brownian)\n"
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
	"#47: time\n";
	//
	fout_rheo << fout_rheo_col_def << endl;
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

void
Simulation::setDefaultParameters(){
	/*
	 * Simulation
	 *
	 * dt: the time step to integrate the equation of motion.
	 *     We need to give a good criterion to give.
	 * dt_mid: the intermediate time step for the mid-point
	 *     algortithm. dt/dt_mid = dt_ratio
	 *     Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 *    ASD code from Brady has dt_ratio=150
	 *
	 */
	sys.Pe_switch = 5;
	sys.dt_max = 1e-4;
	sys.dt_lowPeclet = 1e-4;
	/*
	 * integration_method:
	 * 0 Euler's Method,
	 * 1 predictor-corrector,
	 */
	int _integration_method = 1;
	/*
	 * Stokes drag coeffient
	 */
	double _sd_coeff = 1;
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 ???
	 */
	int _lubrication_model = 2;
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	if (sys.friction_model != 2) {
		sys.friction_model = 1;
	}
	sys.rolling_friction = false;
	/*
	 * Shear flow
	 *  shear_rate: shear rate
	 *  strain(): total strain (length of simulation)
	 *
	 */
	sys.shear_strain_end = 10;
	/*
	 * Lubrication force
	 * lub_max: reduced large cutoff distance for lubrication
	 * I think lub_max = 2.5 and 3 generate different results.
	 * We should give suffiently larger value.
	 * The value 3 or 3.5 should be better (To be checked.)
	 */
	double _lub_max = 2.5;
	/*
	 * gap_nondim_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	sys.lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxation_factor:
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	sys.contact_relaxation_time = 1e-3;
	sys.contact_relaxation_time_tan = 0;
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	sys.kn = 10000;
	sys.kt = 6000;
	sys.kr = 6000;
	sys.kn_lowPeclet = 10000;
	sys.kt_lowPeclet = 6000;
	sys.kr_lowPeclet = 6000;
	sys.kn_kt_adjustment = false;
	strain_interval_knkt_adjustment = 5;
	sys.overlap_target = 0.05;
	sys.disp_tan_target = 0.05;
	sys.max_kn = 1000000;
	/*
	 * repulsive force parameter
	 * Short range repulsion is assumed.
	 * cf_amp_dl0: cf_amp_dl at shearrate = 1
	 */
	if (sys.repulsiveforce) {
		sys.set_repulsiveforce_length(0.05);
	} else {
		sys.set_repulsiveforce_length(0);
	}
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	double _mu_static = 1;
	/*
	 * Output interval:
	 * strain_interval_output_data is for outputing rheo_...
	 * strain_interval_output is for outputing int_... and par_...
	 */
	sys.strain_interval_output_data = 0.02;
	sys.strain_interval_output = 0.1;
	/*
	 *  Data output
	 */
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
	origin_zero_flow = true;
	/*
	 * position and interaction data
	 */
	out_data_particle = true;
	out_data_interaction = true;
	sys.set_sd_coeff(_sd_coeff);
	sys.set_integration_method(_integration_method);
	sys.set_lubrication_model(_lubrication_model);
	sys.set_lub_max(_lub_max);
	sys.set_mu_static(_mu_static);
}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	int n1, n2;
	double volume_fraction_;
	double lx, ly, lz;
	double vf1, vf2;
	char buf;
	getline(file_import, import_line[0]);
	getline(file_import, import_line[1]);
	stringstream ss(import_line[1]);
	ss >> buf >> n1 >> n2 >> volume_fraction_ >> lx >> ly >> lz >> vf1 >> vf2;
	volume_fraction = volume_fraction_;
	double x_, y_, z_, a_;
	vector<vec3d> initial_position;
	while (file_import >> x_ >> y_ >> z_ >> a_) {
		initial_position.push_back(vec3d(x_, y_, z_));
		radius.push_back(a_);
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void
Simulation::prepareSimulationName(){
	ostringstream ss_simu_name;
	string::size_type pos_ext_position = filename_import_positions.find(".dat");
	string::size_type pos_ext_parameter = filename_parameters.find(".txt");
	string::size_type pos_ext_sequence = filename_sequence.find(".dat");
	cout << filename_sequence << endl;
	cout << filename_sequence.substr(0, pos_ext_sequence) << endl << endl;

	ss_simu_name << filename_import_positions.substr(0, pos_ext_position);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(0, pos_ext_parameter);
	
	if (control_var == "strain") {
		if (sys.dimensionless_shear_rate == -1) {
			ss_simu_name << "_srinf" ; // shear rate infinity
		} else {
			ss_simu_name << "_sr" << sys.dimensionless_shear_rate;
		}
	}
	if (control_var == "stress") {
		if (user_sequence) {
			ss_simu_name << "_st_" << filename_sequence.substr(0, pos_ext_sequence);
		} else {
			ss_simu_name << "_st" << sys.target_stress;
		}
	}
	sys.simu_name = ss_simu_name.str();
}

void
Simulation::evaluateData(){
	sys.analyzeState();
	sys.calcStress();
	sys.calcLubricationForce();
	/* NOTE:
	 *
	 * The total stress DID not include the contact GU terms,
	 * because we consider that the relative motion is not expected hard spheres
	 * and artificial in the soft-sphere contact model.
	 * [Aug 15, 2013]
	 * In the contact model, force is divided into two parts (spring and dash-pot).
	 * In physics, the total force is important.
	 * Therefore, both should be included for the stress calculation.
	 *
	 */
	total_contact_stressXF = sys.total_contact_stressXF_normal+sys.total_contact_stressXF_tan;
	total_stress = sys.total_hydro_stress;
	total_stress += total_contact_stressXF;
	total_stress += sys.total_contact_stressGU; // added (Aug 15 2013)
	if (sys.repulsiveforce) {
		total_repulsive_stress = sys.total_repulsive_stressXF+sys.total_repulsive_stressGU;
		total_stress += total_repulsive_stress;
	}
	if (sys.brownian) {
		total_stress += sys.total_brownian_stressGU;
	}
	/*
	 * Viscosity is only the increment of stress (=del_eta).
	 * The total viscosity should be
	 * eta_r = eta/eta_0 = 1 + del_eta.
	 */
	viscosity = total_stress.getStressXZ()+5*volume_fraction/(12*M_PI);
	normalstress_diff_1 = total_stress.getNormalStress1();
	normalstress_diff_2 = total_stress.getNormalStress2();
	particle_pressure = total_stress.getParticlePressure();
	viscosity_hydro = sys.total_hydro_stress.getStressXZ();
	normalstress_diff_1_hydro = sys.total_hydro_stress.getNormalStress1();
	normalstress_diff_2_hydro = sys.total_hydro_stress.getNormalStress2();
	viscosity_cont_XF = total_contact_stressXF.getStressXZ();
	normalstress_diff_1_cont_XF = total_contact_stressXF.getNormalStress1();
	normalstress_diff_2_cont_XF = total_contact_stressXF.getNormalStress2();
	particle_pressure_cont = total_contact_stressXF.getParticlePressure();
	viscosity_friction = sys.total_contact_stressXF_tan.getStressXZ();
	normalstress_diff_1_friction = sys.total_contact_stressXF_tan.getNormalStress1();
	normalstress_diff_2_friction = sys.total_contact_stressXF_tan.getNormalStress2();
	viscosity_cont_GU = sys.total_contact_stressGU.getStressXZ();
	normalstress_diff_1_cont_GU = sys.total_contact_stressGU.getNormalStress1();
	normalstress_diff_2_cont_GU = sys.total_contact_stressGU.getNormalStress2();
	if (sys.repulsiveforce) {
		viscosity_repulsive_XF = sys.total_repulsive_stressXF.getStressXZ();
		normalstress_diff_1_repulsive_XF = sys.total_repulsive_stressXF.getNormalStress1();
		normalstress_diff_2_repulsive_XF = sys.total_repulsive_stressXF.getNormalStress2();
		particle_pressure_repulsive = sys.total_repulsive_stressXF.getParticlePressure();
		viscosity_repulsive_GU = sys.total_repulsive_stressGU.getStressXZ();
		normalstress_diff_1_repulsive_GU = sys.total_repulsive_stressGU.getNormalStress1();
		normalstress_diff_2_repulsive_GU = sys.total_repulsive_stressGU.getNormalStress2();
	}
	if (sys.brownian) {
		viscosity_brownian = sys.total_brownian_stressGU.getStressXZ();
		normalstress_diff_1_brownian = sys.total_brownian_stressGU.getNormalStress1();
		normalstress_diff_2_brownian = sys.total_brownian_stressGU.getNormalStress2();
	}
}

void
Simulation::outputStressTensorData(){
	fout_st << sys.get_shear_strain() << ' ';
	fout_st << 6*M_PI*viscosity << ' ';
	/* total_stress = sys.total_hydro_stress;
	 * + total_contact_stressXF + total_repulsive_stress;
	 */

	// As it is, the output stress lacks a 6pi factor (as the viscosity)
	total_stress.outputStressTensor(fout_st); // (3,4,5,6,7,8)
	sys.total_hydro_stress.outputStressTensor(fout_st); // (9,10,11,12,13,14)
	total_contact_stressXF.outputStressTensor(fout_st); // (15,16,17,18,19,20)
	sys.total_contact_stressGU.outputStressTensor(fout_st); // (21,22,23,24,25,26)
	total_repulsive_stress.outputStressTensor(fout_st); // (27,28,29,30,31,32)
	sys.total_brownian_stressGU.outputStressTensor(fout_st); // (33,34,35,36,37,38)
	fout_st << sys.dimensionless_shear_rate/6/M_PI << ' '; // 39
	fout_st << endl;
}

void
Simulation::outputRheologyData(){
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
	 */
	fout_rheo << sys.get_shear_strain() << ' '; //1
	fout_rheo << 6*M_PI*viscosity << ' '; //2
	fout_rheo << 6*M_PI*normalstress_diff_1 << ' '; //3
	fout_rheo << 6*M_PI*normalstress_diff_2 << ' '; //4
	/*
	 * Hydrodynamic contribution means
	 * stresslet_hydro_GU_i+stresslet_ME_i from vel_hydro
	 * vel_hydro is obtained with GE for the rhs.
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
	fout_rheo << sys.min_gap_nondim << ' '; //28
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
	fout_rheo << sys.get_kn() << ' '; //44
	fout_rheo << sys.get_kt() << ' '; //45
	fout_rheo << sys.get_dt() << ' '; //46
	fout_rheo << sys.get_time() << ' ' ; //47
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow) {
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

void
Simulation::outputDataHeader(ofstream &fout){
	fout << "# LF_DEM version " << VERSION << endl;
	fout << "# np " << sys.get_np() << endl;
	fout << "# VF " << volume_fraction << endl;
	fout << "# Lx " << sys.get_lx() << endl;
	fout << "# Ly " << sys.get_ly() << endl;
	fout << "# Lz " << sys.get_lz() << endl;
}

void
Simulation::outputConfigurationData(){
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
	if (origin_zero_flow) {
		for (int i=0; i<np; i++) {
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0) {
				vel[i].x -= sys.get_lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (out_data_particle) {
		fout_particle << "# " << sys.get_shear_strain() << ' ';
		fout_particle << sys.shear_disp << ' ';
		fout_particle << sys.dimensionless_shear_rate << endl;
		for (int i=0; i<np; i++) {
			vec3d &p = pos[i];
			vec3d &v = vel[i];
			vec3d &o = sys.ang_velocity[i];
			double lub_xzstress = sys.lubstress[i].getStressXZ();
			double contact_xzstressGU = sys.contactstressGU[i].getStressXZ();
			double brownian_xzstressGU = 0;
			if (sys.brownian) {
				brownian_xzstressGU = sys.brownianstressGU[i].getStressXZ();
			}
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << ' ' << p.x << ' ' << p.y << ' ' << p.z; //3, 4, 5: position
			fout_particle << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
			fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions
			fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions
			fout_particle << ' ' << 6*M_PI*brownian_xzstressGU; //14: xz stress contributions
			if (sys.twodimension) {
				fout_particle << ' ' << sys.angle[i]; // 15
			}
			fout_particle << endl;
		}
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			cnt_interaction++;
		}
	}
	if (out_data_interaction) {
		fout_interaction << "# " << sys.get_shear_strain();
		fout_interaction << ' ' << cnt_interaction << endl;
		for (int k=0; k<sys.nb_interaction; k++) {
			if (sys.interaction[k].is_active()) {
				unsigned short i, j;
				sys.interaction[k].get_par_num(i, j);
				vec3d nr_vec = sys.interaction[k].get_nvec();
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
				fout_interaction << sys.interaction[k].get_gap_nondim() << ' '; // 7
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
				fout_interaction << sys.interaction[k].get_f_repulsive_norm() << ' '; // 16
				fout_interaction << 6*M_PI*stress_contact.getStressXZ() << ' '; // 17
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << ' ';
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << ' ';
				fout_interaction << endl;
			}
		}
	}
}

void
Simulation::outputFinalConfiguration(){
	ofstream fout_finalconfig;
	string filename_final_configuration = "./after_relax/"+filename_import_positions;
	fout_finalconfig.open(filename_final_configuration.c_str());
	fout_finalconfig << import_line[0] << endl;
	fout_finalconfig << import_line[1] << endl;
	int np = sys.get_np();
	for (int i = 0; i < np; i++) {
		fout_finalconfig << sys.position[i].x << ' ';
		fout_finalconfig << sys.position[i].y << ' ';
		fout_finalconfig << sys.position[i].z << ' ';
		fout_finalconfig << sys.radius[i] << endl;
	}
}

