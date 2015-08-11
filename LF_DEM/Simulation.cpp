//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include "Simulation.h"
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <cctype>

using namespace std;

Simulation::Simulation():
sys(System(p,events)),
shear_rate_expectation(-1),
internal_unit_scales("hydro"),
target_stress_input(0)
{
	unit_longname["h"] = "hydro";
	unit_longname["r"] = "repulsive";
	unit_longname["b"] = "thermal";
	unit_longname["c"] = "cohesive";
	unit_longname["cl"] = "critical_load";
	unit_longname["m"] = "magnetic";
	unit_longname["ft"] = "ft";

	kill = false;
};

Simulation::~Simulation()
{
	if (fout_data.is_open()) {
		fout_data.close();
	}
	if (fout_particle.is_open()) {
		fout_particle.close();
	}
	if (fout_interaction.is_open()) {
		fout_interaction.close();
	}
};


bool Simulation::keepRunning()
{
	if (time_end == -1) {
		return (sys.get_shear_strain() < strain_end-1e-8) && !kill;
	}
	else {
		return (sys.get_time() < time_end-1e-8) && !kill;
	}
}

void Simulation::setupEvents(){
	if (p.event_handler == "shear_jamming") {
		sys.eventLookUp = &System::eventShearJamming;
		return;
	}
	sys.eventLookUp = NULL;
}

void Simulation::handleEvents(){
	if (p.event_handler == "shear_jamming") {
		for (const auto &ev : events) {
			if (ev.type == "negative_shear_rate") {
				cout << " negative rate " << endl;
				p.disp_max /= 1.1;
			}
		}
		if(p.disp_max < 1e-6){
			cout << "jammed" << endl;
//			kill = true;
			p.cross_shear = true;
		}
	}
	events.clear();
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
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale);
	if (sys.cohesion) {
		sys.new_contact_gap = 0.02; //@@ To be changed to a better way.
	} else {
		sys.new_contact_gap = 0;
	}
//	int jammed = 0;
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;
	/******************** OUTPUT INITIAL DATA ********************/
	evaluateData(); //
	outputData(); // new
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	setupEvents();

	double next_output_data = 0;
	double next_output_config = 0;

	while (keepRunning()) {
		if (time_interval_output_data == -1) {
			next_output_data += strain_interval_output_data;
			sys.timeEvolution("strain", next_output_data);
		} else {
			next_output_data +=  time_interval_output_data;
			sys.timeEvolution("time", next_output_data);
		}
		handleEvents();

		/******************** OUTPUT DATA ********************/
		evaluateData();
		outputData(); // new
		outputConfigurationBinary();
		if (time_interval_output_config == -1) {
			if (sys.get_shear_strain() >= next_output_config-1e-8) {
				outputConfigurationData();
				next_output_config += strain_interval_output_config;
			}
		} else {
			if (sys.get_time() >= next_output_config-1e-8) {
				outputConfigurationData();
				next_output_config +=  time_interval_output_config;
			}
		}
		/*****************************************************/
		if (time_end != -1) {
			cout << "time: " << sys.get_time() << " / " << time_end << " , strain: " << sys.get_shear_strain() << endl;
		} else {
			cout << "time: " << sys.get_time() << " , strain: " << sys.get_shear_strain()  << " / " << strain_end << endl;
		}

		// if (!sys.zero_shear
		// 	&& abs(sys.get_shear_rate()) < p.rest_threshold){
		// 	cout << "shear jamming " << jammed << endl;
		// 	jammed ++;
		// 	if (jammed > 10) {
		// 		cout << "shear jamming";
		// 		break;
		// 	}
		// } else {
		// 	jammed = 0;
		// }
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

void Simulation::simulationInverseYield(string in_args,
										vector<string> &input_files,
										bool binary_conf,
										double dimensionless_number,
										string input_scale,
										string control_variable)
{
	user_sequence = false;
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale);

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
	outputData(); // new
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	double next_output_data = 0;
	double next_output_config = 0;

	while (keepRunning()) {
		if (time_interval_output_data == -1) {
			next_output_data += strain_interval_output_data;
			sys.timeEvolution("strain", next_output_data);
		} else{
			next_output_data +=  time_interval_output_data;
			sys.timeEvolution("time", next_output_data);
		}

		/******************** OUTPUT DATA ********************/
		evaluateData();
		outputData(); // new
		outputConfigurationBinary();
		if (time_interval_output_config == -1) {
			if (sys.get_shear_strain() >= next_output_config-1e-8) {
				outputConfigurationData();
				next_output_config += strain_interval_output_config;
			}
		} else {
			if (sys.get_time() >= next_output_config-1e-8) {
				outputConfigurationData();
				next_output_config +=  time_interval_output_config;
			}
		}
		/*****************************************************/

		cout << "time: " << sys.get_time() << " / " << p.time_end << endl;
		if (!sys.zero_shear
			&& abs(sys.get_shear_rate()) < p.rest_threshold) {
			cout << "shear jamming " << jammed << endl;
			jammed ++;
			if (jammed > 20) {
				sys.set_shear_rate(1);
				cout << "target_stress = " << target_stress_input << endl;
				target_stress_input *= 0.95;
				sys.target_stress = target_stress_input/6/M_PI;
				sys.updateUnscaledContactmodel();
				cout << "new target_stress = " << target_stress_input << endl;
				jammed = 0;
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

void Simulation::simulationMagnetic(string in_args,
									vector<string> &input_files,
									bool binary_conf,
									double dimensionless_number,
									string input_scale,
									string control_variable)
{
	/* Monolayer: Particles are confined in y = 0 plane.
	 *
	 *
	 */
	user_sequence = false;
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale);
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	double strain_output_config = 0;
	double time_output_data = 0;
	double time_output_config = 0;
	/******************** OUTPUT INITIAL DATA ********************/
	evaluateData();
	outputData(); // new
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	double external_magnetic_field_norm = 2*sqrt(2*dimensionless_number); //@@@@@@ TO BE CHECKED!!
	if (p.magnetic_field_type == 0) {
		// Field direction is fixed
		sys.external_magnetic_field.set(0, external_magnetic_field_norm, 0);
	} else if (p.magnetic_field_type == 1) {
		sys.external_magnetic_field.set(external_magnetic_field_norm*sin(sys.angle_external_magnetic_field),
										external_magnetic_field_norm*cos(sys.angle_external_magnetic_field),
										0);
		exit(1);
	} else if (p.magnetic_field_type == 2) {
		sys.external_magnetic_field.set(external_magnetic_field_norm*cos(sys.angle_external_magnetic_field),
										0,
										external_magnetic_field_norm*sin(sys.angle_external_magnetic_field));
		exit(1);
	}
	sys.setMagneticMomentExternalField();
	// Main simulation loop
	while (keepRunning()) {
		time_output_data = cnt_simu_loop*time_interval_output_data;
		time_output_config = cnt_config_out*time_interval_output_config;
		sys.timeEvolution(time_output_data);
		cnt_simu_loop ++;
		/******************** OUTPUT DATA ********************/
		evaluateData();
		outputData();
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
		cout << "time: " << sys.get_time() << " / " << time_end << endl;
	}

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



void Simulation::catchSuffixedValue(string type, string keyword, string value_str, double *value_ptr){
	InputValue inv;
	inv.type = type;
	inv.name = keyword;
	inv.value = value_ptr;

	string numeral, suffix;
	bool caught_suffix = true;
	caught_suffix = getSuffix(value_str, numeral, suffix);
	suffix = unit_longname[suffix];
	*(inv.value) = atof(numeral.c_str());
	inv.unit = suffix;
	input_values.push_back(inv);

	if (!caught_suffix) {
		errorNoSuffix(keyword);
	}
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
	conf_export.open(conf_filename.c_str(), ios::binary | ios::out);
	conf_export.write((char*)&np, sizeof(int));
	conf_export.write((char*)&volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&lx, sizeof(double));
	conf_export.write((char*)&ly, sizeof(double));
	conf_export.write((char*)&lz, sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.x), sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.y), sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	conf_export.close();
}


double Simulation::getRate(){
/**
	\brief The shear rate in the input units
*/
	if (control_var == "rate") {
		return input_rate;
	} else if (control_var == "stress") {
		return sys.get_shear_rate();
	} else {
		return 1;
	}
}

void Simulation::evaluateData()
{
	/**
	 \brief Get rheological data from the System class.

	 In this method we keep the internal units. There is no conversion to output units at this stage

	 */

	sys.analyzeState();
	sys.calcStress();
	sys.calcLubricationForce();
}


void Simulation::outputData()
{
	/**
	 \brief Output data.

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

	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;
//	cerr << internal_unit_scales << " " << output_unit_scales << endl;

	if (dimensionless_numbers.find(dimless_nb_label) == dimensionless_numbers.end()) {
		cerr << " Error : don't manage to convert from \"" << internal_unit_scales << "\" units to \"" << output_unit_scales << "\" units to output data." << endl; exit(1);
	}
	outdata.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);

	if (p.magnetic_type == 0) {
		outdata.init(36, output_unit_scales);
	} else {
		outdata.init(40, output_unit_scales);
	}

	double sr = sys.get_shear_rate();
	unsigned int shear_stress_index;
	if (!p.cross_shear) {
		shear_stress_index = 2;
	} else {
		shear_stress_index = 3;
	}
	double shear_stress = 6*M_PI*(sys.einstein_stress+sys.total_stress.elm[shear_stress_index]);

	outdata.entryData(1, "time", "time", sys.get_time());
	outdata.entryData(2, "shear strain", "none", sys.get_shear_strain());
	outdata.entryData(3, "shear rate", "rate", sys.get_shear_rate());

	outdata.entryData(5, "viscosity", "viscosity", shear_stress/sr);
	outdata.entryData(6, "Viscosity(lub)", "viscosity", sys.total_hydro_stress.elm[shear_stress_index]/sr);
	outdata.entryData(7, "Viscosity(xF_contact part)", "viscosity", sys.total_contact_stressXF.elm[shear_stress_index]/sr);
	outdata.entryData(8, "Viscosity(GU_contact part)", "viscosity", sys.total_contact_stressGU.elm[shear_stress_index]/sr);
	if (sys.repulsiveforce) {
		outdata.entryData(9, "Viscosity(repulsive force XF)", "viscosity", sys.total_repulsive_stressXF.elm[shear_stress_index]/sr);
		outdata.entryData(10, "Viscosity(repulsive force GU)", "viscosity", sys.total_repulsive_stressGU.elm[shear_stress_index]/sr);
	}
	if (sys.brownian) {
		outdata.entryData(11, "Viscosity(brownian)", "viscosity", sys.total_brownian_stressGU.elm[shear_stress_index]/sr);
	}
	/*
	 * Stress
	 */
	outdata.entryData(14, "shear stress", "stress", shear_stress);
	outdata.entryData(15, "N1 viscosity", "viscosity", sys.total_stress.getNormalStress1()/sr);
	outdata.entryData(16, "N2 viscosity", "viscosity", sys.total_stress.getNormalStress2()/sr);
	outdata.entryData(17, "particle pressure", "stress", sys.total_stress.getParticlePressure());
	outdata.entryData(18, "particle pressure contact", "stress", sys.total_contact_stressXF.getParticlePressure());
	/* energy
	 */
	outdata.entryData(21, "energy", "none", sys.get_total_energy());
	/* maximum deformation of contact bond
	 */
	outdata.entryData(22, "min gap", "none", sys.min_reduced_gap);
	outdata.entryData(23, "max gap(cohesion)", "none", sys.max_contact_gap);
	outdata.entryData(24, "max tangential displacement", "none", sys.max_disp_tan);
	outdata.entryData(25, "max rolling displacement", "none", sys.max_disp_rolling);
	/* contact number
	 */
	outdata.entryData(26, "contact number", "none", sys.getContactNumber());
	outdata.entryData(27, "frictional contact number", "none", sys.getFrictionalContactNumber());
	outdata.entryData(28, "number of interaction", "none", sys.get_nb_of_active_interactions());
	/* maximum velocity
	 */
	outdata.entryData(29, "max velocity", "velocity", sys.max_velocity);
	outdata.entryData(30, "max angular velocity", "velocity", sys.max_ang_velocity);
	/* simulation parameter
	 */
	outdata.entryData(31, "dt", "time", sys.dt);
	outdata.entryData(32, "kn", "none", p.kn);
	outdata.entryData(33, "kt", "none", p.kt);
	outdata.entryData(34, "kr", "none", p.kr);
	outdata.entryData(35, "shear displacement", "none", sys.shear_disp.x);
	if (p.magnetic_type != 0) {
		outdata.entryData(37, "magnetic energy", "none", sys.magnetic_energy);
		outdata.entryData(38, "magnetic field angle", "none", sys.angle_external_magnetic_field);
	}
	outdata.exportFile(fout_data);

	/****************************   Stress Tensor Output *****************/
	outdata_st.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);

   	outdata_st.init(8, output_unit_scales);

	outdata_st.entryData(1, "time", "time", sys.get_time());
	outdata_st.entryData(2, "shear strain", "none", sys.get_shear_strain());
	outdata_st.entryData(3, "shear rate", "rate", sys.get_shear_rate());
	outdata_st.entryData(4, "total stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_stress);
	outdata_st.entryData(5, "hydro stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_hydro_stress);
	outdata_st.entryData(6, "contact stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_contact_stressXF+sys.total_contact_stressGU);
	outdata_st.entryData(7, "repulsive stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_repulsive_stress);
	outdata_st.entryData(8, "brownian stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_brownian_stressGU);
	outdata_st.exportFile(fout_st);
}


vec3d Simulation::shiftUpCoordinate(double x, double y, double z)
{
	if (p.origin_zero_flow) {
		z += sys.Lz_half();
		if (z > sys.Lz_half()) {
			x -= sys.shear_disp.x;
			y -= sys.shear_disp.y;
			if (x < -sys.Lx_half()) {
				x += sys.get_lx();
			}
			if (y < -sys.Ly_half()) {
				y += sys.get_ly();
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
	for (int i=0; i<np; i++) {
		vel[i] = sys.velocity[i];
		if (p.origin_zero_flow) {
			if (pos[i].z < 0) {
				vel[i].x -= sys.get_shear_rate()*sys.get_lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (p.out_data_particle) {
		cout << "   out config: " << sys.get_shear_strain() << endl;
		fout_particle << "# " << sys.get_shear_strain() << ' ';
		fout_particle << sys.shear_disp.x << ' ';
		fout_particle << getRate() << ' ';
		fout_particle << target_stress_input << ' ';
		fout_particle << sys.get_time() << ' ';
		if (p.magnetic_type != 0) {
			fout_particle << sys.angle_external_magnetic_field;
		}
		fout_particle << endl;

		unsigned int shear_stress_index;
		if (!p.cross_shear) {
			shear_stress_index = 2;
		} else {
			shear_stress_index = 3;
		}
		for (int i=0; i<np; i++) {
			const vec3d &r = pos[i];
			const vec3d &v = vel[i];
			const vec3d &o = sys.ang_velocity[i];
			double lub_xzstress = sys.lubstress[i].elm[shear_stress_index];
			double contact_xzstressGU = sys.contactstressGU[i].elm[shear_stress_index];
			double brownian_xzstressGU = 0;
			if (sys.brownian) {
				brownian_xzstressGU = sys.brownianstressGU[i].elm[shear_stress_index];
			}
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << ' ' << r.x << ' ' << r.y << ' ' << r.z; //3, 4, 5: position
			fout_particle << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
			fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions
			fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions
			fout_particle << ' ' << 6*M_PI*brownian_xzstressGU; //14: xz stress contributions
			if (p.magnetic_type == 0) {
				if (sys.twodimension) {
					fout_particle << ' ' << sys.angle[i]; // 15
				}
			} else {
				fout_particle << ' ' << sys.magnetic_moment[i].x;
				fout_particle << ' ' << sys.magnetic_moment[i].y;
				fout_particle << ' ' << sys.magnetic_moment[i].z;
				fout_particle << ' ' << sys.magnetic_susceptibility[i];
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
		unsigned int shear_stress_index;
		if (!p.cross_shear) {
			shear_stress_index = 2;
		} else {
			shear_stress_index = 3;
		}
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
				fout_interaction << 6*M_PI*stress_contact.elm[shear_stress_index] << ' '; // 17
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
	fout_finalconfig << header_imported_configulation[0] << endl;
	fout_finalconfig << header_imported_configulation[1] << endl;
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
