//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>
#include <cmath>
#include <stdexcept>
#ifdef USE_DSFMT
#include <time.h>
#endif
#define DELETE(x) if(x){delete [] x; x = NULL;}
#ifndef USE_DSFMT
#define GRANDOM ( r_gen->randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.
#define RANDOM ( rand_gen.rand() ) // RNG uniform [0,1]
#endif
#ifdef USE_DSFMT
#define GRANDOM  ( sqrt( -2.0 * log( 1.0 - dsfmt_genrand_open_open(&r_gen) ) ) * cos(2.0 * 3.14159265358979323846264338328 * dsfmt_genrand_close_open(&r_gen) ) ) // RNG gaussian with mean 0. and variance 1.
#define RANDOM ( dsfmt_genrand_close_open(&rand_gen) ) // RNG uniform [0,1]
#endif

using namespace std;

System::System(ParameterSet& ps, list <Event>& ev):
events(ev),
p(ps),
brownian(false),
friction(false),
rolling_friction(false),
repulsiveforce(false),
magnetic(false),
cohesion(false),
critical_load(false),
lowPeclet(false),
twodimension(false),
zero_shear(false),
target_stress(0),
init_strain_shear_rate_limit(0),
init_shear_rate_limit(999),
new_contact_gap(0),
magnetic_rotation_active(false),
eventLookUp(NULL)
{
	amplitudes.repulsion = 0;
	amplitudes.sqrt_temperature = 0;
	amplitudes.contact = 0;
	amplitudes.cohesion = 0;
	amplitudes.magnetic = 0;
	amplitudes.critical_normal_force = 0;
	max_sliding_velocity = 0;
	max_contact_gap = 0;
	max_disp_rolling = 0;
}

vec3d System::randUniformSphere(double r)
{
	double z = 2*RANDOM-1;
	double phi = 2*M_PI*RANDOM;
	double sin_theta = sqrt(1-z*z);
	return vec3d(r*sin_theta*cos(phi), r*sin_theta*sin(phi), r*z);
}

#ifdef USE_DSFMT
unsigned long
System::hash( time_t t, clock_t c )
{
	/**
		\brief Utility function to start up the DSFMT RNG with a nice seed.

	 From MersenneTwister v1.0 by Richard J. Wagner
	 comments below are from the original code.

	 Get a unsigned long from t and c
	 Better than unsigned long(x) in case x is floating point in [0,1]
	 Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
	*/

	static unsigned long differ = 0; // guarantee time-based seeds will change
	unsigned long h1 = 0;
	unsigned char *pp = (unsigned char *) &t;
	for( size_t i = 0; i < sizeof(t); ++i )
	{
		h1 *= UCHAR_MAX + 2U;
		h1 += pp[i];
	}
	unsigned long h2 = 0;
	pp = (unsigned char *) &c;
	for( size_t j = 0; j < sizeof(c); ++j )
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += pp[j];
	}
	return ( h1 + differ++ ) ^ h2;
}
#endif

System::~System()
{
	DELETE(position);
	DELETE(radius);
	DELETE(radius_squared);
	DELETE(radius_cubed);
	if (twodimension) {
		DELETE(angle);
	}
	DELETE(velocity);
	DELETE(ang_velocity);
	DELETE(na_velocity);
	DELETE(na_ang_velocity);
	if (p.integration_method == 1) {
		DELETE(velocity_predictor);
		DELETE(ang_velocity_predictor);
	}
	DELETE(vel_contact);
	DELETE(ang_vel_contact);
	DELETE(vel_hydro);
	DELETE(ang_vel_hydro);
	DELETE(contact_force);
	DELETE(contact_torque);
	DELETE(lubstress);
	DELETE(contactstressGU);
	if (!p.out_particle_stress.empty()) {
		DELETE(contactstressXF);
	}
	DELETE(interaction);
	DELETE(interaction_list);
	DELETE(interaction_partners);
	if (brownian) {
		DELETE(vel_brownian);
		DELETE(ang_vel_brownian);
		DELETE(brownian_force);
		DELETE(brownianstressGU);
		DELETE(brownianstressGU_predictor);
	}
	if (repulsiveforce) {
		DELETE(repulsive_force);
		DELETE(repulsivestressGU);
		if (!p.out_particle_stress.empty()) {
			DELETE(repulsivestressXF);
		}
		DELETE(vel_repulsive);
		DELETE(ang_vel_repulsive);
	}
	if (magnetic) {
		DELETE(magnetic_moment);
		DELETE(magnetic_force);
		DELETE(magnetic_torque);
		DELETE(vel_magnetic);
		DELETE(ang_vel_magnetic);
		DELETE(magneticstressGU);
	}
};

void System::allocateRessources()
{
	linalg_size = 6*np;
	double interaction_volume;
	if (twodimension) {
		interaction_volume = M_PI*pow(p.interaction_range, 2);
		double particle_volume = M_PI;
		maxnb_interactionpair_per_particle = interaction_volume/particle_volume;
	} else {
		interaction_volume = (4*M_PI/3)*pow(p.interaction_range, 3);
		double particle_volume = 4*M_PI/3;
		maxnb_interactionpair_per_particle = 1*interaction_volume/particle_volume;
	}
	maxnb_interactionpair = maxnb_interactionpair_per_particle*np;
	radius_cubed = new double [np];
	radius_squared = new double [np];
	if (p.lubrication_model == 0) {
		// Stokes-drag simulation
		stokesdrag_coeff_f = new double [np];
		stokesdrag_coeff_f_sqrt = new double [np];
		stokesdrag_coeff_t = new double [np];
		stokesdrag_coeff_t_sqrt = new double [np];
	}
	// Configuration
	if (twodimension) {
		angle = new double [np];
	}
	// Velocity
	velocity = new vec3d [np];
	ang_velocity = new vec3d [np];
	na_velocity = new vec3d [np];
	na_ang_velocity = new vec3d [np];
	if (p.integration_method == 1) {
		velocity_predictor = new vec3d [np];
		ang_velocity_predictor = new vec3d [np];
	}
	vel_contact = new vec3d [np];
	ang_vel_contact = new vec3d [np];
	vel_hydro = new vec3d [np];
	ang_vel_hydro = new vec3d [np];
	if (repulsiveforce) {
		vel_repulsive = new vec3d [np];
		ang_vel_repulsive = new vec3d [np];
	}
	if (brownian) {
		vel_brownian = new vec3d [np];
		ang_vel_brownian = new vec3d [np];
	}
	if (magnetic) {
		vel_magnetic = new vec3d [np];
		ang_vel_magnetic = new vec3d [np];
	}
	// Forces and Stress
	contact_force = new vec3d [np];
	contact_torque = new vec3d [np];
	lubstress = new StressTensor [np];
	contactstressGU = new StressTensor [np];
	if (!p.out_particle_stress.empty()) {
		contactstressXF = new StressTensor [np];
	}
	if (repulsiveforce) {
		repulsive_force = new vec3d [np];
		repulsivestressGU = new StressTensor [np];
		if (!p.out_particle_stress.empty()) {
			repulsivestressXF = new StressTensor [np];
		}
	}
	if (brownian) {
		brownian_force = new double [linalg_size];
		brownianstressGU = new StressTensor [np];
		brownianstressGU_predictor = new StressTensor [np];
	}
	if (magnetic) {
		magnetic_force = new vec3d [np];
		magnetic_torque = new vec3d [np];
		magneticstressGU = new StressTensor [np];
	}
	//
	interaction = new Interaction [maxnb_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];

	//
	if (p.auto_determine_knkt) {
		kn_avg = new Averager<double>(p.memory_strain_avg);
		kt_avg = new Averager<double>(p.memory_strain_avg);
		overlap_avg = new Averager<double>(p.memory_strain_avg);
		max_disp_tan_avg = new Averager<double>(p.memory_strain_avg);
	}
	if (brownian) {
		if (lowPeclet) {
			double stress_avg_relaxation_parameter = 10*p.time_interval_output_data; // 0 --> no average
			stress_avg = new Averager<StressTensor>(stress_avg_relaxation_parameter);
		}
	}
	movable.resize(np);
	for (int i=0; i<np; i++) {
		movable[i] = true;
	}
}

void System::setInteractions_GenerateInitConfig()
{
	calcInteractionRange = &System::calcLubricationRange;
	for (int k=0; k<maxnb_interactionpair; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	nb_interaction = 0;
	shear_strain = 0;
	shear_disp.reset();
	vel_difference.reset();
	initializeBoxing();
	checkNewInteraction();
}

void System::allocatePositionRadius()
{
	position = new vec3d [np];
	radius = new double [np];
}

void System::setConfiguration(const vector <vec3d>& initial_positions,
							  const vector <double>& radii,
							  double lx_, double ly_, double lz_)
{
	/**
		\brief Set positions of the particles for initialization.
	 */
	string indent = "  System::\t";
	set_np(initial_positions.size());
	setBoxSize(lx_, ly_, lz_);
	allocatePositionRadius();
	for (int i=0; i<np; i++) {
		position[i] = initial_positions[i];
		radius[i] = radii[i];
	}
	if (ly_ == 0) {
		twodimension = true;
	} else {
		twodimension = false;
	}
	if (twodimension) {
		// @@@ this is arbitrary. Is it really necessary?
		// @@@ In particular, the stress calculation uses the system volume, and most users expect the volume to be lx*lz in this case.
		// @@@-----------------------------------------------------------------------
		// @@@ I see. I though monolayer system including spheres is easier to switch between 2D and 3D simulations.
		// @@@ But, I agree keeping such useless dimension looks bit awkward.
		// @@@ Since I considered 2D version is useful just for visual demos, I did not care about this part.
		// @@@ Please check consistency of 2D rheology.
		// @@@ I don't expect clyinders instead of spheres.
		// @@@ Probably We can remove the 2.5 phi term vfrom both 2D and 3D versions.
		//
		/* [note]
		 * The depth of mono-layer is the diameter of the largest particles.e
		 * The sample needs to be labeled from smaller particles to larger particles
		 * in the configuration file.
		 */
		double largest_diameter = 2*radius[np-1];
		setSystemVolume(largest_diameter);
	} else {
		setSystemVolume();
	}
	double particle_volume = 0;
	for (int i=0; i<np; i++) {
		particle_volume += (4*M_PI/3)*pow(radius[i],3);
	}
	volume_fraction = particle_volume/system_volume;
	cout << indent << "volume_fraction = " << volume_fraction << endl;
}

void System::setContacts(const vector <struct contact_state>& cs)
{
	/**
		\brief Set a list of contacts with their state variables.

		Used to restart the simulation from a given state.
	 */
	for (const auto& c : cs) {
		for (int k=0; k<nb_interaction; k++) {
			unsigned short p0, p1;
			interaction[k].get_par_num(p0, p1);
			if (p0 == c.p0 && p1 == c.p1) {
				interaction[k].contact.setState(c);
			}
		}
	}
}

void System::getContacts(vector <struct contact_state>& cs)
{
	/**
		\brief Get the list of contacts with their state variables.

		Used to output a configuration including contact info. Useful if you want to restart from exact same configuration.
	 */
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			cs.push_back(interaction[k].contact.getState());
		}
	}
}

void System::setInducedMagneticMoment()
{
	for (int i=0; i<np; i++) {
		magnetic_moment[i] = magnetic_susceptibility[i]*external_magnetic_field;
	}
}

void System::setMagneticMomentZero()
{
	for (int i=0; i<np; i++) {
		magnetic_moment[i].reset();
	}
}

void System::setMagneticConfiguration(const vector <vec3d>& magnetic_moment_,
									  const vector <double>& magnetic_susceptibility_)
{
	magnetic_susceptibility.resize(np);
	magnetic_moment = new vec3d [np];
	magnetic_pair.resize(np-1);
	sq_magnetic_interaction_range = pow(p.magnetic_interaction_range, 2);
	time_update_magnetic_pair = 0;
	if (p.magnetic_type == 1) {
		/* Each particle has magnetic dipole moment.
		 * Ferromagnetism.
		 */
		for (int i=0; i<np; i++) {
			magnetic_moment[i] = magnetic_moment_[i];
		}
		throw runtime_error("This is not implemented yet");
		//num_magnetic_particles = i_magnetic;
	} else if (p.magnetic_type == 2) {
		/* Particle can have magnetic moment when external magnetic field is applied.
		 * Paramagnetism
		 * Magnetic susceptibility.
		 * Magnetic moments are fixed as long as external field is unchanged.
		 */
		for (int i=0; i<np; i++) {
			magnetic_susceptibility[i] = magnetic_susceptibility_[i];
		}
	}
}

void System::randomSelectFixParticles()
{
	if (p.fixed_particle_fraction == 0) {
		return;
	}
	int num_fixed_particle = (int)(np*p.fixed_particle_fraction);
	int counter = 0;
	do {
		int picked_particle = (int)(drand48()*np);
		if (movable[picked_particle] == true) {
			movable[picked_particle] = false;
			counter++;
		}
	} while (counter < num_fixed_particle);
}

void System::updateUnscaledContactmodel()
{
	if (abs(target_stress) != 0) {
		/* What is the reasonable way
		 * when the target stress is changed during a simulation?
		 *
		 * [temporally change]
		 * For destressing tests, the spring constants are fixed at the previous values.
		 */
		if (!cohesion) {
			p.kn = kn_master*abs(target_stress);
			p.kt = kt_master*abs(target_stress);
			p.kr = kr_master*abs(target_stress);
		} else {
			p.kn = kn_master;
			p.kt = kt_master;
			p.kr = kr_master;
		}
		cout << " kn " << p.kn << "  kn_master " << kn_master << " target_stress "  << target_stress << endl;
	}
	lub_coeff_contact = 4*p.kn*p.contact_relaxation_time;
	if (lowPeclet) {
		lub_coeff_contact *= p.Pe_switch;
	}
	if (p.lubrication_model > 0) {
		if (p.lubrication_model == 1) {
			log_lub_coeff_contact_tan_lubrication = 0;
			log_lub_coeff_contact_tan_dashpot = 0;
		} else if (p.lubrication_model == 2) {
			log_lub_coeff_contact_tan_lubrication = log(1/p.lub_reduce_parameter);
			/* [Note]
			 * We finally do not introduce a dashpot for the sliding mode.
			 * This is set in the parameter file, i.e. p.contact_relaxation_time_tan = 0
			 * So log_lub_coeff_contact_tan_dashpot = 0;
			 */
			log_lub_coeff_contact_tan_dashpot = 6*p.kt*p.contact_relaxation_time_tan;
		} else {
			throw runtime_error("Error: lubrication_model>2 ???");
		}
	}
	log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].contact.setInteractionData();
		}
	}
}

void System::setupSystem(string control)
{
	/**
		\brief Initialize the system class for the simulation.
	 */

	/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 * @ We have to consider p.contact_relaxation_time in Brownian case.
	 * @ The resistance coeffient affects Brownian force.
	 * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 *
	 * @@@ --> I agree. I want to make clear the best way to handle contacting hard-sphere Brownian aprticles.
	 * @@@     The contact force parameters kn and contact_relaxation_time may need to depend on dt in Brownian simulation.
	 */
	string indent = "  System::\t";
	cout << indent << "Setting up System... " << endl;
	if (control != "magnetic") {
		if (control == "rate") {
			rate_controlled = true;
		}
		if (control == "stress") {
			rate_controlled = false;
		}
		stress_controlled = !rate_controlled;
	} else {
		rate_controlled = false;
		stress_controlled = false;
		magnetic = true;
	}
	if (p.integration_method == 0) {
		timeEvolutionDt = &System::timeEvolutionEulersMethod;
	} else if (p.integration_method == 1) {
		timeEvolutionDt = &System::timeEvolutionPredictorCorrectorMethod;
	} else {
		ostringstream error_str;
		error_str << indent << "integration_method = " << p.integration_method << endl << indent << "The integration method is not impremented yet." << endl;
		throw runtime_error(error_str.str());
	}
	if (p.lubrication_model == 0) {
		/* Stokes drag simulation
		 * Resistance matrix is constant.
		 */
	} else if (p.lubrication_model == 1) {
		buildLubricationTerms = &System::buildLubricationTerms_squeeze;
	} else if (p.lubrication_model == 2) {
		buildLubricationTerms = &System::buildLubricationTerms_squeeze_tangential;
	} else {
		throw runtime_error(indent+"lubrication_model >= 3 is not implemented yet.\n");
	}
	if (p.interaction_range == -1) {
		/* If interaction_range is not indicated,
		 * interaction object is created at the lubrication cutoff.
		 */
		calcInteractionRange = &System::calcLubricationRange;
		p.interaction_range = 2+p.lub_max_gap;
	} else {
		calcInteractionRange = &System::calcInteractionRangeDefault;
	}
	if (p.friction_model == 0) {
		cout << indent+"friction model: no friction" << endl;
		p.mu_static = 0;
		friction = false;
	} else if (p.friction_model == 1) {
		cout << indent+"friction model: Coulomb" << endl;
		friction = true;
	} else if (p.friction_model == 2 || p.friction_model == 3) {
		cout << indent+"friction model: Coulomb + Critical Load" << endl;
		friction = true;
	} else if (p.friction_model == 5) {
		cout << indent+"friction_model: Max tangential force" << endl;
		friction = true;
	} else if (p.friction_model == 6) {
		cout << indent+"friction_model: Coulomb law + Max tangential force" << endl;
		friction = true;
	} else {
		throw runtime_error(indent+"Error: unknown friction model\n");
	}
	if (p.mu_rolling > 0) {
		rolling_friction = true;
		if (friction == false) {
			throw runtime_error(indent+"Error: Rolling friction without sliding friction?\n");
		}
	}
	if (p.lub_max_gap >= 1) {
		throw runtime_error(indent+"lub_max_gap must be smaller than 1\n");
	}
	if (p.repulsive_length <= 0) {
		repulsiveforce = false;
		p.repulsive_length = 0;
	}
	costheta_shear = cos(p.theta_shear);
	sintheta_shear = sin(p.theta_shear);

	// Memory
	allocateRessources();
	//
	for (int k=0; k<maxnb_interactionpair; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	for (int i=0; i<np; i++) {
		radius_squared[i] = pow(radius[i], 2);
		radius_cubed[i] = pow(radius[i], 3);
		if (twodimension) {
			angle[i] = 0;
		}
		velocity[i].reset();
		na_velocity[i].reset();
		ang_velocity[i].reset();
		na_ang_velocity[i].reset();
		vel_contact[i].reset();
		ang_vel_contact[i].reset();
		vel_hydro[i].reset();
		ang_vel_hydro[i].reset();
		if (repulsiveforce) {
			vel_repulsive[i].reset();
			ang_vel_repulsive[i].reset();
		}
		if (magnetic) {
			vel_magnetic[i].reset();
			ang_vel_magnetic[i].reset();
		}
	}
	shear_strain = 0;
	nb_interaction = 0;
	if (p.unscaled_contactmodel) {
		kn_master = p.kn;
		kt_master = p.kt;
		kr_master = p.kr;
		// cout << indent+" kn " << p.kn << "  kn_master " << kn_master << " target_stress "  << target_stress << endl;
	}
	if (p.contact_relaxation_time < 0) {
		// 1/(h+c) --> 1/c
		lub_coeff_contact = 1/p.lub_reduce_parameter;
	} else {
		/* t = beta/kn
		 *  beta = t*kn
		 * lub_coeff_contact = 4*beta = 4*kn*p.contact_relaxation_time
		 */
		lub_coeff_contact = 4*p.kn*p.contact_relaxation_time;
	}
	/* If a contact is in sliding mode,
	 * lubrication and dashpot forces are activated.
	 * `log_lub_coeff_contactlub' is the parameter for lubrication during dynamic friction.
	 *
	 */
	if (p.lubrication_model == 0) {
		/* Stokes drag simulation
		 */
	} else if (p.lubrication_model == 1) {
		log_lub_coeff_contact_tan_lubrication = 0;
		log_lub_coeff_contact_tan_dashpot = 0;
		log_lub_coeff_contact_tan_total = 0;
	} else if (p.lubrication_model == 2) {
		log_lub_coeff_contact_tan_lubrication = log(1/p.lub_reduce_parameter);
		/* [Note]
		 * We finally do not introduce a dashpot for the sliding mode.
		 * This is set in the parameter file, i.e. p.contact_relaxation_time_tan = 0
		 * So log_lub_coeff_contact_tan_dashpot = 0;
		 */
		log_lub_coeff_contact_tan_dashpot = 6*p.kt*p.contact_relaxation_time_tan;
		log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	} else {
		throw runtime_error("lubrication_model must be smaller than 3\n");
	}

	if (p.unscaled_contactmodel) {
		updateUnscaledContactmodel();
	}
	// cout << "lub_coeff_contact = " << lub_coeff_contact << endl;
	// cout << "1/lub_reduce_parameter = " << 1/p.lub_reduce_parameter << endl;
	// cout << "log_lub_coeff_contact_tan_lubrication = " << log_lub_coeff_contact_tan_total << endl;
	// cout << "log_lub_coeff_contact_tan_dashpot = " << log_lub_coeff_contact_tan_dashpot << endl;
	if (brownian) {
#ifdef DEV
		/* In developing and debugging phases,
		 * we give a seed to generate the same series of random number.
		 * DEV is defined as a preprocessor option in the Makefile
		 */
#ifndef USE_DSFMT
		r_gen = new MTRand(17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
#endif
#ifdef USE_DSFMT
		dsfmt_init_gen_rand(&r_gen, 17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
		dsfmt_init_gen_rand(&rand_gen, 17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
#endif
#endif

#ifndef DEV
#ifndef USE_DSFMT
		r_gen = new MTRand;
#endif
#ifdef USE_DSFMT
		dsfmt_init_gen_rand(&r_gen, hash(std::time(NULL), clock()) ) ; // hash of time and clock trick from MersenneTwister v1.0 by Richard J. Wagner
		dsfmt_init_gen_rand(&rand_gen, hash(std::time(NULL), clock()) ) ; // hash of time and clock trick from MersenneTwister v1.0 by Richard J. Wagner
#endif
#endif
	}
	if (magnetic) {
		if (p.magnetic_type == 1) {
			magnetic_rotation_active = true;
		} else if (p.magnetic_type == 2) {
			magnetic_rotation_active = false;
		} else {
			throw runtime_error("magnetic_type needs to be 1 or 2\n");
		}
	}
	if (p.time_init_relax > 0) {
		time = -p.time_init_relax;
		time_in_simulation_units = -p.time_init_relax*(*ratio_unit_time);
	} else {
		time = 0;
		time_in_simulation_units = 0;
	}
	total_num_timesteps = 0;

	vel_difference.reset();
	if (!p.cross_shear) {
		vel_difference.x = shear_rate*lz;
	} else {
		vel_difference.x = costheta_shear*shear_rate*lz;
		vel_difference.y = sintheta_shear*shear_rate*lz;
	}
	if (p.lubrication_model > 0) {
		stokes_solver.init(np);
	}
	dt = p.dt;
	initializeBoxing();
	checkNewInteraction();
	/* Pre-calculation
	 */
	/* einstein_stress may be affected by sd_coeff.
	 * However, the reason to set value of sd_coeff is not very certain for the moment.
	 * This is why we limit sd_coeff dependence only the diagonal constant.
	 */
	einstein_viscosity = (1+2.5*volume_fraction)/(6*M_PI); // 6M_PI because  6\pi eta_0/T_0 = F_0/L_0^2. In System, stresses are in F_0/L_0^2

	if (p.lubrication_model > 0) {
		resistance_matrix_dblock.resize(np);
		for (int i=0; i<np; i++) {
			resetDBlock(resistance_matrix_dblock[i]);
		}
	}
	for (int i=0; i<np; i++) {
		double FUvalue = p.sd_coeff*radius[i];
		double TWvalue = p.sd_coeff*radius_cubed[i]*4.0/3;
		if (p.lubrication_model == 0) {
			// Stokes drag simulation
			stokesdrag_coeff_f[i] = FUvalue;
			stokesdrag_coeff_f_sqrt[i] = sqrt(FUvalue);
			stokesdrag_coeff_t[i] = TWvalue;
			stokesdrag_coeff_t_sqrt[i] = sqrt(TWvalue);
		} else {
			resistance_matrix_dblock[i].col0[0] = FUvalue;
			resistance_matrix_dblock[i].col1[0] = FUvalue;
			resistance_matrix_dblock[i].col2[0] = FUvalue;
			resistance_matrix_dblock[i].col3[0] = TWvalue;
			resistance_matrix_dblock[i].col4[0] = TWvalue;
			resistance_matrix_dblock[i].col5[0] = TWvalue;
		}
	}
	angle_output = false;
	if (twodimension) {
		if (magnetic == false ||
			magnetic_rotation_active) {
			angle_output = true;
		}
	}
	cout << indent << "Setting up System... [ok]" << endl;
}

void System::initializeBoxing()
{
	/**
		\brief Initialize the boxing system.

		Initialize the BoxSet instance using as a minimal Box size the maximal interaction range between any two particles in the System.
	 */
	double range;
	double max_range = 0;
	for (int i=0; i < np-1; i++) { // N^2 init, sorry :(
		for (int j=i+1; j < np; j++) {
			range = (this->*calcInteractionRange)(i, j);
			if (range > max_range) {
				max_range = range;
			}
		}
	}
	boxset.init(max_range, this);
	for (int i=0; i<np; i++) {
		boxset.box(i);
	}
	boxset.update();
}

void System::timeStepBoxing(const double strain_increment)
{
	/**
		\brief Apply a strain step to the boxing system.
	 */
	if (!zero_shear) {
		shear_strain += strain_increment;
		if (!p.cross_shear) {
			shear_disp.x += strain_increment*lz;
			int m = (int)(shear_disp.x/lx);
			if (shear_disp.x < 0) {
				m--;
			}
			shear_disp.x = shear_disp.x-m*lx;
		} else {
			shear_disp.x += costheta_shear*strain_increment*lz;
			shear_disp.y += sintheta_shear*strain_increment*lz;
			int m = (int)(shear_disp.x/lx);
			if (shear_disp.x < 0) {
				m--;
			}
			shear_disp.x = shear_disp.x-m*lx;
			m = (int)(shear_disp.y/ly);
			if (shear_disp.y < 0) {
				m--;
			}
			shear_disp.y = shear_disp.y-m*ly;
		}
	}
	boxset.update();
}

void System::eventShearJamming(){
	/**
	 \brief Create an event when the shear rate is negative
	*/
	if (shear_rate < 0) {
		Event ev;
		ev.type = "negative_shear_rate";
		events.push_back(Event(ev));
	}
}

void System::timeEvolutionEulersMethod(bool calc_stress,
									   const string& time_or_strain,
									   const double& value_end)
{
	/**
	 \brief One full time step, Euler's method.

	 This method is never used when running a Brownian simulation.
	 */
	in_predictor = true;
	in_corrector = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	setMagneticForceToParticle();
	if (p.lubrication_model == 0) {
		computeVelocitiesStokesDrag();
	} else {
		computeVelocities(calc_stress);
	}
	if (calc_stress) {
		calcStressPerParticle();
	}
	timeStepMove(time_or_strain, value_end);

	if (eventLookUp != NULL){
		(this->*eventLookUp)();
	}

}

/****************************************************************************************************
 ******************************************** Mid-Point Scheme ***************************************
 ****************************************************************************************************/

void System::timeEvolutionPredictorCorrectorMethod(bool calc_stress,
												   const string& time_or_strain,
												   const double& value_end)
{
	/**
	 \brief One full time step, predictor-corrector method.

	 This method is always used when running a Brownian simulation.

	 ### non-Brownian Case

	 Simple mid-point method to solve at dt^2 order
	 \f$ \bm{A}(\bm{U}-\bm{U}_{\infty}) = \bm{F} \f$
	 where \f$\bm{A} \f$ (in Jeffrey notations \cite
	 jeffrey_calculation_1992, that is \f$\bm{R}_{\mathrm{FU}}\f$ in Bossis and Brady
	 \cite brady_stokesian_1988 notations) is the current resistance matrix
	 stored in the stokes_solver, and \f$\bm{F} \f$ are the forces included by the parameter files.

	 - 1st step:
	 - \f$ \bm{U}^{-} = \bm{A}^{-1}( \bm{X}(t) ) \bm{F} ( \bm{X}(t) )  \f$
	 - \f$ \bm{X}' = \bm{X}(t) + \bm{U}^{-}dt \f$

	 - 2nd step:
	 - \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}^{-} ) \bm{F} ( \bm{X}^{-} )  \f$
	 - \f$ \bm{X}(t + dt) = \bm{X}(t) + \frac{1}{2}(\bm{U}^{+}+\bm{U}^{-})dt =  \bm{X}' + \frac{1}{2}(\bm{U}^{+}-\bm{U}^{-})dt \f$

	 ### Brownian Case

	 This routine implements a two-step algorithm for the dynamics with Brownian motion,
	 initially derived by Fixman \cite fixman_1978.
	 The basis of this algorithm is exposed in \cite Ball_1997 and \cite banchio_accelerated_2003.

	 The equation of motion is:
	 \f$ \bm{A}(\bm{U}-\bm{U}_{\infty}) = \bm{F} + \bm{F}_\mathrm{B} + kT \bm{A} \nabla \bm{A}^{-1} \f$
	 with
	 \f$ \langle \bm{F}_\mathrm{B} \rangle = 0\f$
	 and
	 \f$\langle \bm{F}_\mathrm{B} \bm{F}_\mathrm{B}\rangle = \frac{2kT}{dt} \bm{A}\f$,
	 and \f$\bm{F} \f$ are all the other forces included by the parameter files.

	 Reminding that we obtain the Cholesky decomposition \f$ \bm{A} = \bm{L} \bm{L}^T \f$ in the stokes_solver,
	 we obtain X_B in a 2-step algorithm [ B & M 1997 ]:
		- 1st step:
	 + generate a Brownian force \f$ \bm{F}_\mathrm{B}= \sqrt{\frac{2}{dt}} \bm{L} \psi \f$ with \f$ \langle \psi \rangle = 0 \f$ and \f$ \langle \psi \psi \rangle = 1\f$
	 + \f$ \bm{U}^{-} = \bm{A}^{-1}( \bm{X}(t) )( \bm{F}_\mathrm{B} + \bm{F} ( \bm{X}(t) ) )  \f$
	 + \f$ \bm{X}' = \bm{X}(t) + \bm{U}^{-}dt \f$
		- 2nd step:
	 + \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}^{-} ) ( \bm{F}_\mathrm{B} + \bm{F} ( \bm{X}^{-} ) )  \f$ (\b same \f$\bm{F}_\mathrm{B}\f$ as in the first step)
	 + \f$ \bm{X}(t + dt) = \bm{X}(t) + \frac{1}{2}(\bm{U}^{+}+\bm{U}^{-})dt =  \bm{X}' + \frac{1}{2}(\bm{U}^{+}-\bm{U}^{-})dt \f$

	 */

	/* predictor */
	in_predictor = true;
	in_corrector = false;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	setMagneticForceToParticle();
	if (p.lubrication_model > 0) {
		computeVelocities(calc_stress);
	} else {
		computeVelocitiesStokesDrag();
	}
	if (calc_stress) {
		calcStressPerParticle(); // stress compornents
	}
	timeStepMovePredictor(time_or_strain, value_end);
	/* corrector */
	in_predictor = false;
	in_corrector = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	setMagneticForceToParticle();
	if (p.lubrication_model > 0) {
		computeVelocities(calc_stress);
	} else {
		computeVelocitiesStokesDrag();
	}
	if (calc_stress) {
		calcStressPerParticle(); // stress compornents
	}
	if (lowPeclet) {
		// Comupute total stress every time steps for better averaging
		calcStress();
	}
	timeStepMoveCorrector();
}

void System::adaptTimeStep()
{
	/**
	\brief Adapt the time step so that the maximum relative displacement is p.disp_max .
	*/
	if (max_velocity > 0 || max_sliding_velocity > 0) { // small density system can have na_velocity=0
		if (max_velocity > max_sliding_velocity) {
			dt = p.disp_max/max_velocity;
		} else {
			dt = p.disp_max/max_sliding_velocity;
		}
	} else {
		dt = p.disp_max/shear_rate;
	}
}

void System::adaptTimeStep(const string& time_or_strain,
						   const double& value_end)
{
	/**
	\brief Adapt the time step so that (a) the maximum relative displacement is p.disp_max, and (b) time or strain does not get passed the end value.
	*/
	adaptTimeStep();
	if (time_or_strain == "strain") {
		if (fabs(dt*shear_rate) > value_end-fabs(get_shear_strain())) {
			dt = fabs((value_end-fabs(get_shear_strain()))/shear_rate);
		}
	} else {
		if (dt > (value_end-get_time())) {
			dt = (value_end-get_time());
		}
	}
}

void System::timeStepMove(const string& time_or_strain,
						  const double& value_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, Euler method step.
	 */

	/* Adapt dt to get desired p.disp_max	 */
	if (!p.fixed_dt) {
		adaptTimeStep(time_or_strain, value_end);
	}

	time += dt;
	time_in_simulation_units += dt*(*ratio_unit_time);
	total_num_timesteps ++;
	/* evolve PBC */
	double strain_increment = 0;
	if (!zero_shear) {
		strain_increment = dt*shear_rate;
	}

	timeStepBoxing(strain_increment);
	/* move particles */
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (magnetic_rotation_active) {
		for (int i=0; i<np; i++) {
			magnetic_moment[i] += cross(ang_velocity[i], magnetic_moment[i])*dt;
			// @@ Normalize
			// double norm = magnetic_moment[i].norm();
			//magnetic_moment[i] *= magnetic_moment_norm[i]
		}
	}
	if (angle_output) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	checkNewInteraction();
	updateInteractions();
}

void System::timeStepMovePredictor(const string& time_or_strain,
								   const double& value_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, predictor step.
	 */
	if (!brownian) { // adaptative time-step for non-Brownian cases
		if (!p.fixed_dt) {
			adaptTimeStep(time_or_strain, value_end);
		}
	}

	time += dt;
	time_in_simulation_units += dt*(*ratio_unit_time);
	total_num_timesteps ++;
	/* evolve PBC */
	double strain_increment = 0;
	if (!zero_shear) {
		strain_increment = dt*shear_rate;
	}

	/* The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	timeStepBoxing(strain_increment);

	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (angle_output) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	if (magnetic_rotation_active) {
		for (int i=0; i<np; i++) {
			magnetic_moment[i] += cross(ang_velocity[i], magnetic_moment[i])*dt;
		}
	}
	updateInteractions();
	/*
	 * Keep V^{-} to use them in the corrector.
	 */
	for (int i=0; i<np; i++) {
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
}

void System::timeStepMoveCorrector()
{
	/**
	 \brief Moves particle positions according to previously computed velocities, corrector step.
	 */
	for (int i=0; i<np; i++) {
		velocity[i] = 0.5*(velocity[i]+velocity_predictor[i]); // real velocity, in predictor and in corrector
		ang_velocity[i] = 0.5*(ang_velocity[i]+ang_velocity_predictor[i]);
	}
	for (int i=0; i<np; i++) {
		displacement(i, (velocity[i]-velocity_predictor[i])*dt);
	}
	if (angle_output) {
		for (int i=0; i<np; i++) {
			angle[i] += (ang_velocity[i].y-ang_velocity_predictor[i].y)*dt; // no cross_shear in 2d
		}
	}
	if (magnetic_rotation_active) {
		for (int i=0; i<np; i++) {
			magnetic_moment[i] += cross((ang_velocity[i]-ang_velocity_predictor[i]), magnetic_moment[i])*dt;
			// @ Need to be normalized every time step(?)
		}
	}
	checkNewInteraction();
	updateInteractions();
}

bool System::keepRunning(const string& time_or_strain,
						 const double& value_end)
{
	bool keep_running;
	if (time_or_strain == "strain") {
		keep_running = (fabs(get_shear_strain()) < value_end-1e-8) && events.empty();
	} else {
		keep_running = (get_time() < value_end-1e-8) && events.empty();
	}
	return keep_running;
}

void System::timeEvolution(const string& time_or_strain,
						   const double& value_end)
{
	/**
	 \brief Main time evolution routine: evolves the system until time_end

	 This method essentially loops the appropriate one time step
	 method method, according to the Euler vs predictor-corrector or
	 strain rate vs stress controlled choices. On the last time step,
	 the stress is computed.
	 (In the case of low Peclet simulations, the stress is computed at every time step.)
	 r
	 \param time_end Time to reach.
	 */
	static bool firsttime = true;
	in_predictor = false;
	in_corrector = false;
	if (firsttime) {
		double dt_bak = dt; // to avoid stretching contact spring
		dt = 0;
		checkNewInteraction();
		in_predictor = true;
		updateInteractions();
		in_predictor = false;
		dt = dt_bak;
		firsttime = false;
	}
	bool calc_stress = false;
	if (lowPeclet) {
		calc_stress = true;
	}

	while (keepRunning(time_or_strain, value_end)) {
		(this->*timeEvolutionDt)(calc_stress, time_or_strain, value_end); // no stress computation except at low Peclet
	};
	if (events.empty()) {
		calc_stress = true;
		(this->*timeEvolutionDt)(calc_stress, time_or_strain, value_end); // last time step, compute the stress
	}
	if (p.auto_determine_knkt
		&& shear_strain > p.start_adjust){
		adjustContactModelParameters();
	}
}

void System::checkNewInteraction()
{
	/**
	 \brief Checks if there are new pairs of interacting particles. If so, creates and sets up the corresponding Interaction objects.

	 To be called after particle moved.
	 */
	vec3d pos_diff;
	int zshift;
	double sq_dist;
	for (int i=0; i<np-1; i++) {
		for (const int& j : boxset.neighborhood(i)) {
			if (j > i) {
				if (interaction_partners[i].find(j) == interaction_partners[i].end()) {
					pos_diff = position[j]-position[i];
					periodize_diff(pos_diff, zshift);
					sq_dist = pos_diff.sq_norm();
					double scaled_interaction_range = (this->*calcInteractionRange)(i, j);
					double sq_dist_lim = scaled_interaction_range*scaled_interaction_range;
					if (sq_dist < sq_dist_lim) {
						int interaction_new;
						if (deactivated_interaction.empty()) {
							// add an interaction object.
							interaction_new = nb_interaction;
							nb_interaction ++;
						} else {
							// fill a deactivated interaction object.
							interaction_new = deactivated_interaction.front();
							deactivated_interaction.pop();
						}
						// new interaction
						if (nb_interaction >= maxnb_interactionpair) {
							throw runtime_error("Too many interactions.\n"); // @@@ at some point we should lift this limitation
						}
						interaction[interaction_new].activate(i, j, scaled_interaction_range);
					}
				}
			}
		}
	}
	if (magnetic) {
		if (get_time() > time_update_magnetic_pair) {
			updateMagneticPair();
			time_update_magnetic_pair += p.timeinterval_update_magnetic_pair;
		}
	}
}

void System::updateMagneticPair()
{
	for (int i=0; i<np-1; i++) {
		magnetic_pair[i].clear();
	}
	vec3d pos_diff;
	for (int i=0; i<np-1; i++) {
		for (int j=i+1; j<np; j++) {
			pos_diff = position[j]-position[i];
			periodize_diff(pos_diff);
			double sq_dist = pos_diff.sq_norm();
			if (sq_dist < sq_magnetic_interaction_range) {
				magnetic_pair[i].push_back(j);
			}
		}
	}
}

void System::updateInteractions()
{
	/**
	 \brief Updates the state of active interactions.

	 To be called after particle moved.
	 Note that this routine does not look for new interactions (this is done by System::checkNewInteraction), it only updates already known active interactions.
	 It however desactivate interactions removes interactions that became inactive (ie when the distance between particles gets larger than the interaction range).

	 */
	double sq_max_sliding_velocity = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			bool deactivated = false;
			interaction[k].updateState(deactivated);
			if (deactivated) {
				deactivated_interaction.push(k);
			}
			if (interaction[k].is_contact()) {
				double sq_sliding_velocity = interaction[k].relative_surface_velocity.sq_norm();
				if (sq_sliding_velocity > sq_max_sliding_velocity) {
					sq_max_sliding_velocity = sq_sliding_velocity;
				}
			}
		}
	}
	max_sliding_velocity = sqrt(sq_max_sliding_velocity);
	if (magnetic) {
		if (p.magnetic_type == 2) {
			updateMagneticInteractions();
		} else {
			exit(1); // not yet
		}
	}
}

void System::updateMagneticInteractions()
{
	/**

	 Magnetic force is noramlized with

	 \f$ F_M^{ast} = 3*mu_f*m*n/4*pi*(2a)^4 \f$

	 \hat{F}_M = - (16/r**4)*( (m.n)m + (m.n)m + (m.m)n - 5(m.n)(m.n)n)
	 (where r and m are dimensionless distance and dimensionless magnetic moment.)

	 The nondimentinoal force in the themral unit is given by
	 \tilde{F}_M = \hat{F}_M Pe_M
	 */
	magnetic_force_stored.clear();
	vec3d pos_diff;
	vec3d nvec;
	vec3d force0;
	if (p.magnetic_field_type == 0) {
		// External field is vertical.
		// p.magnetic_type needs to be 2.
		for (int p0=0; p0<np-1; p0++) {
			for (const int p1: magnetic_pair[p0]) {
				pos_diff = position[p1]-position[p0];
				periodize_diff(pos_diff);
				double r_sq = pos_diff.sq_norm();
				nvec = pos_diff/sqrt(r_sq);
				double m0_m1 = magnetic_moment[p0].y*magnetic_moment[p1].y;
				force0 = -(16*m0_m1)/(r_sq*r_sq)*nvec;
				force0 *= amplitudes.magnetic; // Peclet number
				magnetic_force_stored.push_back(make_pair(force0, make_pair(p0, p1)));
			}
		}
	} else {
		for (int p0=0; p0<np-1; p0++) {
			for (const int p1: magnetic_pair[p0]) {
				pos_diff = position[p1]-position[p0];
				periodize_diff(pos_diff);
				double r_sq = pos_diff.sq_norm();
				nvec = pos_diff/sqrt(r_sq);
				double m0_n = dot(magnetic_moment[p0], nvec);
				double m1_n = dot(magnetic_moment[p1], nvec);
				double m0_m1 = dot(magnetic_moment[p0], magnetic_moment[p1]);
				force0 = -16*(m1_n*magnetic_moment[p0]+m0_n*magnetic_moment[p1]+(m0_m1-5*m0_n*m1_n)*nvec)/(r_sq*r_sq);
				force0 *= amplitudes.magnetic; // Peclet number
				magnetic_force_stored.push_back(make_pair(force0, make_pair(p0, p1)));
			}
		}
	}
}

void System::buildHydroTerms(bool build_res_mat, bool build_force_GE)
{
	/**
	 \brief Builds the hydrodynamic resistance matrix and hydrodynamic driving force.

	 @param build_res_mat Build the resistance matrix
	 \f$R_{\mathrm{FU}}\f$ (in Bossis and Brady \cite
	 brady_stokesian_1988 notations) and set it as the current
	 resistance in the StokesSolver. If false, the resistance in the
	 StokesSolver is untouched, it is not reset.
	 @param build_force_GE Build the \f$R_{\mathrm{FE}}:E_{\infty}\f$ force
	 and \b add it to the right-hand-side of the StokesSolver
	 */
	if (build_res_mat) {
		// create a new resistance matrix in stokes_solver
		nb_of_active_interactions = nb_interaction-deactivated_interaction.size();
		stokes_solver.resetResistanceMatrix(nb_of_active_interactions,
											resistance_matrix_dblock);
		/* [note]
		 * The resistance matrix is reset with resistance_matrix_dblock,
		 * which is calculated at the beginning.
		 */
		// add GE in the rhs and lubrication terms in the resistance matrix
		(this->*buildLubricationTerms)(true, build_force_GE);
		stokes_solver.completeResistanceMatrix();
	} else {
		// add GE in the rhs
		(this->*buildLubricationTerms)(false, build_force_GE);
	}
}

/* We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
 * This method computes:
 *  - elements of the resistance matrix if 'mat' is true
 *       (only terms diverging as 1/h if lubrication_model == 1, terms in 1/h and log(1/h) for lubrication_model>1 )
 *  - vector Gtilde*Einf if 'rhs' is true (default behavior)
 */
void System::buildLubricationTerms_squeeze(bool mat, bool rhs)
{
	bool shearrate_is_1 = true;
	if (shear_rate != 1) {
		shearrate_is_1 = false;
	}
	for (int i=0; i<np-1; i ++) {
		for (auto& inter : interaction_list[i]){
			int j = inter->partner(i);
			if (j > i) {
				if (inter->lubrication.is_active()) { // Range of interaction can be larger than range of lubrication
					if (mat) {
						const vec3d& nr_vec = inter->nvec;
						inter->lubrication.calcXFunctions();
						stokes_solver.addToDiagBlock(nr_vec, i,
													 inter->lubrication.scaledXA0(), 0, 0, 0);
						stokes_solver.addToDiagBlock(nr_vec, j,
													 inter->lubrication.scaledXA3(), 0, 0, 0);
						stokes_solver.setOffDiagBlock(nr_vec, j,
													  inter->lubrication.scaledXA2(), 0, 0, 0, 0);
					}
					if (rhs) {
						double GEi[3];
						double GEj[3];
						inter->lubrication.calcGE(GEi, GEj); // G*E_\infty term
						if (shearrate_is_1 == false) {
							for (int u=0; u<3; u++) {
								GEi[u] *= shear_rate;
								GEj[u] *= shear_rate;
							}
						}
						stokes_solver.addToRHSForce(i, GEi);
						stokes_solver.addToRHSForce(j, GEj);
					}
				}
			}
		}
		stokes_solver.doneBlocks(i);
	}
	stokes_solver.doneBlocks(np-1);
	stokes_solver.doneBlocks(np);
}

void System::buildLubricationTerms_squeeze_tangential(bool mat, bool rhs)
{
	bool shearrate_is_1 = true;
	if (shear_rate != 1) {
		shearrate_is_1 = false;
	}
	for (int i=0; i<np-1; i ++) {
		for (auto& inter : interaction_list[i]){
			int j = inter->partner(i);
			if (j > i) {
				if (inter->lubrication.is_active()) { // Range of interaction can be larger than range of lubrication
					if (mat) {
						const vec3d& nr_vec = inter->nvec;
						inter->lubrication.calcXYFunctions();
						stokes_solver.addToDiagBlock(nr_vec, i,
													 inter->lubrication.scaledXA0(),
													 inter->lubrication.scaledYA0(),
													 inter->lubrication.scaledYB0(),
													 inter->lubrication.scaledYC0());
						stokes_solver.addToDiagBlock(nr_vec, j,
													 inter->lubrication.scaledXA3(),
													 inter->lubrication.scaledYA3(),
													 inter->lubrication.scaledYB3(),
													 inter->lubrication.scaledYC3());
						stokes_solver.setOffDiagBlock(nr_vec, j,
													  inter->lubrication.scaledXA1(),
													  inter->lubrication.scaledYA1(),
													  inter->lubrication.scaledYB2(),
													  inter->lubrication.scaledYB1(),
													  inter->lubrication.scaledYC1());
					}
					if (rhs) {
						double GEi[3];
						double GEj[3];
						double HEi[3];
						double HEj[3];
						inter->lubrication.calcGEHE(GEi, GEj, HEi, HEj); // G*E_\infty term
						if (shearrate_is_1 == false) {
							for (int u=0; u<3; u++) {
								GEi[u] *= shear_rate;
								GEj[u] *= shear_rate;
								HEi[u] *= shear_rate;
								HEj[u] *= shear_rate;
							}
						}
						stokes_solver.addToRHSForce(i, GEi);
						stokes_solver.addToRHSForce(j, GEj);
						stokes_solver.addToRHSTorque(i, HEi);
						stokes_solver.addToRHSTorque(j, HEj);
					}
				}
			}
		}
		stokes_solver.doneBlocks(i);
	}
}

void System::generateBrownianForces()
{
	/**
	 \brief Generates a Brownian force realization and sets is as the RHS of the stokes_solver.

	 The generated Brownian force \f$F_B\f$ satisfies
	 \f$ \langle F_\mathrm{B} \rangle = 0\f$
	 and
	 \f$\langle F_\mathrm{B} F_\mathrm{B}\rangle = \frac{2kT}{dt} A\f$
	 where \f$A \f$ (in Jeffrey notations \cite
	 jeffrey_calculation_1992, that is \f$R_{\mathrm{FU}}\f$ in Bossis and Brady
	 \cite brady_stokesian_1988 notations) is the current resistance matrix
	 stored in the stokes_solver.

	 Note that it \b sets the rhs of the solver as \f$ rhs = F_B \f$.

	 \f$ F_B\f$ is also stored in sys->brownian_force.
	 */
	double sqrt_2_dt_amp = sqrt(2/dt)*amplitudes.sqrt_temperature;
	for (int i=0; i<linalg_size; i++) {
		brownian_force[i] = sqrt_2_dt_amp*GRANDOM; // \sqrt(2kT/dt) * random vector A (force and torque)
	}
	if (p.lubrication_model > 0) {
		/* L*L^T = RFU
		 */
		stokes_solver.setRHS(brownian_force);
		stokes_solver.compute_LTRHS(brownian_force); // F_B = \sqrt(2kT/dt) * L^T * A
	} else {
		/*
		 *  F_B = \sqrt(2kT/dt) * L^T * A
		 *  U_B = (RFU)^{-1} F_B
		 *  In Stokes drag simulation
		 *  F_B = \sqrt(2kT/dt) * sqrt(RFU) * A
		 *  U_B = (RFU)^{-1} F_B
		 *	    = (RFU)^{-1} \sqrt(2kT/dt) * sqrt(RFU) * A
		 *      = \sqrt(2kT/dt) * A / sqrt(RFU)
		 *  In order to reduce trivial calculations,
		 *  Here, sqrt(RFU) is not included in F_B
		 *  F_B = \sqrt(2kT/dt) * A
		 *  In the function computeBrownianVelocities(),
		 *  U_B = F_B / sqrt(RFU)
		 */
	}
}

void System::setContactForceToParticle()
{
	for (int i=0; i<np; i++) {
		contact_force[i].reset();
		contact_torque[i].reset();
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			interaction[k].contact.addUpContactForceTorque();
		}
	}
}

void System::setRepulsiveForceToParticle()
{
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsive_force[i].reset();
		}
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].repulsion.addUpForce();
			}
		}
	}
}

void System::setMagneticForceToParticle()
{
	if (p.magnetic_type == 1) {
		throw runtime_error("unfinished @ setMagneticForceToParticle\n");
		if (external_magnetic_field.is_zero() ||
			p.magnetic_type == 2) {
			for (int i=0; i<np; i++) {
				magnetic_force[i].reset();
				magnetic_torque[i].reset();
			}
		} else {
			for (int i=0; i<np; i++) {
				magnetic_force[i].reset();
				magnetic_torque[i] = cross(magnetic_moment[i], external_magnetic_field);
			}
		}
		for (const auto& mf : magnetic_force_stored) {
			magnetic_force[mf.second.first] += mf.first;
			magnetic_force[mf.second.second] -= mf.first;
		}
	} else if (p.magnetic_type == 2) {
		for (int i=0; i<np; i++) {
			magnetic_force[i].reset();
		}
		for (const auto& mf : magnetic_force_stored) {
			magnetic_force[mf.second.first] += mf.first;
			magnetic_force[mf.second.second] -= mf.first;
		}
	}
}

void System::buildContactTerms(bool set_or_add)
{
	// sets or adds ( set_or_add = t or f resp) contact forces to the rhs of the stokes_solver.
	if (set_or_add) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, contact_force[i]);
			if (friction) {
				stokes_solver.setRHSTorque(i, contact_torque[i]);
			}
		}
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, contact_force[i]);
			if (friction) {
				stokes_solver.addToRHSTorque(i, contact_torque[i]);
			}
		}
	}
}

void System::buildRepulsiveForceTerms(bool set_or_add)
{
	// sets or adds ( set_or_add = t or f resp) repulsive forces to the rhs of the stokes_solver.
	if (set_or_add) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, repulsive_force[i]);
		}
		stokes_solver.resetRHStorque();
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, repulsive_force[i]);
		}
	}
}

void System::buildMagneticForceTerms(bool set_or_add)
{
	if (p.magnetic_type == 1) {
		if (set_or_add) {
			for (int i=0; i<np; i++) {
				stokes_solver.setRHSForce(i, magnetic_force[i]);
				stokes_solver.setRHSTorque(i, magnetic_torque[i]);
			}
		} else {
			for (int i=0; i<np; i++) {
				stokes_solver.addToRHSForce(i, magnetic_force[i]);
				stokes_solver.addToRHSTorque(i, magnetic_torque[i]);
			}
		}
	} else if (p.magnetic_type == 2) {
		if (set_or_add) {
			for (int i=0; i<np; i++) {
				stokes_solver.setRHSForce(i, magnetic_force[i]);
			}
		} else {
			for (int i=0; i<np; i++) {
				stokes_solver.addToRHSForce(i, magnetic_force[i]);
			}
		}
	}
}

void System::computeMaxNAVelocity()
{
	/**
	 \brief Compute the maximum non-affine velocity

	 Note: it does \b not compute the velocities, just takes the maximum.
	 */
	double sq_max_na_velocity = 0;
	double sq_na_velocity;
	for (int i=0; i<np; i++) {
		sq_na_velocity = na_velocity[i].sq_norm();
		if (sq_na_velocity > sq_max_na_velocity) {
			sq_max_na_velocity = sq_na_velocity;
		}
	}
	//	double sq_na_ang_velocity;
	//	for (int i=0; i<np; i++) {
	//		sq_na_ang_velocity = na_ang_velocity[i].sq_norm()*radius_squared[i];
	//		if (sq_max_na_velocity < sq_na_ang_velocity) {
	//			sq_max_na_velocity = sq_na_ang_velocity;
	//		}
	//	}
	max_velocity = sqrt(sq_max_na_velocity);
}

void System::computeVelocityComponents()
{
	/**
	 \brief Compute velocities component by component.
	 */
	if (!zero_shear) {
		buildHydroTerms(true, true); // build matrix and rhs force GE
	} else {
		buildHydroTerms(true, false); // zero shear-rate
	}
	stokes_solver.solve(vel_hydro, ang_vel_hydro); // get V_H
	buildContactTerms(true); // set rhs = F_C
	stokes_solver.solve(vel_contact, ang_vel_contact); // get V_C
	if (repulsiveforce) {
		buildRepulsiveForceTerms(true); // set rhs = F_repulsive
		stokes_solver.solve(vel_repulsive, ang_vel_repulsive); // get V_repulsive
	}
	if (magnetic) {
		buildMagneticForceTerms(true);
		stokes_solver.solve(vel_magnetic, ang_vel_magnetic); // get V_repulsive
	}
}

void System::computeShearRate()
{
	/**
	 \brief Compute the shear rate under stress control conditions.
	 */
	calcStressPerParticle();
	calcStress();
	double shearstress_con;
	shearstress_con = shearStressComponent(total_contact_stressXF+total_contact_stressGU, p.theta_shear);
	double shearstress_hyd = target_stress-shearstress_con; // the target_stress minus all the other stresses
	double shearstress_rep = 0;
	if (repulsiveforce) {
		shearstress_rep = shearStressComponent(total_repulsive_stressXF+total_repulsive_stressGU, p.theta_shear);
		shearstress_hyd -= shearstress_rep;
	}
	// the total_hydro_stress is computed above with shear_rate=1, so here it is also the viscosity.
	double viscosity_hyd = einstein_viscosity+shearStressComponent(total_hydro_stress, p.theta_shear);
	shear_rate = shearstress_hyd/viscosity_hyd;
	if (shear_strain < init_strain_shear_rate_limit) {
		if (shear_rate > init_shear_rate_limit) {
			shear_rate = init_shear_rate_limit;
		}
	}
	for (int i=0; i<np; i++) {
		vel_hydro[i] *= shear_rate;
		ang_vel_hydro[i] *= shear_rate;
	}
}

void System::computeVelocities(bool divided_velocities)
{
	/**
	 \brief Compute velocities in the current configuration.

	 \param divided_velocities Divide the velocities in components
	 (hydro, contacts, Brownian, ...). (Note that in Brownian
	 simulations the Brownian component is always computed explicitely, independently of the values of divided_velocities.)
	 */
	stokes_solver.resetRHS();
	if (divided_velocities || stress_controlled) {
		if (stress_controlled) {
			shear_rate = 1;
		}
		computeVelocityComponents();
		if (stress_controlled) {
			computeShearRate();
		}
		for (int i=0; i<np; i++) {
			na_velocity[i] = vel_hydro[i]+vel_contact[i];
			na_ang_velocity[i] = ang_vel_hydro[i]+ang_vel_contact[i];
		}
		if (repulsiveforce) {
			for (int i=0; i<np; i++) {
				na_velocity[i] += vel_repulsive[i];
				na_ang_velocity[i] += ang_vel_repulsive[i];
			}
		}
		if (magnetic) {
			for (int i=0; i<np; i++) {
				na_velocity[i] += vel_magnetic[i];
				na_ang_velocity[i] += ang_vel_magnetic[i];
			}
		}
	} else {
		if (!zero_shear) {
			buildHydroTerms(true, true); // build matrix and rhs force GE
		} else {
			buildHydroTerms(true, false); // zero shear-rate
		}
		// for most of the time evolution
		buildContactTerms(false); // add rhs += F_C
		if (repulsiveforce) {
			buildRepulsiveForceTerms(false); // add rhs += F_repulsive
		}
		if (magnetic) {
			buildMagneticForceTerms(false);
		}
		stokes_solver.solve(na_velocity, na_ang_velocity); // get V
		// @@@ I may need to use the force version: stokes_solver.solve(na_velocity) for magnetic work.
	}
	if (brownian) {
		if (in_predictor) {
			/* generate new F_B only in predictor
			 * Resistance matrix is used.
			 */
			generateBrownianForces();
		}
		computeBrownianVelocities();
		for (int i=0; i<np; i++) {
			na_velocity[i] += vel_brownian[i];
			na_ang_velocity[i] += ang_vel_brownian[i];
		}
	}
	/*
	 * The max velocity is used to find dt from max displacement
	 * at each time step.
	 */
	if (!p.fixed_dt && in_predictor) {
		computeMaxNAVelocity();
	}
	adjustVelocitiesLeesEdwardsPeriodicBoundary();
	stokes_solver.solvingIsDone();
}

void System::computeVelocitiesStokesDrag()
{
	/**
	 \brief Compute velocities in Stokes-drag simulation.

	 Note: Velocities of particles are simply proportional to the total forces acting on respective particles.
	 When the contact model includes dashpots, Stokes-drag simulation cannot be used.
	 */
	for (int i=0; i<np; i++) {
		na_velocity[i] = contact_force[i]/stokesdrag_coeff_f[i];
		na_ang_velocity[i] = contact_torque[i]/stokesdrag_coeff_t[i];
	}
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			na_velocity[i] += repulsive_force[i]/stokesdrag_coeff_f[i];
		}
	}
	if (magnetic) {
		if (p.magnetic_type == 1) {
			for (int i=0; i<np; i++) {
				na_velocity[i] += magnetic_force[i]/stokesdrag_coeff_f[i];
				na_ang_velocity[i] += magnetic_torque[i]/stokesdrag_coeff_t[i];
			}
		} else if (p.magnetic_type == 2) {
			for (int i=0; i<np; i++) {
				na_velocity[i] += magnetic_force[i]/stokesdrag_coeff_f[i];
			}
		}
	}
	if (brownian) {
		if (in_predictor) {
			/* generate new F_B only in predictor
			 * Resistance matrix is used.
			 */
			generateBrownianForces();
		}
		computeBrownianVelocities();
		for (int i=0; i<np; i++) {
			na_velocity[i] += vel_brownian[i];
			na_ang_velocity[i] += ang_vel_brownian[i];
		}
	}
	adjustVelocitiesLeesEdwardsPeriodicBoundary();
}

void System::computeBrownianVelocities()
{
	if (p.lubrication_model > 0) {
		stokes_solver.setRHS(brownian_force); // set rhs = F_B (force and torque)
		stokes_solver.solve(vel_brownian, ang_vel_brownian); // get V_B
	} else {
		/* See the comment given in generateBrownianForces()
		 */
		for (int i=0; i<np; i++) {
			int i6 = 6*i;
			vel_brownian[i].x = brownian_force[i6]/stokesdrag_coeff_f_sqrt[i];
			vel_brownian[i].y = brownian_force[i6+1]/stokesdrag_coeff_f_sqrt[i];
			vel_brownian[i].z = brownian_force[i6+2]/stokesdrag_coeff_f_sqrt[i];
			ang_vel_brownian[i].x = brownian_force[i6+3]/stokesdrag_coeff_t_sqrt[i];
			ang_vel_brownian[i].y = brownian_force[i6+4]/stokesdrag_coeff_t_sqrt[i];
			ang_vel_brownian[i].z = brownian_force[i6+5]/stokesdrag_coeff_t_sqrt[i];
		}
	}
	if (twodimension) {
		rushWorkFor2DBrownian();
	}
}

void System::adjustVelocitiesLeesEdwardsPeriodicBoundary()
{
	for (int i=0; i<np; i++) {
		velocity[i] = na_velocity[i];
		ang_velocity[i] = na_ang_velocity[i];
	}
	if (!zero_shear) {
		if (!p.cross_shear) {
			for (int i=0; i<np; i++) {
				velocity[i].x += shear_rate*position[i].z;
				ang_velocity[i].y += 0.5*shear_rate;
			}
			vel_difference.x = shear_rate*lz;
		} else {
			for (int i=0; i<np; i++) {
				velocity[i].x += costheta_shear*shear_rate*position[i].z;
				velocity[i].y += sintheta_shear*shear_rate*position[i].z;
				ang_velocity[i].y += 0.5*costheta_shear*shear_rate;
				ang_velocity[i].x -= 0.5*sintheta_shear*shear_rate;
			}
			vel_difference.x = costheta_shear*shear_rate*lz;
			vel_difference.y = sintheta_shear*shear_rate*lz;
		}
	}
}

void System::rushWorkFor2DBrownian()
{
	/* [note]
	 * Native 2D simulation is not implemented yet.
	 * As a quick implementation, the velocity elements for extra dimension are set to zero
	 */
	if (p.monolayer) {
		/* Particle (3D sphere) cannot move along y-direction.
		 * All other degrees of freedom exist.
		 */
		for (int i=0; i<np; i++) {
			vel_brownian[i].y = 0; // @@ To be checked
		}
	} else {
		/* Particle (2D disk) can rotate only along y-axis.
		 */
		for (int i=0; i<np; i++) {
			vel_brownian[i].y = 0; // @@ To be checked
			ang_vel_brownian[i].x = 0;
			ang_vel_brownian[i].z = 0;
		}
	}
}

void System::displacement(int i, const vec3d& dr)
{
	if (movable[i]) {
		position[i] += dr;
		int z_shift = periodize(position[i]);
		/* Note:
		 * When the position of the particle is periodized,
		 * we need to modify the velocity, which was already evaluated.
		 * The position and velocity will be used to calculate the contact forces.
		 */
		if (z_shift != 0) {
			velocity[i] += z_shift*vel_difference;
		}
	} else {
		velocity[i].reset();
	}
	boxset.box(i);
}

// [0,l]
int System::periodize(vec3d& pos)
{
	/* Lees-Edwards boundary condition
	 *
	 */
	int z_shift = 0;
	if (pos.z >= lz) {
		pos.z -= lz;
		pos -= shear_disp;
		z_shift = -1;
	} else if (pos.z < 0) {
		pos.z += lz;
		pos += shear_disp;
		z_shift = 1;
	}
	while (pos.x >= lx) {
		pos.x -= lx;
	}
	while (pos.x < 0) {
		pos.x += lx;
	}
	if (!twodimension) {
		while (pos.y >= ly) {
			pos.y -= ly;
		}
		while (pos.y < 0) {
			pos.y += ly;
		}
	}
	return z_shift;
}

// periodize + give z_shift= number of boundaries crossed in z-direction
void System::periodize_diff(vec3d& pos_diff, int& zshift)
{
	/* Lees-Edwards boundary condition
	 * The displacement of the second particle along z direction
	 * is zshift * lz;
	 */
	zshift = 0;
	if (pos_diff.z > lz_half) {
		pos_diff.z -= lz;
		pos_diff -= shear_disp;
		zshift = -1;
	} else if (pos_diff.z < -lz_half) {
		pos_diff.z += lz;
		pos_diff += shear_disp;
		zshift = 1;
	}
	while (pos_diff.x > lx_half) {
		pos_diff.x -= lx;
	}
	while (pos_diff.x < -lx_half) {
		pos_diff.x += lx;
	}
	if (!twodimension) {
		while (pos_diff.y > ly_half) {
			pos_diff.y -= ly;
		}
		while (pos_diff.y < -ly_half) {
			pos_diff.y += ly;
		}
	}
}

void System::periodize_diff(vec3d& pos_diff)
{
	/* Used in simple periodic boundary condition
	 * (not Lees-Edwards boundary condition)
	 */
	if (pos_diff.z > lz_half) {
		pos_diff.z -= lz;
	} else if (pos_diff.z < -lz_half) {
		pos_diff.z += lz;
	}
	if (pos_diff.x > lx_half) {
		pos_diff.x -= lx;
	} else if (pos_diff.x < -lx_half) {
		pos_diff.x += lx;
	}
	if (!twodimension) {
		if (pos_diff.y > ly_half) {
			pos_diff.y -= ly;
		} else if (pos_diff.y < -ly_half) {
			pos_diff.y += ly;
		}
	}
}

void System::setSystemVolume(double depth)
{
	string indent = "  System::\t";
	if (twodimension) {
		system_volume = lx*lz*depth;
		cout << indent << "lx = " << lx << " lz = " << lz << " ly = " << depth << endl;
	} else {
		system_volume = lx*ly*lz;
		cout << indent << "lx = " << lx << " lz = " << lz << " ly = " << ly << endl;
	}
}

void System::adjustContactModelParameters()
{
	/**
	 * This method tries to determine
	 * the values of the contact parameters kn, kt
	 * required to reach given maximum normal / tangential spring stretches
	 * (= overlap / tangential displacement ) in steady state.
	 *
	 * The algorithm is a simple adaptative scheme trying to resp. increase or decrease
	 * the spring constants when the stretches are too large or too small.
	 * It works reasonably well away from transitions.
	 * Close to discontinuous shear thickening, one cannot expect satisfying results,
	 * as gigantic stress fluctuations confuse the algorithm.
	 *
	 * The stretch values used are average of maximal values
	 * weighted with an exponential memory kernel:
	 *    \f$k_{avg}(\gamma_n) = C_n \sum_{i=1}^n k(\gamma_i)e^{(\gamma_n-\gamma_i)/\gamma_{avg}} \f$
	 * where \f$\gamma_{avg}\f$ is the user-defined parameter ParameterSet::memory_strain_avg.
	 *
	 * The kn and kt are bounded by user-defined parameters ParameterSet::min_kn, ParameterSet::max_kn, ParameterSet::min_kt, ParameterSet::max_kn.
	 *
	 * The target stretches are given by ParameterSet::overlap_target and ParameterSet::disp_tan_target.
	 *
	 * Additionally, this routine estimates the time step dt.
	 */

	analyzeState();

	double overlap = -min_reduced_gap;
	overlap_avg->update(overlap, shear_strain);
	max_disp_tan_avg->update(max_disp_tan, shear_strain);
	kn_avg->update(p.kn, shear_strain);
	kt_avg->update(p.kt, shear_strain);

	static double previous_shear_strain = 0;
	double deltagamma = (shear_strain-previous_shear_strain);
	double kn_target = kn_avg->get()*overlap_avg->get()/p.overlap_target;
	double dkn = (kn_target-p.kn)*deltagamma/p.memory_strain_k;

	p.kn += dkn;
	if (p.kn < p.min_kn) {
		p.kn = p.min_kn;
	}
	if (p.kn > p.max_kn) {
		p.kn = p.max_kn;
	}
	double kt_target = kt_avg->get()*max_disp_tan_avg->get()/p.disp_tan_target;
	double dkt = (kt_target-p.kt)*deltagamma/p.memory_strain_k;
	p.kt += dkt;
	if (p.kt < p.min_kt) {
		p.kt = p.min_kt;
	}
	if (p.kt > p.max_kt) {
		p.kt = p.max_kt;
	}

	adaptTimeStep();

	previous_shear_strain = shear_strain;

	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() && interaction[k].contact.state>0 ) {
			interaction[k].contact.setSpringConstants();
		}
	}
}

void System::calcLubricationForce()
{
	/*
	 * Calculate lubrication force to output
	 */
	stokes_solver.resetRHS();
	if (!zero_shear) {
		buildHydroTerms(true, true);
	} else {
		buildHydroTerms(true, false); // no GE
	}
	setContactForceToParticle();
	buildContactTerms(false);
	if (repulsiveforce) {
		setRepulsiveForceToParticle();
		buildRepulsiveForceTerms(false);
	}
	stokes_solver.solve(na_velocity, na_ang_velocity);
	stokes_solver.solvingIsDone();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].lubrication.calcLubricationForce();
		}
	}
}

double System::calcInteractionRangeDefault(const int& i, const int& j)
{
	return p.interaction_range*0.5*(radius[i]+radius[j]);
}

double System::calcLubricationRange(const int& i, const int& j)
{
	double rad_ratio = radius[i]/radius[j];
	if (rad_ratio < 2 && rad_ratio > 0.5) {
		return (2+p.lub_max_gap)*0.5*(radius[i]+radius[j]);
	} else {
		double minradius = (radius[i]<radius[j] ? radius[i] : radius[j]);
		return radius[i]+radius[j]+p.lub_max_gap*minradius;
	}
}
