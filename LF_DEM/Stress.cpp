//
//  Stress.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"

using namespace std;

void System::stressReset()
{
	/**
	   \brief Sets stresses arrays to zero.

	   To be called by System::calcStressPerParticle()
	*/
	for (int i=0; i<np; i++) {
		lubstress[i].reset();
		contactstressGU[i].reset();
	}
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsivestressGU[i].reset();
		}
	}
	if (brownian) {
		for (int i=0; i<np; i++) {
			brownianstressGU[i].reset();
		}
	}
	if (magnetic) {
		for (int i=0; i<np; i++) {
			magneticstressGU[i].reset();
		}
	}
}

void
System::calcStressPerParticle()
{
	/**
	   This method computes the stresses per particle, split by components (hydro, contact, ...).
	   
	   From velocities \f$ V_{\mathrm{I}}\f$ associated with
	   interaction \f$\mathrm{I}\f$, this method gets the stresses \f$ - GV_{\mathrm{I}} \f$. (This corresponds to
	   \f$- GU_{\mathrm{I}} - H\Omega_{\mathrm{I}} \f$ in Jeffrey
	   notations \cite jeffrey_calculation_1992, and
	   \f$- R_{\mathrm{SU}} U_{\mathrm{I}} \f$ in Bossis and Brady 
	   \cite brady_stokesian_1988 notations.)

	   For the hydrodynamic component, it also gets the \f$ M
	   E_{\infty}\f$ term (\f$R_{\mathrm{SE}} E_{\infty}\f$ is B&B
	   notations), so that all in all \f$ S_{\mathrm{H}} =
	   - GV_{\mathrm{H}} + M E_{\infty}\f$.

	   For the point forces, it also gets the \f$ -xF_{\mathrm{I}} \f$ term, so that \f$ S_{\mathrm{I}} =
	   - GV_{\mathrm{I}} - xF_{\mathrm{I}} \f$.

	   For the Brownian forces, it computes (in B&B notations) \f$ S_{\mathrm{B}} =
	   - kT \nabla\dot (R_{\mathrm{SU}}.R_{\mathrm{FU}}^{-1}) \f$ with the mid-step algorithm of Banchio and Brady 
	   \cite banchio_accelerated_2003 with \f$ n=1 \f$.

	   In the Brownian mode, because of the mid-point scheme for the Brownian stress, you
	   should be careful when calling this method from outside of the
	   System::timeEvolutionPredictorCorrectorMethod method.
	*/
	stressReset();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (p.lubrication_model == 1) {
				interaction[k].lubrication.calcXFunctionsStress();
			} else if (p.lubrication_model == 2) {
				interaction[k].lubrication.calcXYFunctionsStress();
			} else {
				cerr << "lubrication_model = " << p.lubrication_model << endl;
				cerr << "lubrication_model = 3 is not implemented" << endl;
				exit(1);
			}
			interaction[k].lubrication.addHydroStress(); // R_SE:Einf-R_SU*v
			interaction[k].contact.calcContactStress(); // - rF_cont
			if (repulsiveforce) {
				interaction[k].repulsion.calcStressXF(); // - rF_rep
			}
		}
	}
	if (brownian) {
		if (in_predictor) {
			for (int i=0; i<np; i++) {
				brownianstressGU_predictor[i] = brownianstressGU[i];
			}
		} else {
			for (int i=0; i<np; i++) {
				/*
				 * [ Banchio & Brady 2003 ] [ Ball & Melrose 1997 ]
				 */
				brownianstressGU[i] = 0.5*(brownianstressGU[i]-brownianstressGU_predictor[i]);
			}
		}
	}
	if (magnetic) {
		vec3d pos_diff;
		vec3d nvec;
		StressTensor magstressXF;
		magneticstressXF.clear();
		for (const auto & mf : magnetic_force_stored) {
			pos_diff = position[mf.second.second]-position[mf.second.first];
			periodize_diff(pos_diff);
			double r = pos_diff.norm();
			nvec = pos_diff/r;
			magstressXF.set(nvec, mf.first);
			magneticstressXF.push_back(magstressXF);
		}
	}
}

void
System::calcStress()
{
	// Lubrication stress
	total_hydro_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
	}
	total_hydro_stress /= system_volume;
	// Stress from contact force
	// GU contribution
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += contactstressGU[i];
	}
	total_contact_stressGU /= system_volume;
	// XF contribution
	total_contact_stressXF_normal.reset();
	total_contact_stressXF_tan.reset();

	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			total_contact_stressXF_normal += interaction[k].contact.getContactStressXF_normal();
			total_contact_stressXF_tan += interaction[k].contact.getContactStressXF_tan();
			if (p.out_particle_stress) {
				StressTensor sc = interaction[k].contact.getContactStressXF();
				unsigned short i,j;
				interaction[k].get_par_num(i,j);
				double r_ij = interaction[k].get_ro();
				contactstressXF[i] += (radius[i]/r_ij)*sc;
				contactstressXF[j] += (radius[j]/r_ij)*sc;
				if(contactstressXF[i].elm[2]>1e6){
					cout << " big value " << endl;
					getchar();
				}
			}
		}
	}
	total_contact_stressXF_normal /= system_volume;
	total_contact_stressXF_tan /= system_volume;
	total_contact_stressXF = total_contact_stressXF_normal+total_contact_stressXF_tan;
	// Stress from repulsive force
	if (repulsiveforce) {
		// XF contribution
		total_repulsive_stressXF.reset();
		for (int k=0; k<nb_interaction; k++) {
			total_repulsive_stressXF += interaction[k].repulsion.getStressXF();
			if (p.out_particle_stress) {
				/* NOTE: 
					As the repulsive force is not a contact force, there is an ambiguity defining the stress per particle. Here we make the choice of attributing 1/2 of the interaction stress to each particle.
				*/
				StressTensor sc = 0.5*interaction[k].repulsion.getStressXF();
				unsigned short i,j;
				interaction[k].get_par_num(i,j);
				repulsivestressXF[i] += sc; 
				repulsivestressXF[j] += sc;
			}
		}
		total_repulsive_stressXF /= system_volume;
		// GU contribution
		total_repulsive_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_repulsive_stressGU += repulsivestressGU[i];
		}
		total_repulsive_stressGU /= system_volume;
	}
	// Stress from Brownian force
	if (brownian) {
		// GU contribution
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += brownianstressGU[i];
		}
		total_brownian_stressGU /= system_volume;
	}
	// Stress from magnetic force
	if (magnetic) {
		// XF contribution
		total_magnetic_stressXF.reset();
		for (const auto &ms : magneticstressXF) {
			total_magnetic_stressXF += ms;
		}
		total_magnetic_stressXF /= system_volume;
		// GU contribution
		total_magnetic_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_magnetic_stressGU += magneticstressGU[i];
		}
		total_magnetic_stressGU /= system_volume;
	}
	//
	total_stress = total_hydro_stress;
	total_stress += total_contact_stressXF;
	total_stress += total_contact_stressGU; // added (Aug 15 2013)
	if (repulsiveforce) {
		total_repulsive_stress = total_repulsive_stressXF+total_repulsive_stressGU;
		total_stress += total_repulsive_stress;
	}
	if (brownian) {
		total_stress += total_brownian_stressGU;
		if (lowPeclet) { // take an averaged stress instead of instantaneous //@@
			stress_avg->update(total_stress, get_time());
			total_stress = stress_avg->get();
		}
	}
	if (magnetic) {
		total_magnetic_stress = total_magnetic_stressXF+total_magnetic_stressGU;
		total_stress += total_magnetic_stress;
	}
	einstein_stress = einstein_viscosity*shear_rate; // should we include that in the hydro_stress definition?
}
