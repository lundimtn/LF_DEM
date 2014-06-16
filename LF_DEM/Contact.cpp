//
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Contact.h"
#include "Interaction.h"

Contact::Contact()
{}

Contact::Contact(const Contact& obj)
{
	disp_tan = obj.disp_tan;
	state = obj.state;
}

void
Contact::init(System *sys_, Interaction *interaction_){
	sys = sys_;
	interaction = interaction_;
	state = 0;
	if (sys->friction_model == 1) {
		frictionlaw = &Contact::frictionlaw_standard;
	} else if (sys->friction_model == 2) {
		frictionlaw = &Contact::frictionlaw_criticalload;
	} else if (sys->friction_model == 3) {
		frictionlaw = &Contact::frictionlaw_criticalload_mu_inf;
	} else {
		frictionlaw = &Contact::frictionlaw_null;
	}
}

void
Contact::getInteractionData(){
	interaction->get_par_num(p0, p1);
	double &ro_12 = interaction->ro_12;
	kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
	kt_scaled = ro_12*sys->get_kt(); // F = kt_scaled * disp_tan <-- disp is not scaled
	mu = sys->get_mu_static();
}

/*
 * This is just used in the parameter determination code.
 */
void
Contact::updateContactModel(){
	if (interaction->active) {
		double &ro_12 = interaction->ro_12;
		kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
		kt_scaled = ro_12*sys->get_kt(); // F = kt_s
		if (state > 0) {
			interaction->lubrication.setResistanceCoeff(sys->get_lub_coeff_contact(),
														sys->get_log_lub_coeff_dynamicfriction());
		}
	}
}

void
Contact::activate(){
	// r < a0 + a1
	state = 1;
	disp_tan.reset();
	interaction->lubrication.setResistanceCoeff(sys->get_lub_coeff_contact(),
												sys->get_log_lub_coeff_dynamicfriction());
}

void
Contact::deactivate(){
	// r > a0 + a1
	state = 0;
	disp_tan.reset();
	f_contact_normal_norm = 0;
	f_contact_normal.reset();
	f_contact_tan.reset();
}

/*********************************
 *                                *
 *	   Contact Forces Methods    *
 *                                *
 *********************************/
void
Contact::incrementTangentialDisplacement(){
	if (sys->in_predictor) {
		prev_disp_tan = disp_tan;
	}
	disp_tan = prev_disp_tan+interaction->relative_surface_velocity*sys->get_dt();
}

/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * rvec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Contact::calcContactInteraction(){
	/* gap_nondim is negative, therefore it is allways positive. */
	f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim();
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
	f_contact_tan = kt_scaled*disp_tan;
	(this->*frictionlaw)();
	// if( f_contact_tan.y != 0){
	// 	cout << "weird tan" << endl;
	// 	cout <<  f_contact_tan.x << " " << f_contact_tan.y << " " << f_contact_tan.z << endl;
	// 	cout <<  interaction->nvec.x << " " << interaction->nvec.y << " " << interaction->nvec.z << endl;
	// 	getchar();
	// }
	// if( f_contact_normal.y != 0){
	// 	cout << "weird normal" << endl;
	// 	cout <<  f_contact_normal.x << " " << f_contact_normal.y << " " << f_contact_normal.z << endl;
	// 	cout <<  interaction->nvec.x << " " << interaction->nvec.y << " " << interaction->nvec.z << endl;
	// }
}

void
Contact::frictionlaw_standard(){
	/* [!!!!! NOTE !!!!!]
	 * 
	 * In Brownian case,
	 * "lubforce_normal" significantly affects the friction law.
	 * The difference seems to be subtle in non-Brownian case.
	 *
	 */
	double supportable_tanforce = mu*f_contact_normal_norm;
	double sq_f_tan = f_contact_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // sliding
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		f_contact_tan = kt_scaled*disp_tan;
	} else {
		state = 2; // static friction
	}
	return;
}

void
Contact::frictionlaw_criticalload(){
	/* Since gap_nondim < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 * 
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		supportable_tanforce *= mu;
		double sq_f_tan = f_contact_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			state = 3; // sliding
			disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
			f_contact_tan = kt_scaled*disp_tan;
		} else {
			state = 2; // static friction
		}
	}
	return;
}

void
Contact::frictionlaw_criticalload_mu_inf(){
	/* Since gap_nondim < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		// Never rescaled.
		// [note]
		// The tangential spring constant may be rescaled to control maximum strain
		state = 2; // static friction
	}
	return;
}

void
Contact::frictionlaw_null(){;}

void
Contact::addUpContactForceTorque(){
	if (state > 0) {
		sys->contact_force[p0] += f_contact_normal;
		sys->contact_force[p1] -= f_contact_normal;
		if (state >= 2) {
			sys->contact_force[p0] += f_contact_tan;
			sys->contact_force[p1] -= f_contact_tan;
			vec3d t_ij = cross(interaction->nvec, f_contact_tan);
			sys->contact_torque[p0] += interaction->a0*t_ij;
			sys->contact_torque[p1] += interaction->a1*t_ij;
		}
	}
}

void
Contact::addContactStress(){
	/*
	 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
	 * Fc_normal = -Fc_normal_norm*nvec;
	 * This force acts on particle 1.
	 * stress1 is a0*nvec[*]force.
	 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
	 */
	/* When we compose stress tensor,
	 * even individual leverl, this part calculates symmetric tensor.
	 * This symmetry is expected in the average ensumble.
	 * I'm not sure this is allowed or not.
	 */
	if (state > 0) {
		/*
		 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
		 * Fc_normal = -Fc_normal_norm*nvec;
		 * This force acts on particle 1.
		 * stress1 is a0*nvec[*]force.
		 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
		 * stress1 + stress2 = (a1+a2)*nvec[*]force
		 */
		contact_stresslet_XF_normal.set(interaction->rvec, f_contact_normal);
		contact_stresslet_XF_tan.set(interaction->rvec, f_contact_tan);
	} else {
		contact_stresslet_XF_normal.reset();
		contact_stresslet_XF_tan.reset();
	}
}
