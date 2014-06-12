//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"


void
System::calcStressPerParticle(){
	/////////////////////////////////////////////////
	// from the velocities V_H, V_C, V_Coll, V_B, 
	// compute stresses R_SV * V
	// and then add rF stress for F_C and F_Coll
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (lubrication_model == 1) {
				interaction[k].lubrication.calcXFunctionsStress();
			} else if (lubrication_model == 2) {
				interaction[k].lubrication.calcXYFunctionsStress();
			} else if (lubrication_model == 3) {
				if (interaction[k].is_contact()){
					interaction[k].lubrication.calcXYFunctionsStress();
				} else {
					interaction[k].lubrication.calcXFunctionsStress();
				}
			}
			interaction[k].lubrication.addHydroStress(); // - R_SU * v
			interaction[k].contact.addContactStress(); //  - rF_cont
			interaction[k].addColloidalStress(); //  - rF_colloid
		}
	}
}



void
System::calcStress(){
	//	calcStressesHydroContact();
	total_hydro_stress.reset();
	total_contact_stressGU.reset();
	total_colloidal_stressGU.reset();
	total_contact_stressXF_normal.reset();
	total_contact_stressXF_tan.reset();
	total_colloidal_stressXF.reset();
	if (brownian) {
		total_brownian_stressGU.reset();
	}	
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
		total_contact_stressGU += contactstressGU[i];
		total_colloidal_stressGU += colloidalstressGU[i];
		if (brownian) {
			total_brownian_stressGU += brownianstressGU[i];
		}
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			total_colloidal_stressXF += interaction[k].getColloidalStressXF();
		}
		if (interaction[k].is_contact()) {
			total_contact_stressXF_normal += interaction[k].contact.getContactStressXF_normal();
			total_contact_stressXF_tan += interaction[k].contact.getContactStressXF_tan();
		}
	}
	total_hydro_stress /= System_volume();
	total_contact_stressGU /= System_volume();
	total_contact_stressXF_normal /= System_volume();
	total_contact_stressXF_tan /= System_volume();
	total_colloidal_stressGU /= System_volume();
	total_colloidal_stressXF /= System_volume();
	if (brownian) {
		total_brownian_stressGU /= System_volume();
	}
}
