//
//  Contact.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Contact
 \brief Contact object, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Contact__
#define __LF_DEM__Contact__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"
using namespace std;
class System;
class Interaction;

class Contact{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	Interaction *interaction;
	void (Contact::*frictionlaw)();
	unsigned short p0;
	unsigned short p1;
	//===== forces and stresses ==================== //
	double kt_scaled;
	double kr_scaled;
	double kn_scaled;
	double mu_static;
	double mu_dynamic;
	double mu_rolling;
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	StressTensor contact_stresslet_XF_normal; //stress tensor of normal contact force
	StressTensor contact_stresslet_XF_tan; //stress tensor of frictional contact force
	vec3d f_contact_normal; // normal contact force
	vec3d f_contact_tan; // tangential contact force
	vec3d f_rolling;
	double ft_max; // friction_model = 5;
	void incrementTangentialDisplacement();
	void incrementRollingDisplacement();

public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	//======= internal state =====================//
	Contact(){};
	Contact(const Contact& obj);
	void init(System *sys_, Interaction *int_);
	void getInteractionData();
	void activate();
	void deactivate();
	vec3d disp_tan; // tangential displacement
	vec3d disp_rolling;
	vec3d prev_disp_tan; // useful for predictor-corrector method: disp_tan in the previous time step
	vec3d prev_disp_rolling;
	vec3d slid_direction;
	void incrementDisplacements();
	int state;
	/* state:
	 * 0 No contact
	 * 1 Friction is not activated (critical load model)
	 * 2 Static friction
	 * 3 Sliding
	 * -2 Switching dynamic to static
	 */
	double f_contact_normal_norm; // normal contact force
	void frictionlaw_criticalload();
	void frictionlaw_criticalload_mu_inf();
	void frictionlaw_standard();
	void frictionlaw_test();
	void frictionlaw_ft_max();
	void frictionlaw_null();
	//===== forces/stresses  ========================== //
	void calcContactInteraction();
	void addUpContactForceTorque();
	inline double get_f_contact_normal_norm()
	{
		return f_contact_normal_norm;
	}
	vec3d get_f_contact_tan()
	{
		return f_contact_tan;
	}
	inline double get_f_contact_tan_norm()
	{
		return f_contact_tan.norm();
	}
	void calcContactStress();
	StressTensor getContactStressXF()
	{
		return contact_stresslet_XF_normal+contact_stresslet_XF_tan;
	}
	StressTensor getContactStressXF_normal()
	{
		return contact_stresslet_XF_normal;
	}
	StressTensor getContactStressXF_tan()
	{
		return contact_stresslet_XF_tan;
	}
};
#endif /* defined(__LF_DEM__Contact__) */
