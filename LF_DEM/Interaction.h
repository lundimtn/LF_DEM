//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
//#define RECORD_HISTORY 1

#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "System.h"
#include "StressTensor.h"
using namespace std;
class System;

class Interaction{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	//======= relative position/velocity data  =========//
	double r; // center-center distance
	double a0; // radii
	double a1; // second raddi > a0
	double ro; // ro = a0+a1;
	double ro_half; // = ro/2
	int zshift;
	double gap_nondim; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	double lub_coeff; // = 1/(gap + lub_reduce_parameter)
	vec3d r_vec; // normal vector
	vec3d contact_velocity;
	vec3d disp_tan; // tangential displacement
	vec3d disp_tan_predictor; // tangential displacement
	//===== forces and stresses ==================== //
	double r_lub_max;  // max distance for lubrication
	vec3d lubforce_i; // lubforce_j = - lubforce_i
	double kn_scaled;
	double kt_scaled;
	double colloidalforce_amplitude;
	//===== observables  ========================== //
	double strain_lub_start; // the strain when lubrication object starts.
	double strain_contact_start; // the strain at h=0.
	double duration; // entire lifetime
	double duration_contact; // enture duraction for h < 0
	double max_stress; // Maximum value of stress in the all history of this object.
	int cnt_sliding;  // to count the number of slips.
#ifdef RECORD_HISTORY
	vector <double> gap_history;
	vector <double> overlap_history;
	vector <double> disp_tan_sq_history;
	void outputHistory();
#endif
	/*********************************
	 *       Private Methods         *
	 *********************************/
	
	//======= particles data  ====================//
	double lambda; // a1/a0
	double invlambda; // a0/a1
	
	//======= relative position/velocity  ========//
//	void set_r(const double &_r);
	void setNormalVectorDistanceGap();
	void calcContactVelocity();
	void calcContactVelocity_predictor();
	void calcContactVelocity_corrector();
	
	//======= internal state switches  ===========//
	void activate_contact();
	void deactivate_contact();
	
	//=======   ===========//
	void outputSummary();
	
	//===== forces and stresses computations =====//
	double Fc_normal_norm; // normal contact force
	double F_colloidal_norm;
	vec3d Fc_normal; // normal contact force
	vec3d Fc_tan; // tangential contact force
	vec3d F_colloidal;
	StressTensor colloidal_stresslet_XF; //stress tensor of colloidal force
	StressTensor contact_stresslet_XF_normal; //stress tensor of normal contact force
	StressTensor contact_stresslet_XF_tan; //stress tensor of frictional contact force
	void calcContactInteraction();
	void checkBreakupStaticFriction();
	
protected:
public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
 	Interaction(){};
	~Interaction(){};
	void init(System *sys_);
	//======= state updates  ====================//
	/* Update the follow items:
	 * - r_vec, zshift, _r, and nr_vec
	 * - contact_velocity_tan
	 * - disp_tan
	 * - Fc_normal and Fc_tan
	 * - check breakup of static friction
	 * - State (deactivation, contact)
	 */
	void updateState(bool &deactivated);
	void activate(int i, int j);
	void deactivate();
	
	//======= particles data  ====================//
	int par_num[2];
	inline int
	partner(int i){
		return (i == par_num[0] ? par_num[1] : par_num[0]);
	}
	
	int label;

	inline double get_a0(){
		return a0;
	}
	
	inline double get_a1(){
		return a1;
	}
	
	inline void Ro(double val){
		ro = val;
		ro_half = 0.5*ro;
	}; // ro = a0 + a1

	inline double Ro(){
		return ro;
	}
	//======= relative position/velocity  ========//
	vec3d nr_vec; // vector center to center
	inline double R(){return r;}
	inline double Gap_nondim(){return gap_nondim;}
	
	//======= internal state =====================//
	bool active;
	bool contact;
	//======= Data ===============================//
	
	//	double total_stress_xz;
	//	double stress_xz_integration;
	//=============  Resistance Matrices ====================/
	double XA[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	void GE(double *GEi, double *GEj);
	void calcXA();
	void calcXG();
	void calcXM();
	//===== forces/stresses  ========================== //
	void addUpContactForceTorque();
	void addUpColloidalForce();
	void evaluateLubricationForce();
	double getContactVelocity();
	double getNormalVelocity();
	double getPotentialEnergy();
	inline double getFcNormal(){return Fc_normal_norm;}
	inline vec3d getFcTan(){return Fc_tan;}
	inline double getFcTan_norm(){return Fc_tan.norm();}
	inline double getColloidalForce(){return F_colloidal_norm;}
	inline double disp_tan_norm(){return disp_tan.norm();}
	inline double getLubForce(){return -dot(lubforce_i, nr_vec);}
	//	inline double lubStresslet(int i){return lubstresslet.elm[i];}
	void addHydroStress();
	void addContactStress();
	void addColloidalStress();
	StressTensor getColloidalStressXF(){return colloidal_stresslet_XF;}
	StressTensor getContactStressXF(){return contact_stresslet_XF_normal+contact_stresslet_XF_tan;}
	StressTensor getContactStressXF_normal(){return contact_stresslet_XF_normal;}
	StressTensor getContactStressXF_tan(){return contact_stresslet_XF_tan;}
	void pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
							   StressTensor &stresslet_i, StressTensor &stresslet_j);
	void pairVelocityStresslet(double* &vel_array, StressTensor &stresslet_i, StressTensor &stresslet_j);
	void pairStrainStresslet(StressTensor &stresslet_i, StressTensor &stresslet_j);
	void integrateStress();
	//=========== observables ===============================//
};
#endif /* defined(__LF_DEM__Interaction__) */
