//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "System.h"
using namespace std;
class System;

class Interaction{
private:
	System *sys;
	vec3d contact_velocity;
	vec3d contact_velocity_tan;
	vec3d unit_contact_velocity_tan;
	vec3d xi; // tangential displacement

	void calcDistanceNormalVector();
	void assignDistanceNormalVector(vec3d, double, int);
	void calcContactVelocity();

	void incrementContactTangentialDisplacement();
//	vec3d t_tangent;
	int pd_z;
	double sqnorm_contact_velocity;

	
	// state switch
	void deactivate();
	void activate_contact();
	void deactivate_contact();
protected:
	void calcStaticFriction();
	void calcDynamicFriction();
	void calcNormalVector();
public:
 	Interaction(){};
	~Interaction(){};
	void init(System *sys_);

	void activate(int i, int j, vec3d pos_diff, double distance, int zshift);
	bool update(); // after particles dispacement
	bool active;
	int particle_num[2];
	double r; // center-center distance // done
	double a0, a1;
	double ro; // ro = a0 + a1
	double lambda, invlambda;  // a1/a0 , a0/a1
	vec3d r_vec; // vector center to center
	vec3d nr_vec; // normal vector
	int partner(int);


	bool contact;
	bool static_friction;

	double f_normal;
	vec3d f_tangent;

	double valNormalForce();
	void calcContactInteraction();
	void calcContactInteractionNoFriction();
	void calcContactStress();
	void addLubricationStress();
	void addContactStress();

};

#endif /* defined(__LF_DEM__Interaction__) */
