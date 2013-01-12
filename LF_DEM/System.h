//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__
#define CHOLMOD 1
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <string>
//#include <Accelerate/Accelerate.h>
#include "Interaction.h"
#include "cholmod.h"
#include "vec3d.h"
//#include "ContactForce.h"
#include "BrownianForce.h"
#include "BoxSet.h"

using namespace std;
class Simulation;
class Interaction;
class BrownianForce;
class BoxSet;

class System{
private:
	int n3;
	int maxnum_interactionpair;

	queue<int> deactivated_interaction;
	void buildLubricationTerms();
	void buildLubricationTerms_new();

	void buildBrownianTerms();
	void buildContactTerms();
	void addStokesDrag();
	void XA(double iksi, double lambda, double invlambda, double &XAii, double &XAij, double &XAji, double &XAjj);
	void XG(double iksi, double lambda, double invlambda, double &XGii, double &XGij, double &XGji, double &XGjj);

#ifdef CHOLMOD
	cholmod_sparse *sparse_res;
	cholmod_dense *v;
	cholmod_dense *v_nonBrownian, *v_Brownian_init, *v_Brownian_mid; 
	cholmod_dense *rhs_b, *brownian_force;
	int max_lub_int;
	int stype;
	int sorted;
	int packed;
	int xtype;
	vector <int> rows;
	double *diag_values;
	vector <double> *off_diag_values;
	int *ploc;
	void fillSparseResmatrix();
	void allocateSparseResmatrix();
	void addToDiag(double *nvec, int ii, double alpha);
	void appendToColumn(double *nvec, int jj, double alpha);
#else
	double *res;
	int nrhs;
	int *ipiv;
	int lda;
	int ldb;
	int info;
	double *strain_term;
	int lwork;
	double *work;
	char UPLO;
#endif

	BoxSet* boxset;
	void print_res();

	double _lx;
	double _ly;
	double _lz;
	double _lx2; // =lx/2
	double _ly2; // =ly/2
	double _lz2; // =lz/2

protected:
public:
    /* For DEMsystem
     */
	System(){};
	~System();
	int n; // number of particles
	int ts; // time steps
	int dimension;
	vec3d *position;
	double *radius;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *ang_velocity;
	vec3d *force;
	vec3d *torque;
	double **lubstress; // S_xx S_xy S_xz S_yz S_yy
	double **contactstress; // S_xx S_xy S_xz S_yz S_yy
	double mean_lub_stress[5];
	double mean_contact_stress[5];
	double kn;
	double kt;
	double eta;
	double lub_max;
	double sq_lub_max;
	double mu_static; // static friction coefficient.
	double mu_dynamic;// dynamic friction coefficient.
	double dynamic_friction_critical_velocity;
	double sq_critical_velocity;
	bool lubrication;
	bool friction;
	bool brownian;
	double diag_stokes_drag;
	bool poly;
	Interaction *interaction;
	int num_interaction;
	/*
	 * Leading term of lubrication force is 1/h
	 * h = r - 2a.
	 * 1/h when h > h_cutoff
	 * 1/h_cutoff when h <= h_cutoff
	 */
	double h_cutoff;
	BrownianForce *fb;
	/*************************************************************/
	void lx(double length){
	  _lx=length;
	  _lx2=0.5*_lx;
	}
	void ly(double length){
	  _ly=length;
	  _ly2=0.5*_ly;
	}
	void lz(double length){
	  _lz=length;
	  _lz2=0.5*_lz;
	}
	inline double lx(){
	  return _lx;
	}
	inline double ly(){
	  return _ly;
	}
	inline double lz(){
	  return _lz;
	}
	inline double lx2(){
	  return _lx2;
	}
	inline double ly2(){
	  return _ly2;
	}
	inline double lz2(){
	  return _lz2;
	}

	void set_n(int _n){
	  n=_n;
	  n3=3*n;
	}

	double shear_disp;
	double shear_rate;
	double kb_T;
	double volume_fraction;
	double vel_difference;
	double dt;
	double dt_mid, dt_ratio;
	bool draw_rotation_2d;
	vector <int> lubparticle;
	vector <double> lubparticle_vec[3];
	string simu_name;
	
//	void prepareSimulationName();
	void prepareSimulation();
	void allocateRessources();
	void timeEvolution(int time_step);
	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions();

	void calcContactForces();
	double sq_distance(int i, int j);

	double distance(int i, int j);
	double lubricationForceFactor(int i, int j);
	void displacement(int i, const double &dx_, const double &dy, const double &dz);
	void periodize(vec3d*);
	void periodize_diff(vec3d*);
	void periodize_diff(vec3d*, int*);
	void updateVelocity();
	void updateVelocityLubrication();
	void updateVelocityLubricationBrownian();
	void deltaTimeEvolution();
	void forceReset();
	void torqueReset();
	void stressReset();
	void calcStress();
	int numpart(){
		return n;
	}

#ifdef CHOLMOD
	cholmod_factor *L ;
	cholmod_common c ;
#endif

	void lubricationStress(int i, int j);
	void initializeBoxing();


	// interactions
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;


};
#endif /* defined(__LF_DEM__State__) */
