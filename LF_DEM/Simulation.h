//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include "System.h"

class Simulation{
private:
	System sys;
	string filename_;
	bool kn_kt_adjustment;
	double shear_strain_end;
	vector<vec3d> initial_position;
	vector<double> radius;
	string filename_import_positions;
	string filename_parameters;
	double strain_interval_output_data;
	double strain_interval_output;
	/*
	 * Resultant data
	 */
	double viscosity;
	double normalstress_diff_1;
	double normalstress_diff_2;
	double viscosity_hydro;
	double normalstress_diff_1_hydro;
	double normalstress_diff_2_hydro;
	double viscosity_cont_XF;
	double normalstress_diff_1_cont_XF;
	double normalstress_diff_2_cont_XF;
	double viscosity_friction; // Fc_tan contribution.
	double normalstress_diff_1_friction;
	double normalstress_diff_2_friction;
	double viscosity_cont_GU;
	double normalstress_diff_1_cont_GU;
	double normalstress_diff_2_cont_GU;
	double viscosity_col_XF;
	double normalstress_diff_1_col_XF;
	double normalstress_diff_2_col_XF;
	double viscosity_col_GU;
	double normalstress_diff_1_col_GU;
	double normalstress_diff_2_col_GU;
	double viscosity_brownian;
	double normalstress_diff_1_brownian;
	double normalstress_diff_2_brownian;
	/*
	 * For output data.
	 */
	ofstream fout_rheo;
	ofstream fout_particle;
	ofstream fout_interaction;
	bool out_data_particle;
	bool out_data_interaction;
	bool origin_zero_flow;
	/*
	 */
	void timeEvolution();
	void evaluateData();
	/*
	 * For inputs
	 */
	void setDefaultParameters();
	void readParameterFile();
	void openOutputFiles();
	void prepareSimulationName();
	void autoSetParameters(const string &keyword,
						   const string &value);
	void importInitialPositionFile();
	/*
	 * For outputs
	 */
	void outputDataHeader(ofstream &fout);
	void outputRheologyData();
	void outputConfigurationData();
	vec3d shiftUpCoordinate(double x, double y, double z);
public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void SimulationMain(int argc, const char * argv[]);
	void RelaxationZeroShear(vector<vec3d> &position_,
							  vector<double> &radius_,
							  double lx_, double ly_, double lz_, double volume_fraction_);
};
#endif /* defined(__LF_DEM__Simulation__) */
