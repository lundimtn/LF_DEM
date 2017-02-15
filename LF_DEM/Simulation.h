/**
 \class Simulation
 \brief Class launching the simulation by setting up the System class and performing predefined shear evolution
 \author Ryohei Seto
 \author Romain Mari
 */
//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012-2016 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <string>
#include <ctime>
#include <map>
#include <set>
#include <algorithm>
#include "global.h"
#include "System.h"
#include "ParameterSet.h"
#include "DimensionalValue.h"
#include "OutputData.h"
#include "Events.h"
#include "Timer.h"

class Simulation
{
private:
	System sys;
	ParameterSet p_initial;
	std::map <std::string, double> force_ratios; // pairs: (force_type_1/force_type_2, force_value_1/force_value_2)
	std::map <std::string, Dimensional::DimensionalValue<double>> dimensional_input_params;
	std::string header_imported_configulation[2];
	ControlVariable::ControlVariable control_var;
	double shear_rate_expectation;
	double strain_end;
	double time_end;
	/*
	 * Resultant data
	 */
	Dimensional::Unit::Unit internal_units;
	Dimensional::Unit::Unit output_units;
	double target_stress_input;
	double input_rate;
	double dimensionless_rate;
	time_t time_strain_0;
	time_t time_strain_1;
	time_t time_strain_end;
	int timestep_1;
	int timestep_end;
	/*
	 * For output data.
	 */
	std::ofstream fout_time;
	std::ofstream fout_input;
	OutputData outdata;
	OutputData outdata_st;
	OutputData outdata_pst;
	OutputData outdata_par;
	OutputData outdata_int;

	/*
	 * For inputs
	 */

	void setupOptionalSimulation(std::string indent);
	std::vector<Sym2Tensor> getParticleStressGroup(std::string group);
	void setConfiguration(bool binary_conf, std::string filename_import_positions);

	Dimensional::UnitSystem units;
public:
	/* For DEMsystem*/
	Simulation();
	~Simulation();
	void simulationSteadyShear(std::string in_args,
	                           std::vector<std::string>& input_files,
	                           bool binary_conf,
	                           ControlVariable::ControlVariable control_variable,
	                           Dimensional::DimensionalValue<double> control_value,
	                           std::string simu_identifier);
	// void simulationfinedSequence(std::string seq_type, std::string in_args, std::vector<std::string> &input_files, bool binary_conf, std::string control_variable);

	void simulationInverseYield(std::string in_args,
	                           std::vector<std::string>& input_files,
	                           bool binary_conf,
	                           ControlVariable::ControlVariable control_variable,
	                           Dimensional::DimensionalValue<double> control_value,
	                           std::string simu_identifier);

	void setupSimulation(std::string in_args,
	                     std::vector<std::string>& input_files,
	                     bool binary_conf,
	                     Dimensional::DimensionalValue<double> control_value,
	                     std::string simu_identifier);
	TimeKeeper initTimeKeeper();
	ParameterSet p;
	bool keepRunning();
	// void timeEvolution(double& next_output_data);
	void generateOutput(const std::set<std::string> &output_events, int& binconf_counter);
	/*********** Events  ************/
	std::list <Event> events;
	void setupEvents();
	void handleEvents();
	System &getSys()
	{
		return sys;
	}

	void assertParameterCompatibility();
	void setDefaultParameters(Dimensional::DimensionalValue<double> control_value);
	void readParameterFile(const std::string& filename_parameters);
	void openOutputFiles(std::string simu_name);
	std::string prepareSimulationName(bool binary_conf,
	                                  const std::string& filename_import_positions,
	                                  const std::string& filename_parameters,
	                                  const std::string& simu_identifier,
	                                  Dimensional::DimensionalValue<double> control_value);
	void echoInputFiles(std::string in_args,
	                    std::vector<std::string>& input_files);
	void autoSetParameters(const std::string& keyword,
	                       const std::string& value);
	void contactForceParameter(std::string filename);
	void contactForceParameterBrownian(std::string filename);
	void importPreSimulationData(std::string filename);
	void resolveTimeOrStrainParameters(const std::map <std::string, Dimensional::DimensionalValue<double>> &);
	std::map<std::string,std::string> getConfMetaData(const std::string &, const std::string &);
	std::string getMetaParameter(std::map<std::string,std::string> &, std::string &, const std::string &);
	std::string getMetaParameter(std::map<std::string,std::string> &, std::string &);
	void exportForceAmplitudes();
	void setLowPeclet();
	Dimensional::Unit::Unit pickInternalUnitsRateControl();
	void setupNonDimensionalization(Dimensional::DimensionalValue<double> control_value);
	/*
	 * For outputs
	 */
	void createDataHeader(std::stringstream& data_header);
	void outputDataHeader(std::ofstream& fout);
	void getSnapshotHeader(std::stringstream& snapshot_header);
	void outputData();
	void outputConfigurationData();
	void outputFinalConfiguration(const std::string&);
	void outputIntFileTxt();
	void outputParFileTxt();
	void outputConfigurationBinary();
	void outputConfigurationBinary(std::string);
	double getRate();
	vec3d shiftUpCoordinate(double x, double y, double z);
	void outputComputationTime();
	bool kill;
	bool force_to_run;
	bool long_file_name;
	bool diminish_output;
	/*********** Events  ************/
	void handleEventsShearJamming();
	void handleEventsFragility();
	std::string gitVersion();
	std::string simu_name;
	void timeEvolutionUntilNextOutput(const TimeKeeper &tk);
	void printProgress();
};
#endif /* defined(__LF_DEM__Simulation__) */
