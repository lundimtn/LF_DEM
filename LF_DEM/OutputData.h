//
//  OutputData.h
//  LF_DEM
//
//  Created by Ryohei Seto on 6/1/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class OutputData
 \brief Utility class to output data line-by-line, primitively dimension-aware.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_OutputData_h
#define LF_DEM_OutputData_h

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

class OutputData {
private:
	int number_of_data;
	bool first_time;
	vector<string> output_data;
	vector<string> output_data_name;
	vector<string> output_data_type;
	std::map <string, double> converter;
	 
public:
	OutputData():
	first_time(true) {}

	void init(int number_of_data_) {
		if (first_time) {
			number_of_data = number_of_data_;
			output_data.resize(number_of_data);
			output_data_name.resize(number_of_data);
			output_data_type.resize(number_of_data);
			for (string &od : output_data) {
				od = "n";
			}
			for (string &odn : output_data_name) {
				odn = "blank";
			}
			for (string &odt : output_data_type) {
				odt = "none";
			}
		}
	}
	
	void setDimensionlessNumber(double dimensionless_number)
	// dimensionless_number = internal_force_unit/output_force_unit
	{
		converter["none"] = 1;
		converter["viscosity"] = 6*M_PI;
		converter["stress"] = dimensionless_number;
		converter["time"] = 1/dimensionless_number;
		converter["rate"] = dimensionless_number;
		converter["velocity"] = dimensionless_number;
	}
	
	template<typename T>
	void entryData(int num, 
				string name, 
				string type, 
				T value)
	{
		int index = num-1;
		ostringstream str_value;
		str_value << converter[output_data_type[index]]*value;
		if (first_time) {
			if (output_data_name[index] != "blank") {
				cerr << "data["<< index << "] is redefined." << endl;
				exit(1);
			} else {
				output_data_name[index] = name;
			}
		}
		output_data_type[index] = type;
		output_data[index] = str_value.str();
	}
	
	void exportFile(ofstream &fout_data)
	{

		if (first_time) {
			for (int i=0; i<number_of_data; i++) {
				fout_data << "#" << i+1 << ": ";
				fout_data << output_data_name[i];
				fout_data << endl;
			}
			first_time = false;
		}
		
		for (int i=0; i<number_of_data; i++) {
			fout_data << output_data[i] << ' ';
		}
		fout_data << endl;
	}
};
#endif
