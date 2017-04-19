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
#include "DimensionalQty.h"

class OutputData {
public:
	OutputData(): first_time(true), default_precision(6) {}
	~OutputData()
	{
		fout.close();
	}
	void setFile (const std::string& fname,
	              const std::string& data_header,
	              const bool force_to_run);
	void setDefaultPrecision(int precision)
	{
		default_precision = precision;
	}
	void setUnits(Dimensional::UnitSystem units_,
	              Dimensional::Unit::Unit output_unit);
	template<typename T>
	void entryData(std::string name,
	               Dimensional::Dimension dimension,
	               int width,
	               T value,
	               int precision=-1);
	void writeFileHeader();
	void writeColsToFile();
	void writeToFile(std::string header);
	void writeToFile();

private:
	bool first_time;
	Dimensional::Unit::Unit out_unit;
	Dimensional::Unit::Unit internal_unit;
	Dimensional::UnitSystem units;
	std::map <std::string, std::vector<std::string> > output_data;
	std::vector <std::string> insert_order;
	std::vector <std::string> output_data_name;
	std::map <std::string, int> output_data_width;
	std::ofstream fout;
	int default_precision;

	int getLineNumber();
	void initCol(std::string name, int width);
};

template<typename T>
inline void OutputData::entryData(std::string name,
                                  Dimensional::Dimension dimension,
                                  int width,
                                  T value,
                                  int precision)
{
  if (first_time) {
    initCol(name, width);
  }
  int output_precision;
  if (precision > 0) {
    output_precision = precision;
  } else {
    output_precision = default_precision;
  }
  std::ostringstream str_value;
  if (dimension != Dimensional::none) {
    Dimensional::DimensionalQty<T> qty = {dimension, value, internal_unit};
    units.convertFromInternalUnit(qty, out_unit);
    str_value << std::setprecision(output_precision) << qty.value;
  } else {
    str_value << std::setprecision(output_precision) << value;
  }
  output_data[name].push_back(str_value.str());
}

inline void OutputData::setFile(const std::string& fname,
                                const std::string& data_header,
                                const bool force_to_run)
{
	if (force_to_run == false) {
		std::ifstream file_test(fname.c_str());
		if (file_test.good()) {
			file_test.close();
			std::cerr << "The file '" << fname << "' already exists." << std::endl;
			std::cerr << "You need -f option to overwrite." << std::endl;
			exit(1);
		} else {
			file_test.close();
		}
	}
	fout.open(fname.c_str());
	fout << data_header;
}

inline void OutputData::setUnits(Dimensional::UnitSystem units_,
                                 Dimensional::Unit::Unit output_unit)
{
	units = units_;
	auto it = units.getForceTree().cbegin();
	internal_unit = it->second.unit;

	out_unit = output_unit;
}

inline void OutputData::writeToFile(std::string header)
{
	if (first_time) {
		writeFileHeader();
		first_time = false;
	}
	fout << header;
	writeColsToFile();
}

inline void OutputData::writeColsToFile()
{
	int line_nb = getLineNumber();
	for (int i=0; i<line_nb; i++) {
		for (const auto& name : insert_order) {
			const auto &col = output_data[name];
			if (!col.empty()) {
				fout << col[i] << " ";
			} else {
				fout << "n ";
			}
		}
		fout << std::endl;
	}
	for (auto& od : output_data) {
		od.second.clear();
	}
}

inline void OutputData::writeFileHeader()
{
	fout << "# data in " << out_unit << " units." << std::endl;
	int i = 1;
	for (const auto &name : insert_order) {
		int width = output_data_width[name];
		if (width == 1) {
			fout << "#" << i << ": ";
		} else {
			fout << "#" << i << "-" << i+width-1 << ": ";
		}
		i += width;
		fout << name;
		fout << std::endl;
	}
	fout << std::endl;
}

inline void OutputData::writeToFile()
{
	if (first_time) {
		writeFileHeader();
		first_time = false;
	}
	writeColsToFile();
}

inline int OutputData::getLineNumber()
{
	unsigned int line_nb = 0;
	for (const auto& data : output_data) {
		const auto &col = data.second;
		if (line_nb == 0 && col.size() > 0) {
			line_nb = col.size();
		}
		if (col.size() > 0 && col.size() != line_nb) {
			std::cerr << " Error: inconsistent output. Number of lines to output is heterogeneous." << std::endl;
			exit(1);
		}
	}
	return line_nb;
}

inline void OutputData::initCol(std::string name, int width)
{
	if (output_data.find(name) != output_data.end()) {
		return;
	}
	std::vector <std::string> col;
	col.clear();
	output_data[name] = col;
	output_data_name.push_back(name);
	output_data_width[name] = width;
	insert_order.push_back(name);
}

#endif
