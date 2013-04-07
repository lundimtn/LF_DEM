//
//  GenerateInitConfig.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__GenerateInitConfig__
#define __LF_DEM__GenerateInitConfig__
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "vec3d.h"
#include "System.h"

using namespace std;



class GenerateInitConfig{
private:
	System sys;
	
	int dimension;
	
	double volume_fraction;
	double lx_lz, ly_lz;
	double lx, ly, lz;
	double lx2, ly2, lz2;
	double number_ratio;
	double a1;
	double a2;
	
	vec3d *grad;
	vec3d *prev_grad;
	double gradientDescent();
	double computeGradient();
	void moveAlongGradient(vec3d*, int);
	void storeGradient();
	double step_size;
	MTRand rand_gen;
	
	double zeroTMonteCarloSweep();
	int overlapNumber(int);
	double particleEnergy(int);
	void updateInteractions(int);
	
	vector<vec3d> position;
	vector<double> radius;
	int np, np1, np2;
	vec3d dr;
	
	inline vec3d randUniformSphere(double r);
	inline vec3d randUniformCircle(double r);
	double sqContactDistance(int i, int j, double contact_distance);
	
	void putRandom();
	void setParameters(int argc, const char * argv[]);
	void setSystemParameters();
	void outputPositionData();
	void outputPositionData(vector<vec3d> &positions,
							vector<double> &radii);
public:
	GenerateInitConfig(){};
	int generate(int argc, const char * argv[]);
	
};

#endif /* defined(__LF_DEM__GenerateInitConfig__) */
