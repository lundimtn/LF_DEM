//
//  Lubrication.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2014 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Lubrication
 \brief Lubrication interaction, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Lubrication__
#define __LF_DEM__Lubrication__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"

class System;
class Interaction;

class Lubrication{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	Interaction *interaction;
	//======= particles data  ====================//
	unsigned int p0;
	unsigned int p1;
	unsigned int p0_6;
	unsigned int p1_6;
	double range;
	vec3d *nvec;
	double nxnx;
	double nxny;
	double nxnz;
	double nynz;
	double nyny;
	double nznz;
	double nnE;
	double lub_coeff;
	double log_lub_coeff;
	double a0;
	double a1;
	double ro;
	double lambda;
	double invlambda;
	double lambda_square;
	double lambda_cubic;
	double lambda_p_1;
	double lambda_p_1_square;
	double lambda_p_1_cubic;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double YG[4]; // ii ij ji jj
	double YH[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	double YM[4]; // ii ij ji jj
	double cXA[4];
	double cYA[4];
	double cYB[4];
	double cYC[4];
	double cXG[4];
	double cYG[4];
	double cYH[4];
	double cXM[4];
	double cYM[4];
	double ro_12; // = ro/2
	double a0a0_23;
	double a1a1_23;
	double roro_16;
	double a0a0a0_43;
	double a1a1a1_43;
	double rororo_16;
	double a0a0a0_109;
	double a1a1a1_109;
	double rororo_536;
	double g1_XA;
	double g1_inv_XA;
	double g2_YA;
	double g2_inv_YA;
	double g2_YB;
	double g2_inv_YB;
	double g2_YC;
	double g2_inv_YC;
	double g4_YC;
	double g1_XG;
	double g1_inv_XG;
	double g2_YG;
	double g2_inv_YG;
	double g2_YH;
	double g2_inv_YH;
	double g5_YH;
	double g5_inv_YH;
	double g1_XM;
	double g1_inv_XM;
	double g4_XM;
	double g2_YM;
	double g2_inv_YM;
	double g5_YM;

 public:
	Lubrication(Interaction *int_);
	void init(System *sys_);
	bool is_active();
	void getInteractionData();
	void getGeometry();
	void calcLubConstants();
	//===== forces/stresses  ========================== //
    vec3d lubforce_p0; // lubforce_p1 = - lubforce_p0
	void calcPairwiseForce();
	double get_lubforce_normal()
	{
		// positive for compression
        //lubforce_p0.cerr();
		return -dot(lubforce_p0, nvec);
	}
	vec3d get_lubforce_tan()
	{
		return lubforce_p0-dot(lubforce_p0, nvec)*(*nvec);
	}
	vec3d get_lubforce()
	{
		return lubforce_p0;
	}
	void addHydroStress();
	void pairVelocityStresslet(const vec3d& vi, const vec3d& vj,
							   const vec3d& oi, const vec3d& oj,
							   StressTensor& stresslet_i, StressTensor& stresslet_j);
	void pairStrainStresslet(StressTensor& stresslet_i, StressTensor& stresslet_j);
	void updateResistanceCoeff();
	void setResistanceCoeff(double normal_rc, double tangent_rc);
//void setResistanceCoeffTang(double tangent_rc);
	//=============  Resistance Matrices ====================/
	void calcXFunctionsStress();
	void calcXYFunctionsStress();
	std::tuple<vec3d,vec3d> calcGE();
	std::tuple<vec3d,vec3d,vec3d,vec3d> calcGEHE();
	void calcXFunctions();
	void calcXYFunctions();
    // XA
	inline double scaledXA0(){return    a0*XA[0];}
	inline double scaledXA1(){return ro_12*XA[1];}
	inline double scaledXA2(){return ro_12*XA[2];}
	inline double scaledXA3(){return    a1*XA[3];}
    // YA
    inline double scaledYA0(){return    a0*YA[0];}
	inline double scaledYA1(){return ro_12*YA[1];}
	inline double scaledYA2(){return ro_12*YA[2];}
	inline double scaledYA3(){return    a1*YA[3];}
    // YB
	inline double scaledYB0(){return a0a0_23*YB[0];}
	inline double scaledYB1(){return roro_16*YB[1];}
	inline double scaledYB2(){return roro_16*YB[2];}
	inline double scaledYB3(){return a1a1_23*YB[3];}
    // YC
	inline double scaledYC0(){return a0a0a0_43*YC[0];}
	inline double scaledYC1(){return rororo_16*YC[1];}
	inline double scaledYC2(){return rororo_16*YC[2];}
	inline double scaledYC3(){return a1a1a1_43*YC[3];}
    // XG
	inline double scaledXG0(){return a0a0_23*XG[0];}
    inline double scaledXG1(){return roro_16*XG[1];}
    inline double scaledXG2(){return roro_16*XG[2];}
    inline double scaledXG3(){return a1a1_23*XG[3];}
    // YG
    inline double scaledYG0(){return a0a0_23*YG[0];}
	inline double scaledYG1(){return roro_16*YG[1];}
	inline double scaledYG2(){return roro_16*YG[2];}
	inline double scaledYG3(){return a1a1_23*YG[3];}
    // YH
	inline double scaledYH0(){return a0a0a0_43*YH[0];}
	inline double scaledYH1(){return rororo_16*YH[1];}
	inline double scaledYH2(){return rororo_16*YH[2];}
	inline double scaledYH3(){return a1a1a1_43*YH[3];}
    // XM
	inline double scaledXM0(){return a0a0a0_109*XM[0];}
    inline double scaledXM1(){return rororo_536*XM[1];}
    inline double scaledXM2(){return rororo_536*XM[2];}
    inline double scaledXM3(){return a1a1a1_109*XM[3];}
    // YM
	inline double scaledYM0(){return a0a0a0_109*YM[0];}
	inline double scaledYM1(){return rororo_536*YM[1];}
	inline double scaledYM2(){return rororo_536*YM[2];}
	inline double scaledYM3(){return a1a1a1_109*YM[3];}
};
#endif /* defined(__LF_DEM__Lubrication__) */
