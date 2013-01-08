//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Simulation.h"
#include <cmath>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>
Simulation::Simulation(){};
Simulation::~Simulation(){
	if (fout_yap.is_open()){
		fout_yap.close();
	}
	if (fout_rheo.is_open()){
		fout_rheo.close();
	}
};

bool
str2bool(string value){
	if (value == "true"){
		return true;
	} else if (value == "false"){
		return false;
	} else {
		cerr << "The value should be true or false" << endl;
		exit(1);
	}
}

void
Str2KeyValue(string &str_parameter,
			 string &keyword,
			 string &value){
	string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

void
removeBlank(string &str){
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

void
Simulation::AutoSetParameters(const string &keyword,
							  const string &value){
	map<string,int> keylist;
	const int _lub = 1;
	keylist["lub"] = _lub;
	const int _friction = 2;
	keylist["friction"] = _friction;
	const int _brownian = 3;
	keylist["brownian"] = _brownian;
	const int _h_cutoff = 4;
	keylist["h_cutoff"] = _h_cutoff;
	const int _shear_rate = 5;
	keylist["shear_rate"] = _shear_rate;
	const int _kb_T = 6;
	keylist["kb_T"] = _kb_T;
	const int _shear_strain = 7;
	keylist["shear_strain"] = _shear_strain;
	const int _dt = 8;
	keylist["dt"] = _dt;
	const int _dt_ratio = 9;
	keylist["dt_ratio"] = _dt_ratio;
	const int _lub_max = 10;
	keylist["lub_max"] = _lub_max;
	const int _kn = 11;
	keylist["kn"] = _kn;
	const int _kt = 12;
	keylist["kt"] = _kt;
	const int _mu_static = 13;
	keylist["mu_static"] = _mu_static;
	const int _mu_dynamic = 14;
	keylist["mu_dynamic"] = _mu_dynamic;
	const int _dynamic_friction_critical_velocity = 15;
	keylist["dynamic_friction_critical_velocity"] = _dynamic_friction_critical_velocity;
	const int _interval_snapshot = 16;
	keylist["interval_snapshot"] = _interval_snapshot;
	cerr << keyword << ' ' << value  << endl;
	switch(keylist[keyword]){
        case _lub: sys.lub = str2bool(value) ; break;
		case _friction: sys.friction = str2bool(value) ; break;
		case _brownian: sys.brownian = str2bool(value) ; break;
		case _h_cutoff: sys.h_cutoff = atof(value.c_str()); break;
		case _shear_rate: sys.shear_rate = atof(value.c_str()); break;
		case _kb_T: sys.kb_T = atof(value.c_str()); break;
		case _shear_strain: shear_strain = atof(value.c_str()); break;
		case _dt: sys.dt = atof(value.c_str()); break;
		case _dt_ratio: sys.dt_ratio = atof(value.c_str()); break;
		case _lub_max: sys.lub_max = atof(value.c_str()); break;
		case _kn: sys.kn = atof(value.c_str()); break;
		case _kt: sys.kt = atof(value.c_str()); break;
		case _mu_static: sys.mu_static = atof(value.c_str()); break;
		case _mu_dynamic: sys.mu_dynamic = atof(value.c_str()); break;
		case _dynamic_friction_critical_velocity:
			sys.dynamic_friction_critical_velocity = atof(value.c_str()); break;
		case _interval_snapshot: interval_snapshot = atoi(value.c_str()); break;
		default:
			cerr << "The keyword " << keyword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}


void
Simulation::ReadParameterFile(int argc, const char * argv[]){
	filename_import_positions = argv[1];
	filename_parameters = argv[2];
	ifstream fin;
	fin.open(filename_parameters.c_str());
	string keyword;
	string value;
	while (!fin.eof()){
		string line;
		if (!getline( fin, line , ';'))
			break;
		if (fin.eof())
			break;
		string str_parameter;
		removeBlank(line);
		str_parameter = line;
		string::size_type begin_comment;
		string::size_type end_comment;
		do {
			begin_comment = str_parameter.find("/*");
			end_comment = str_parameter.find("*/");
			if (begin_comment > 10000 )
				break;
			str_parameter = str_parameter.substr(end_comment+2);
		}while (true);
		if (begin_comment > end_comment ){
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if( pos_slashslash != string::npos){
			cerr << " // is not syntax for comment out." << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		AutoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void Simulation::SetParametersPostProcess(){
	/* take parameters from import file name.
	 *
	 */
	int i_D = (int)filename_import_positions.find( "D") + 1;
	sys.dimension = atoi( filename_import_positions.substr(i_D, 1).c_str() );
	if (sys.dimension == 2 ){
		// example: D2L10_10vf0.8.dat
		int i_lx = (int)filename_import_positions.find( "L") + 1;
		int j_lx = (int)filename_import_positions.find( "_" );
		int j_lz = (int)filename_import_positions.find( "vf", j_lx);
		int j_vf = (int)filename_import_positions.find( ".dat", j_lz);
		sys.lx = atoi( filename_import_positions.substr(i_lx, j_lx - i_lx).c_str() );
		sys.ly = 0;
		sys.lz = atoi( filename_import_positions.substr(j_lx+1, j_lz - j_lx-1).c_str() );
		sys.volume_fraction = atof( filename_import_positions.substr(j_lz + 2, j_vf - j_lz-2).c_str() );
	} else {
		// example: D3L10_10_10vf0.5.dat
		int i_lx = (int)filename_import_positions.find( "L") + 1;
		int j_lx = (int)filename_import_positions.find( "_", i_lx);
		int j_ly = (int)filename_import_positions.find( "_", j_lx+1);
		int j_lz = (int)filename_import_positions.find( "vf", j_ly+1);
		int j_vf = (int)filename_import_positions.find( ".dat", j_lz);
		sys.lx = atoi( filename_import_positions.substr(i_lx  , j_lx - i_lx).c_str() );
		sys.ly = atoi( filename_import_positions.substr(j_lx+1, j_ly - j_lx-1).c_str() );
		sys.lz = atoi( filename_import_positions.substr(j_ly+1, j_lz - j_ly-1).c_str() );
		sys.volume_fraction = atof( filename_import_positions.substr(j_lz + 2, j_vf-j_lz-2).c_str() );
	}
	
	cerr << "L = " << sys.lx << ' ' << sys.ly << ' ' << sys.lz << endl;
	cerr << "VF = " << sys.volume_fraction << endl;
	/*
	 * Set simulation name and name of output files.
	 */
	sys.prepareSimulationName();
	string yap_filename = "yap_" + sys.simu_name + ".yap";
	string vel_filename = "rheo_" + sys.simu_name + ".dat";
	string vpy_filename = "vpy_" + sys.simu_name + ".dat";
	fout_yap.open(yap_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	fout_vpy.open(vpy_filename.c_str());
	
	sys.sq_lub_max = sys.lub_max*sys.lub_max; // square of lubrication cutoff length.
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	sys.dt_mid = sys.dt/sys.dt_ratio;
	/*
	 * The time steps finishing simulation.
	 */
	ts_max = (int)(shear_strain / sys.dt);
	
	
	/*
	 * For yaplot output data,
	 * rotation of disk (2D) is visualized by cross.
	 */
	if (sys.dimension == 2)
		sys.draw_rotation_2d = true;
	else
		sys.draw_rotation_2d = false;
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
	origin_zero_flow = true;
	/*
	 * The bond width indicates the force strength.
	 */
	yap_force_factor = 0.05;

}

void
Simulation::SetDefaultParameters(int argc, const char * argv[]){
	filename_import_positions = argv[1];
	sys.lub = true;
	sys.brownian = false;
	sys.friction = true;
	sys.h_cutoff = 0.01;
	/*
	 * Simulation parameters
	 */
	/*
	 * Shear rate
	 *
	 */
	sys.shear_rate = 1.0;
	/*
	 * Temperature
	 */
	sys.kb_T = 1.0;
	/*
	 * Simulation terminate at this value
	 */
	shear_strain = 10;
	/*
	 * dt: the time step to integrate the equation of motion.
	 * We need to give a good criterion to give.
	 */
	//	sys.dt = 1e-5;
	sys.dt = 1e-4;
	/*
	 * dt_mid: the intermediate time step for the mid-point 
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	sys.dt_ratio = 100;
	
	/*
	 * Range of lubrication force
	 */
	sys.lub_max = 2.5;
	/*
	 * Contact force parameters
	 *
	 */
	sys.kn = 100; // normal spring constant
	sys.kt = 100; // tangential spring constant
	/*
	 * Particles are spined by the background vorticity.
	 * Small friction coeffient may not stop the sliding between surfaces
	 * of spinning particles, when normal force is small.
	 * We should also estimate the effect of lubrication torque.
	 *
	 */
	sys.mu_static = 1; // static friction coeffient
	sys.mu_dynamic = 0.8; // dynamic friction coeffient
	/*
	 * This is a threshold velocity to swich from dynamic friction to
	 * static friction. But this is a temporal provísional.
	 * There is no reference to give this value.
	 */
	sys.dynamic_friction_critical_velocity = 0.01;
	/*
	 * snapshot for yaplot data.
	 */
	interval_snapshot = 100;

}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open( filename_import_positions.c_str());
	vec3d pos;
	while (true){
		file_import >> pos.x >> pos.y >> pos.z;
		if (file_import.eof())
			break;
		initial_positions.push_back(pos);
	}
	file_import.close();
}

/*
 * Main simulation
 */
void
Simulation::SimulationMain(int argc, const char * argv[]){
	
	if (argc == 2){
		cerr << "Default Parameters" << endl;
		SetDefaultParameters(argc, argv);
	} else {
		cerr << "Read Parameter File" << endl;
		ReadParameterFile(argc, argv);
	}
	SetParametersPostProcess();
	importInitialPositionFile();
	sys.n = (int)initial_positions.size();
	cerr << "N = " << sys.n  << endl;
	sys.prepareSimulation();
	for (int i=0; i < sys.n; i++){
		sys.position[i] = initial_positions[i];
		sys.angle[i] = 0;
	}

	sys.checkNewInteraction();
	//	int count = 0;
	//	while(count++<3){
	double time = 0;
	while(time < ts_max){
		sys.timeEvolution(interval_snapshot);
		outputRheologyData();
		output_yap();
		time += (double)interval_snapshot;
		output_vpython(time);
	}

}

/*
 *
 */
void
Simulation::outputRheologyData(){
	sys.calcStress();
	/*
	 * Output the sum of the normal forces.
	 */
	fout_rheo << sys.dt * sys.ts << ' ';// 1
	fout_rheo << sys.mean_lub_stress[2] + sys.mean_contact_stress[2]  << ' ' ; //2
	fout_rheo << sys.mean_lub_stress[0] << ' ' ; //3
	fout_rheo << sys.mean_lub_stress[1] << ' ' ; //4
	fout_rheo << sys.mean_lub_stress[2] << ' ' ; //5
	fout_rheo << sys.mean_lub_stress[3] << ' ' ; //6
	fout_rheo << sys.mean_lub_stress[4] << ' ' ; //7
	fout_rheo << sys.mean_contact_stress[0] << ' ' ; //8
	fout_rheo << sys.mean_contact_stress[1] << ' ' ; //9
	fout_rheo << sys.mean_contact_stress[2] << ' ' ; //10
	fout_rheo << sys.mean_contact_stress[3] << ' ' ; //11
	fout_rheo << sys.mean_contact_stress[4] << ' ' ; //12
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow){
		z += sys.lz2;
		if (z > sys.lz2){
			x += - sys.shear_disp;
		if ( x < - sys.lx2)
			x += sys.lx;
			z -=  sys.lz;
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::drawLine(char type , vec3d pos, vec3d vec, ofstream &fout){
	fout << type << ' ';
	fout << pos.x << ' '<< pos.y << ' '<< pos.z << ' ';
	fout << pos.x + vec.x << ' '<< pos.y + vec.y << ' '<< pos.z + vec.z << endl;
}

void
Simulation::drawLine2(char type , vec3d pos1, vec3d pos2, ofstream &fout){
	vec3d seg = pos2 - pos1;
	fout << type << ' ';
	if (seg.z > sys.lz2){
		pos2.z -= sys.lz;
		pos2.x -= sys.shear_disp;
		seg = pos2 - pos1;
	} else if (seg.z < -sys.lz2){
		pos2.z += sys.lz;
		pos2.x += sys.shear_disp;
		seg = pos2 - pos1;
	}
		
	while (seg.x > sys.lx2){
		pos2.x -= sys.lx;
		seg = pos2 - pos1;
	}
	while (seg.x < -sys.lx2){
		pos2.x += sys.lx;
		seg = pos2 - pos1;
	}
	
	if (seg.y > sys.ly2){
		pos2.y -= sys.ly;
	} else if (seg.y < -sys.ly2){
		pos2.y += sys.ly;
	}
	
	fout << pos1.x << ' '<< pos1.y << ' '<< pos1.z << ' ';
	fout << pos2.x << ' '<< pos2.y << ' '<< pos2.z << endl;
}

void
Simulation::drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1,
			  ofstream &fout){
	fout << 'l' << ' ';
	fout << x0 << ' ' << y0 << ' ' << z0 << ' ';
	fout << x1 << ' ' << y1 << ' ' << z1 << endl;
}

/* Output data for yaplot visualization.
 *
 */
void
Simulation::output_yap(){
	static bool fasttime = true;
	if (fasttime){
		fasttime = false;
	}else{
		fout_yap << endl;
	}
	/*
	 * yaplot color
	 *
	 * int color_black = 0;
	 * int color_gray = 1;
	 */
	int color_white = 2;
	int color_green = 3;
	int color_yellow = 4;
	int color_orange = 5;
	int color_blue = 6;
	/* Layer 1: Circles for particles
	 */
	fout_yap << "y 1\n";
	fout_yap << "r 1\n";
	fout_yap << "@ " << color_white << endl;
	vec3d pos;
	for (int i=0; i < sys.n; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
								sys.position[i].y - sys.ly2,
								sys.position[i].z - sys.lz2);
		fout_yap << "c " << pos.x << ' ' << pos.y << ' ' << pos.z << endl;
	}

	/* Layer 4: Orientation of particle (2D simulation)
	 */
	
	if (sys.draw_rotation_2d){
		fout_yap << "y 5\n";
		fout_yap << "@ " << color_white << endl;
		for (int i=0; i < sys.n; i++){
			vec3d u(cos(-sys.angle[i]),0,sin(-sys.angle[i]));
			pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
									sys.position[i].y - sys.ly2,
									sys.position[i].z - sys.lz2);
			drawLine('l', pos-u, 2*u, fout_yap);
			u.set(-sin(-sys.angle[i]), 0, cos(-sys.angle[i]));
			drawLine('l', pos-u, 2*u, fout_yap);
		}
	}
	/* Layer 2: Friction
	 */
	if (sys.friction){
		fout_yap << "y 2\n";
		for (int k=0; k < sys.num_interaction; k++){
			if ( sys.interaction[k].contact){
				if (sys.interaction[k].static_friction)
					fout_yap << "@ " << color_green << endl;
				else
					fout_yap << "@ " << color_orange << endl;
				
				int i = sys.interaction[k].particle_num[0];
				pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
										sys.position[i].y - sys.ly2,
										sys.position[i].z - sys.lz2);
				fout_yap << "r " << yap_force_factor*sys.interaction[k].f_tangent.norm()  << endl;
				drawLine('s', pos, -sys.interaction[k].nr_vec, fout_yap);
				int j = sys.interaction[k].particle_num[1];
				pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
										sys.position[j].y - sys.ly2,
										sys.position[j].z - sys.lz2);
				drawLine('s', pos, sys.interaction[k].nr_vec, fout_yap);
			}
		}
	}
	/* Layer 3: Normal
	 * Lubrication + contact force
	 */
	fout_yap << "y 3\n";
	fout_yap << "@ " << color_yellow << endl;
	for (int k=0; k < sys.num_interaction; k++){
		if ( sys.interaction[k].active ){
			int i = sys.interaction[k].particle_num[0];
			int j = sys.interaction[k].particle_num[1];
			fout_yap << "r " << yap_force_factor*sys.interaction[k].valNormalForce() << endl;
			pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
									sys.position[i].y - sys.ly2,
									sys.position[i].z - sys.lz2);
			drawLine('s', pos, -sys.interaction[k].nr_vec, fout_yap);
			pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
									sys.position[j].y - sys.ly2,
									sys.position[j].z - sys.lz2);
			drawLine('s', pos, sys.interaction[k].nr_vec, fout_yap);
		}
	}

	/* Layer 3: Normal
	 * Lubrication + contact force
	 */
//	fout_yap << "y 4\n";
//	fout_yap << "@ " << color_yellow << endl;
//	for (int i=0; i < sys.n; i++){
//		for (int j = i+1; j < sys.n; j++){
//			double f_ij = 0;
//			if (sys.lub == true){
//				f_ij += sys.lubricationForceFactor(i, j);
//			}
//			if (contact_pair[i][j] != -1){
//				f_ij += -fc[contact_pair[i][j]].f_normal;
//			}
//			if (f_ij > 0){
//				fout_yap << "r " << yap_force_factor*f_ij << endl;
//				vec3d pos1 = shiftUpCoordinate(sys.position[i].x - sys.lx2,
//											   sys.position[i].y - sys.ly2,
//											   sys.position[i].z - sys.lz2);
//				vec3d pos2 = shiftUpCoordinate(sys.position[j].x - sys.lx2,
//											   sys.position[j].y - sys.ly2,
//											   sys.position[j].z - sys.lz2);
//				
//				drawLine2('s', pos1, pos2, fout_yap);
//			}
//		}
//	}

		
	/* Layer 6: Box and guide lines
	 */
	fout_yap << "y 6\n";
	fout_yap << "@ " << color_blue << endl;
	if (sys.dimension == 2){
		drawLine(-sys.lx2, 0, -sys.lz2,  sys.lx2, 0, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, 0,  sys.lz2, -sys.lx2, 0, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, 0,  sys.lz2,  sys.lx2, 0,  sys.lz2, fout_yap);
		drawLine( sys.lx2, 0, -sys.lz2,  sys.lx2, 0,  sys.lz2, fout_yap);
		/*
		for (int i = 0 ; i < 10; i ++){
			drawLine(-sys.lx2                , 0, -sys.lz2 + i*0.2*sys.lz2,
					 -sys.lx2 + i*0.2*sys.lz2, 0,                 -sys.lz2, fout_yap);
			drawLine( sys.lx2                   ,        0, sys.lz2 - i*0.2*sys.lz2,
					  sys.lx2 - i*0.2*sys.lz2, 0,  sys.lz2, fout_yap);
		}
		drawLine(-sys.lx2, 0, sys.lz2, sys.lx2 , 0, -sys.lz2, fout_yap);
		 */
	} else {
		drawLine(-sys.lx2, -sys.ly2, -sys.lz2,  sys.lx2, -sys.ly2, -sys.lz2, fout_yap);
		drawLine(-sys.lx2,  sys.ly2, -sys.lz2,  sys.lx2,  sys.ly2, -sys.lz2, fout_yap);

		drawLine(-sys.lx2,  sys.ly2,  sys.lz2,  sys.lx2,  sys.ly2,  sys.lz2, fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, -sys.ly2,  sys.lz2, fout_yap);
		
		drawLine(-sys.lx2, -sys.ly2, -sys.lz2, -sys.lx2, sys.ly2, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2, -sys.lx2, sys.ly2,  sys.lz2, fout_yap);
		drawLine( sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, sys.ly2,  sys.lz2, fout_yap);
		drawLine( sys.lx2, -sys.ly2, -sys.lz2,  sys.lx2, sys.ly2, -sys.lz2, fout_yap);
		
		drawLine( sys.lx2,  sys.ly2,  sys.lz2,  sys.lx2, sys.ly2, -sys.lz2,  fout_yap);
		drawLine(-sys.lx2,  sys.ly2,  sys.lz2, -sys.lx2, sys.ly2, -sys.lz2,  fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2, -sys.lx2, -sys.ly2, -sys.lz2,  fout_yap);
		drawLine( sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, -sys.ly2, -sys.lz2,  fout_yap);
	}
}



/* Output data for vpython visualization.
 *
 */
void
Simulation::output_vpython(double time){
	vec3d pos;
	fout_vpy << "time: " << time << endl;
	for (int i=0; i < sys.n; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
								sys.position[i].y - sys.ly2,
								sys.position[i].z - sys.lz2);
		fout_vpy << i << ' ' << pos.x << ' ' << pos.y << ' ' << pos.z << endl;
	}
	fout_vpy << endl;
}
