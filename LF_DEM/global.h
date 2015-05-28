#ifndef __LF_DEM__global__
#define __LF_DEM__global__

#include <string>
#include "vec3d.h"
#ifndef GIT_VERSION
/*
 * GIT_VERSION stores the output of `git describe --dirty --always`
 * run on the last compilation *made in a git repo*.
 *
 * When in a git repo:
 * - Makefile defines GIT_VERSION := $(shell git describe --dirty --always)
 * - In Xcode, VersionInfo.h is automatically generated
 *   by the following script in the Pre-Action of Build.
 *   -----------------------------------
 *   git=/usr/bin/git
 *   cd ${PROJECT_DIR}/LF_DEM
 *   version=`$git describe --dirty`
 *   echo "#define GIT_VERSION \"$version\"" > VersionInfo.h
 *   -----------------------------------
 *   and the code uses VersionInfo.h
 *
 * When NOT in a git repo, the code tries to import VersionInfo.h. 
 * VersionInfo.h is created upon creation of a source tarball with `make tar` from a git repo.
 * 
 */
#include "VersionInfo.h"
#endif

using namespace std;

inline void removeBlank(string &str)
{
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

inline bool getSuffix(const string &str, string &value, string &suffix)
{
	std::size_t suffix_pos = str.find_first_of("abcdefghijklmnopqrstuvwxyz");
	value = str.substr(0, suffix_pos);
	suffix = str.substr(suffix_pos, str.length());
	if (suffix.empty()) {
		return false;
	}
	return true;
}

inline void errorNoSuffix(string quantity)
{
	cerr << "Error : no unit scale (suffix) provided for " << quantity << endl; exit(1);
}

inline bool str2bool(const string &value)
{
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		cerr << "The value should be true or false" << endl;
		exit(1);
	}
}

inline vec3d str2vec3d(const string &value)
{
	string::size_type l1 = value.find("(", 0);
	if (l1 == string::npos) {
		exit(1);
	}
	string::size_type l2 = value.find(",", l1);
	if (l2 == string::npos) {
		exit(1);
	}
	string::size_type l3 = value.find(",", l2+1);
	if (l3 == string::npos) {
		exit(1);
	}
	string::size_type l4 = value.find(")", l3+1);
	if (l4 == string::npos) {
		exit(1);
	}
	double vx = atof(value.substr(l1+1, l2-l1-1).c_str());
	double vy = atof(value.substr(l2+1, l3-l2-1).c_str());
	double vz = atof(value.substr(l3+1, l4-l3-1).c_str());
	return vec3d(vx,vy,vz);
}

inline void Str2KeyValue(const string &str_parameter,
				  string &keyword,
				  string &value)
{
	string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

#endif /* defined(__LF_DEM__global__) */
