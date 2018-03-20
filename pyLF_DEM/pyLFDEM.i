%module pyLFDEM
%include <typemaps.i>
%apply double &INPUT { double &next_output_data };
%apply double &INPUT { double &next_output_config };
%apply int &INPUT { int &binconf_counter };
%{
/* Put headers and other declarations here */
#include "../LF_DEM/Simulation.h"
#include "../LF_DEM/SystemHelperFunctions.h"
#include "../LF_DEM/Timer.h"
#include "../LF_DEM/global.h"
%}
%include "../LF_DEM/vec3d.h"
%include "../LF_DEM/Sym2Tensor.h"
%include "../LF_DEM/DimensionalQty.h"
%include "../LF_DEM/Configuration.h"
%include <std_string.i>
%include <std_set.i>
%include <std_vector.i>
%include <std_pair.i>

// Instantiate templates used by example
namespace std {
  %template(StringVector) vector<string>;
  %template(Vec3dVector) vector<vec3d>;
  %template(DoubleVector) vector<double>;
  %template(Sym2Vector) vector<Sym2Tensor>;
  %template(PairUIntUInt) pair<unsigned int, unsigned int>;
  %template(PairIntInt) pair<int, int>;
  %template(PairDoubleString) pair<double, string>;
  %template(SetString) set<string>;
}
%template(DoubleDimQty) Dimensional::DimensionalQty<double>;
// imperfect error handling. Errors get passed to Python, but sometimes still seg fault
%include "exception.i"
%exception {
    try {
        $action
    } catch(const std::exception& e) {
	  SWIG_exception(SWIG_RuntimeError, e.what());
    }
}
%include "../LF_DEM/System.h"
%include "../LF_DEM/ParameterSet.h"
%include "../LF_DEM/ParameterSetFactory.h"
%include "../LF_DEM/SystemHelperFunctions.h"
%include "../LF_DEM/Simulation.h"
%include "../LF_DEM/global.h"
%include "../LF_DEM/Timer.h"
