#ifndef __LF_DEM__BoxSet__
#define __LF_DEM__BoxSet__

#include "Box.h"
#include "System.h"
using namespace std;

class BoxSet{
 private:
  double box_xsize;
  double box_ysize;
  double box_zsize;
  
  int x_box_nb;
  int y_box_nb;
  int z_box_nb;
  int box_nb;
  int bottom_box_nb;
  int top_box_nb;
  int bulk_box_nb;
  int topbottom_box_nb;
  bool _is_boxed;

  Box** Boxes;
  Box** BulkBoxes;
  Box** TopBoxes;
  Box** BottomBoxes;
  Box** TopBottomBoxes;

  System *sys;

  int amax, bmax, cmax; // amax = min( x_box_nb, 3), bmax = min( y_box_nb, 3), cmax = min( z_box_nb, 3) 
  

  void updateNeighbors();
  void updateNeighbors(Box*);

  // init methods
  void allocateBoxes();
  void positionBoxes();
  void assignNeighbors();
  void assignNeighborsBulk();
  void assignNeighborsTop();
  void assignNeighborsBottom();
  void assignNeighborsTopBottom();

 protected:
  
 public:

/***** 
	   update()
	   
	   To be called at each time step.
	   It updates the neighborhood relations betwenn boxes.
	   Those relations change at each time step for boxes on top or bottom
*****/
  void update();

/***** 
	   is_boxed()
	   
	   Can be called before calling an other method of BoxSet. 
	   If false, than calls to other method may usually be avoided.
	   They can be performed anyway though, and they are normally safe (and useless)
	   
	   is_boxed() tells if the boxing is effective.
	   If the system size is small, the neighborhood of a box may contain
	   the whole system. If this happens, the boxing consists of only one box (as it 
	   is useless, if not ill defined, to do something else), and is_boxed() returns false.

*****/
  bool is_boxed();

  /*****
	WhichBox(vec3d pos)
	
	return a pointer on the box containg position pos
  *****/
  Box* WhichBox(vec3d);

  BoxSet(double interaction_dist, System *sys_);
  ~BoxSet();

  void printBoxNetwork();

};

#endif /* defined(__LF_DEM__BoxSet__) */
