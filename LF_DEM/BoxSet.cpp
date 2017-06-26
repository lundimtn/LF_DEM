#include <stdexcept>
#include <sstream>
#include "BoxSet.h"
#include "System.h"
using namespace std;

BoxSet::BoxSet(double box_min_size,
               const PBC::PeriodicBoundary &pb,
               std::size_t np)
{
	string indent = "  BoxSet::\t";
	cout << indent << "Setting up Cell List System ... ";

	_box_min_size = box_min_size;
	setPeriodicBoundaryConditions(pb);
}

void BoxSet::setupBoxes(const PBC::PeriodicBoundary &pb) {
	_is_boxed = false;
	for(unsigned i=0; i<box_nb.size(), i++) {
		if (box_nb[i] > 3) {
			_is_boxed = true;
		}
	}
	boxes.resize(box_nb[0]*box_nb[1]*box_nb[2]);

	if (!is_boxed()) {
		box_nb = {1,1,1};
	} else {
		assignNeighbors(pb);
	}
	boxMap.resize(np);
	for (unsigned i=0; i<boxMap.size(); i++) {
		boxMap[i] = NULL;
	}
	cout << " [ok]" << endl;
}


std::array<std::size_t, 3> BoxSet::computeBoxNumber() const
{
	std::array<std::size_t, 3> _box_nb;
	const auto &rhomb_depths = rhomb.getDepths();
	for(unsigned i=0; i<_box_nb.size(), i++) {
		_box_nb[i] = (unsigned)(rhomb_depths[i]/_box_min_size);
		if (_box_nb[i] == 0) {
			_box_nb[i] = 1;
		}
	}
	return _box_nb;
}

bool BoxSet::setPeriodicBoundaryConditions(const PBC::PeriodicBoundary &pb)
{
	rhomb = pb.getRhomboid();
	rhomb_unit_normals = rhomb.getNormals();

	auto new_box_nb = computeBoxNumber();
	for (unsigned i=0; i<new_box_nb.size(); i++) {
		if(new_box_nb[i] != box_nb[i]) {
			box_nb = new_box_nb;
			setupBoxes(pb);  // trigger rebox
			return true;
		}
	}
	return false;
}

void BoxSet::assignNeighbors(const PBC::PeriodicBoundary &pb)
{
	auto box_rhomb = rhomb;
	for (unsigned i=0; i<rhomb.edges.size(); i++) {
		box_rhomb.edges[i] /= box_nb[i];
	}

	for (unsigned i=0; i<box_nb[0]; i++) {
		for (unsigned j=0; j<box_nb[1]; j++) {
			for (unsigned k=0; k<box_nb[2]; k++) {
				vec3d pos = (i+0.5)*box_rhomb.edges[0]
							+ (j+0.5)*box_rhomb.edges[1]
							+ (k+0.5)*box_rhomb.edges[2];
				unsigned box_label = i + j*box_nb[1] + k*box_nb[1]*box_nb[2];
				auto &bx = boxes[box_label];
				vec3d delta;
				int m10p1[] = {-1, 0, 1};
				for (const auto& a : m10p1) {
					delta.x = a*box_rhomb.edges[0];
					for (const auto& b : m10p1) {
						delta.y = b*box_rhomb.edges[1];
						for (const auto& c : m10p1) {
							delta.z = c*box_rhomb.edges[2];
							bx.addNeighborBox(whichBox(pb.periodized(pos+delta)));
						}
					}
				}
			}
		}
	}
}

bool BoxSet::is_boxed() const
{
	return _is_boxed;
}

Box* BoxSet::whichBox(const vec3d &pos) const
{
	unsigned idx [3];
	for (unsigned i=0; i<3; i++) {
		if (box_nb[i] > 1) {
			idx[i] = (unsigned)(dot(pos, rhomb_unit_normals[i])*box_nb[i]/rhomb_depths[i]);
		} else {
			idx[i] = 0;
		}
	}
	unsigned label = (idx[0]*box_nb[1]+idx[1])*box_nb[2]+idx[2];
	if (label >= box_labels.size()) {
		ostringstream error_str;
		error_str  << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	return &box_labels[label];
}

Box* BoxSet::whichBox(unsigned int box_label) const
{
	if (box_label >= box_labels.size()) {
		ostringstream error_str;
		throw runtime_error(error_str.str());
	}
	return &box_labels[box_label];
}

void BoxSet::box(unsigned i, vec3d position_i)
{
	Box* b = whichBox(position_i);
	if (b != boxMap[i]) {
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		b->add(i);
		boxMap[i] = b;
	}
}

//public methods
void BoxSet::update()
{
	for (const auto& bx : boxes) {
		bx.buildNeighborhoodContainer();
	}
}

const vector<unsigned>& BoxSet::neighborhood(unsigned i)
{
	return (boxMap[i])->getNeighborhoodContainer();
}

void BoxSet::printBoxNetwork()
{
	for (const auto& bx : Boxes) {
		const auto& neighbors = bx->getNeighborBox();
		for (const auto& neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->getPosition() << " ";
			cerr << neighbor_box->getPosition() << endl;
		}
	}
}

void BoxSet::printBoxContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printNeighborhoodContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getNeighborhoodContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printBoxMap()
{
	for (int i=0; i<sys->get_np(); i++) {
		cerr << i << " " << boxMap[i]->getPosition() << endl;
	}
}

void BoxSet::yaplotBox(std::ofstream &fout_boxing)
{
	vec3d dx(0.5*box_xsize, 0, 0);
	vec3d dz(0, 0, 0.5*box_zsize);
	vec3d dy(0, 0.01, 0);
	fout_boxing << "y 3\n";
	for (const auto& bx : Boxes) {
		//if (bx->is_active()) {
		vec3d center = bx->getPosition();
		vec3d left_bottom = center - dx - dz + dy;
		vec3d right_bottom = center + dx - dz + dy;
		vec3d right_top = center + dx + dz + dy;
		vec3d left_top = center - dx + dz + dy;
		int box_color = bx->type_neighborhood;
		//int box_color = bx->type;
		if (box_color == 0) {
			fout_boxing << "@ 0" << endl;
		} else if (box_color == 1) {
			fout_boxing << "@ 2" << endl;
		} else if (box_color == 2 || box_color == 3) {
			fout_boxing << "@ 3" << endl;
		} else if (box_color == 4 || box_color == 8) {
			fout_boxing << "@ 5" << endl;
		}  else if (box_color == 6) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 7) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 10) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 11) {
			fout_boxing << "@ 6" << endl;
		} else {
			cerr << "unexpected bx type" << endl;
			cerr << box_color << endl;
			exit(1);
		}
		if (sys->twodimension) {
			fout_boxing << "p 4 "<< left_bottom;
			fout_boxing << ' ' << right_bottom;
			fout_boxing << ' ' << right_top;
			fout_boxing << ' ' << left_top;
			fout_boxing << endl;
		} else {
			if (box_color != 0) {
				fout_boxing << "l "<< left_bottom << ' ' << right_bottom << endl;
				fout_boxing << "l " << right_bottom << ' '<< right_top << endl;
				fout_boxing << "l " << right_top << ' '<< left_top << endl;
				fout_boxing << "l " << left_top << ' '<< left_bottom << endl;
			}
		}
	}
	if (0) {
		fout_boxing << "r 1\n";
		for (const auto& bx : Boxes) {
			if (bx->type != 0 && bx->type_neighborhood >= 1) {
				if (bx->type_neighborhood == 1) {
					fout_boxing << "@ 2" << endl;
				} else if (bx->type_neighborhood == 2 || bx->type_neighborhood == 3) {
					fout_boxing << "@ 3" << endl;
				} else if (bx->type_neighborhood == 4 || bx->type_neighborhood == 8) {
					fout_boxing << "@ 4" << endl;
				} else {
					fout_boxing << "@ 0" << endl;
				}
				fout_boxing << "c " <<  bx->getPosition() << endl;
			}
		}
	}
}
