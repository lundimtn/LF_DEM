#include "Box.h"
using namespace std;


void Box::addNeighborBox(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	neighbors.push_back(neigh_box);
}

void Box::add(unsigned i)
{
	container.insert(i);
}

void Box::remove(unsigned i)
{
	container.erase(i);
}

void Box::buildNeighborhoodContainer()
{
	if (type == 0) {
		return;
	}
	neighborhood_container.clear();
	size_t size = container.size();
	for (const auto& box : neighbors) {
		if (box->type != 0) {
			size += box->getContainer().size();
		}
	}
	neighborhood_container.resize(size);
	int j = 0;
	// own box
	for (const int& k : container) {
		neighborhood_container[j] = k;
		j++;
	}
	for (const auto& box : neighbors) {
		if (box->type != 0) {
			for (const int& k : box->container) {
				neighborhood_container[j] = k;
				j++;
			}
		}
	}
}
