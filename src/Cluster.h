#ifndef _DSTREAM_CLUSTER_H_
#define _DSTREAM_CLUSTER_H_

#include <stdint.h>
#include <vector>

#include "Common.h"

struct GridTuple
{
	float x;
	float y;
	// Do we need to know concrete density or is it enough to know whether grid is dense or transitional
	uint8_t density_status;

	GridTuple() {}

	GridTuple(float x, float y, uint8_t density_status)
		: x(x), y(y), density_status(density_status) {}

	bool operator==(const GridTuple & g)
	{
		return (abs(this->x - g.x) < EPSILON && abs(this->y - g.y) < EPSILON);
	}
};

class Cluster
{
public:
	Cluster() {}

	Cluster(uint32_t label);
	
	void AddElement(GridTuple grid) {
		grids_.push_back(grid);
	}
	
	GridTuple GetElement(int index) {
		return grids_[index];
	}

	GridTuple PopBack();

	void RemoveElement(int index);

	void RemoveElement(float x, float y);

	// Getters

	std::vector<GridTuple>::iterator get_begin_iterator() {
		return grids_.begin();
	}

	std::vector<GridTuple>::iterator get_end_iterator() {
		return grids_.end();
	}

	uint32_t get_size() { return (int)grids_.size(); }

	uint32_t get_label() { return this->label_; }

	void set_label(unsigned int label) {
		this->label_ = label;
	}
private:
	uint32_t label_;
	std::vector<GridTuple> grids_;
};

#endif // !_DSTREAM_CLUSTER_H_

