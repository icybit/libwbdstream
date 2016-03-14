#ifndef STM_D_STREAM_CLUSTER_H_
#define STM_D_STREAM_CLUSTER_H_

#include <tuple>
#include <vector>

#ifndef EPSILON
#define EPSILON 0.00001f
#endif EPSILON

struct GridTuple
{
	float x;
	float y;
	// Do we need to know concrete density or is it enough to know whether grid is dense or transitional
	unsigned __int8 density_status;

	GridTuple() {}

	GridTuple(float x, float y, unsigned __int8 density_status)
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

	Cluster(unsigned int label);
	
	void AddElement(GridTuple grid) {
		grids_.push_back(grid);
	}
	
	GridTuple GetElement(int index) {
		return grids_[index];
	}

	void GetPair(int index, float & x, float & y);

	void MergeClusters(Cluster * cluster);

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

	int get_size() { return (int)grids_.size(); }

	unsigned int get_label() { return this->label_; }

	void set_label(unsigned int label) {
		this->label_ = label;
	}
private:
	unsigned int label_;
	std::vector<GridTuple> grids_;
};

#endif // !STM_D_STREAM_CLUSTER_H_

