#ifndef STM_D_STREAM_CLUSTER_H_
#define STM_D_STREAM_CLUSTER_H_

#include <tuple>
#include <vector>

using namespace std;

struct GridTuple
{
	int x;
	int y;
	float density;
	unsigned __int8 neighbors;

	GridTuple() {}

	GridTuple(int x, int y, float density, unsigned __int8 neighbors)
		: x(x), y(y), density(density), neighbors(neighbors) {}
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

	void GetPair(int index, int & x, int & y);

	void MergeClusters(Cluster * cluster);

	GridTuple PopBack();

	void RemoveElement(int index);

	void RemoveElement(int x, int y);

	// Getters

	int get_size() { return grids_.size(); }

	unsigned int get_label() { return this->label_; }

	unsigned __int8 get_neighbors_count(int index) {
		return grids_[index].neighbors;
	}

	// Setters
	void increase_neighbors(int index) {
		grids_[index].neighbors++;
	}

	void decrease_neighbors(int index) {
		grids_[index].neighbors--;
	}

	void set_label(unsigned int label) {
		this->label_ = label;
	}
private:
	unsigned int label_;
	vector<GridTuple> grids_;
};

#endif // !STM_D_STREAM_CLUSTER_H_

