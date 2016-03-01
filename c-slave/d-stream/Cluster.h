#ifndef STM_D_STREAM_CLUSTER_H_
#define STM_D_STREAM_CLUSTER_H_

#include <tuple>
#include <vector>

struct GridTuple
{
	int x;
	int y;
	float density;

	GridTuple() {}

	GridTuple(int x, int y, float density)
		: x(x), y(y), density(density) {}
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

	void set_label(unsigned int label) {
		this->label_ = label;
	}
private:
	unsigned int label_;
	std::vector<GridTuple> grids_;
};

#endif // !STM_D_STREAM_CLUSTER_H_

