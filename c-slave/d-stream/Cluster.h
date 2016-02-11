#ifndef STM_D_STREAM_CLUSTER_H_
#define STM_D_STREAM_CLUSTER_H_

#include <tuple>
#include <vector>

using namespace std;

class Cluster
{
public:
	Cluster() {}
	Cluster(unsigned int label);

	
	void AddElement(tuple<int, int> grid) {
		grids_.push_back(grid);
	}

	tuple<int, int> GetElement(int index) {
		return grids_[index];
	}

	void GetPair(int index, int & x, int & y);

	void MergeClusters(Cluster * cluster);

	tuple<int, int> PopBack();

	void RemoveElement(int index);

	// Getters

	int get_size() { return grids_.size(); }

	unsigned int get_label() { return this->label_; }

	// Setters

	void set_label(unsigned int label) {
		this->label_ = label;
	}
private:
	unsigned int label_;
	vector<tuple<int, int>> grids_;
};

#endif // !STM_D_STREAM_CLUSTER_H_

