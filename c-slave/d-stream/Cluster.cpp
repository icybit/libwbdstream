#include "Cluster.h"


Cluster::Cluster(unsigned int label)
{
	this->label_ = label;
}

void Cluster::GetPair(int index, int & x, int & y)
{
	x = get<0>(grids_[index]);
	y = get<1>(grids_[index]);
}

void Cluster::MergeClusters(Cluster * cluster)
{
	for (int i = 0; i < cluster->get_size(); i++)
	{
		this->AddElement(cluster->PopBack());
	}
}

tuple<int, int> Cluster::PopBack()
{
	tuple<int, int> last = *(grids_.end() - 1);
	grids_.pop_back();
	return last;
}

void Cluster::RemoveElement(int index)
{
	grids_.erase(grids_.begin() + index);
}
