#include "Cluster.h"


Cluster::Cluster(unsigned int label)
{
	this->label_ = label;
}

void Cluster::GetPair(int index, float & x, float & y)
{
	x = grids_[index].x;
	y = grids_[index].y;
}

void Cluster::MergeClusters(Cluster * cluster)
{
	for (int i = 0; i < cluster->get_size(); i++)
	{
		this->AddElement(cluster->PopBack());
	}
}

GridTuple Cluster::PopBack()
{
	GridTuple last = *(grids_.end() - 1);
	grids_.pop_back();
	return last;
}

void Cluster::RemoveElement(int index)
{
	grids_.erase(grids_.begin() + index);
}

void Cluster::RemoveElement(float x, float y)
{
	int i;
	bool found = false;
	for (i = 0; i < grids_.size(); i++)
	{
		if (abs(grids_[i].x - x) < EPSILON && abs(grids_[i].y - y) < EPSILON)
		{
			found = true;
			break;
		}
	}
	if(found)
		RemoveElement(i);
}
