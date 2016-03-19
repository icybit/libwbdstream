#include "CharacteristicVector.h"

#include <math.h>


CharacteristicVector::CharacteristicVector()
{
	this->time_updated_ = 0;
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->density_ = 0.0f;
	this->changed_ = true;
	this->neighbors_ = 0;
}

CharacteristicVector::CharacteristicVector(unsigned __int64 time_updated)
{
	this->time_updated_ = time_updated;
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->density_ = 1.0f; // Since it's created by hit, add one
	this->changed_ = true;
	this->neighbors_ = 0;
}

void CharacteristicVector::AddNeighbors(int count)
{
	this->neighbors_ += count;
}

void CharacteristicVector::AddRecord(unsigned __int64 time_now)
{
	this->UpdateDensity(time_now);
	this->density_++;
}

void CharacteristicVector::SetStatus(unsigned __int8 status)
{
	if (this->status_ == INITIAL)
	{
		this->set_changed();
		this->status_ = status;
	}
	else if (abs(this->status_ - status) > 1) // It has changed - look for numbering in header file
	{
		this->set_changed();
		if (status == TRANSITIONAL)
		{
			if (status > this->status_)
			{
				this->status_ = TRANSITIONAL_FROM_SPARSE;
			}
			else
			{
				this->status_ = TRANSITIONAL_FROM_DENSE;
			}
		}
		else if (status == SPARSE)
		{
			if (status + 4 < this->status_)
			{
				this->status_ = SPARSE_FROM_DENSE;
			}
			else
			{
				this->status_ = SPARSE_FROM_TRANSITIONAL;
			}
		}
		else if (status == DENSE)
		{
			if (status  - 4 > this->status_)
			{
				this->status_ = DENSE_FROM_SPARSE;
			}
			else
			{
				this->status_ = DENSE_FROM_TRANSITIONAL;
			}
		}
		else if (status == SPORADIC)
		{
			this->set_unchanged();
			this->status_ = status;
		}
	}
	else
	{
		this->status_ = status;
	}
}

void CharacteristicVector::UpdateDensity(unsigned __int64 time_now)
{
	this->density_ = (float)pow(DECAY_FACTOR, (double)(time_now - this->time_updated_)) * this->density_; // (5) - research paper
	this->time_updated_ = time_now;
}

