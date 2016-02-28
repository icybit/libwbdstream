#include "CharacteristicVector.h"

#include <math.h>


CharacteristicVector::CharacteristicVector()
{
	this->time_updated_ = 0;
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->density_ = 0.0f;
	this->changed_ = true;
}

CharacteristicVector::CharacteristicVector(unsigned __int64 time_updated, float density, unsigned int label,
	unsigned __int8 status, bool changed)
{
	this->time_updated_ = time_updated;
	this->status_ = status;
	this->label_ = label;
	this->density_ = density;
	this->changed_ = changed;
}

void CharacteristicVector::AddRecord(unsigned __int64 time_now)
{
	this->UpdateDensity(time_now);
	this->density_++;
}

void CharacteristicVector::SetStatus(unsigned __int8 status)
{
	if (abs(this->status_ - status) > 1 && this->status_ != INITIAL) // It has changed - look for numbering in header file
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
	}
	else if (this->status_ == INITIAL)
	{
		this->status_ = status;
	}
}

void CharacteristicVector::UpdateDensity(unsigned __int64 time_now)
{
	this->density_ = (float)pow(DECAY_FACTOR, (double)(time_now - this->time_updated_)) * density_; // (5) - research paper
	this->time_updated_ = time_now;
}

