#include "CharacteristicVector.h"

#include <math.h>


CharacteristicVector::CharacteristicVector()
{
	this->time_updated_ = 0;
	this->density_ = 0.0f;
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->changed_ = true;
}

CharacteristicVector::CharacteristicVector(unsigned __int64 time_updated)
{
	this->time_updated_ = time_updated;
	this->density_ = 1.0f; // Since it's created by hit, add one
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->changed_ = true;
}

CharacteristicVector::CharacteristicVector(unsigned char * buffer)
{
	unsigned int index = 0;
	memcpy(&this->time_updated_, &buffer[index], sizeof(this->time_updated_));
	index += sizeof(this->time_updated_);
	memcpy(&this->density_, &buffer[index], sizeof(this->density_));
	index += sizeof(this->density_);
	memcpy(&this->label_, &buffer[index], sizeof(this->label_));
	index += sizeof(this->label_);
	memcpy(&this->status_, &buffer[index], sizeof(this->status_));
	index += sizeof(this->status_);
	memcpy(&this->changed_, &buffer[index], sizeof(this->changed_));
	index += sizeof(this->changed_);
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

