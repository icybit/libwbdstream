#include "CharacteristicVector.h"

#include <math.h>


CharacteristicVector::CharacteristicVector()
{
	this->time_created_ = 0;
	this->time_updated_ = 0;
	this->label_ = NO_CLASS;
	this->status_ = false; // non-sporadic
	this->density_ = 0.0f;
}

CharacteristicVector::CharacteristicVector(unsigned __int64 time_created) 
{
	this->time_created_ = time_created;
	this->time_updated_ = time_created;
	this->status_ = false; // non-sporadic
	this->label_ = NO_CLASS;
	this->density_ = 1.0f;
}

void CharacteristicVector::UpdateDensity(unsigned __int64 time_now)
{
	this->density_ = (float)pow(DECAY_FACTOR, (double)(time_now - this->time_updated_)) * density_; // (5) - research paper
	this->time_updated_ = time_now;
}

void CharacteristicVector::AddRecord(unsigned __int64 time_now)
{
	this->UpdateDensity(time_now);
	this->density_++;
}
