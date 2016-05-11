#include "CharacteristicVector.h"
#include "wbdstream/wbcharvectapi.h"

#include <cmath>
#include <cstring>
#include <stdint.h>
#include <stdio.h>

#include "Common.h"

#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)

CharacteristicVector::CharacteristicVector()
{
	this->time_updated_ = 0;
	this->density_ = 0.0f;
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->changed_ = true;
}

CharacteristicVector::CharacteristicVector(uint64_t time_now)
{
	this->time_updated_ = time_now;
	this->density_ = 1.0f; /* Since it's created by hit, add one */
	this->label_ = NO_CLASS;
	this->status_ = INITIAL;
	this->changed_ = true;
}

CharacteristicVector::CharacteristicVector(uint8_t * buffer)
{
	unsigned int index = 0;
	std::memcpy(&this->time_updated_, &buffer[index], sizeof(this->time_updated_));
	index += sizeof(this->time_updated_);
	std::memcpy(&this->density_, &buffer[index], sizeof(this->density_));
	index += sizeof(this->density_);
	std::memcpy(&this->label_, &buffer[index], sizeof(this->label_));
	index += sizeof(this->label_);
	std::memcpy(&this->status_, &buffer[index], sizeof(this->status_));
	index += sizeof(this->status_);
	std::memcpy(&this->changed_, &buffer[index], sizeof(this->changed_));
	index += sizeof(this->changed_);
}

void CharacteristicVector::AddRecord(uint64_t time_now)
{
	this->UpdateDensity(time_now);
	this->density_++;
}

void CharacteristicVector::Merge(uint8_t * buffer)
{
	int index = 0;
	memcpy(&this->label_, &buffer[index], sizeof(this->label_));
	index += sizeof(this->label_);
	memcpy(&this->status_, &buffer[index], sizeof(this->status_));
	this->set_unchanged();
}

void CharacteristicVector::Print()
{
	printf("Time updated: %llu , density: %4.4f, label: %u, status: %d, changed: %d.\n", (unsigned long long)this->time_updated_, this->density_, this->label_, this->status_, this->changed_);
}

void CharacteristicVector::Serialize(uint8_t * buffer)
{		
	int length = 0;
		
	std::memcpy(&buffer[length], &this->time_updated_, sizeof(this->time_updated_));
	length += sizeof(this->time_updated_);

	std::memcpy(&buffer[length], &this->density_, sizeof(this->density_));
	length += sizeof(this->density_);

	std::memcpy(&buffer[length], &this->label_, sizeof(this->label_));
	length += sizeof(this->label_);

	std::memcpy(&buffer[length], &this->status_, sizeof(this->status_));
	length += sizeof(this->status_);
	
	std::memcpy(&buffer[length], &this->changed_, sizeof(this->changed_));
	length += sizeof(this->changed_);
}

void CharacteristicVector::SetStatus(uint8_t status)
{
	if (this->status_ == INITIAL)
	{
		this->set_changed();
		this->status_ = status;
	}
	else if (std::abs(this->status_ - status) > 1)
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

void CharacteristicVector::UpdateDensity(uint64_t time_now)
{
	this->density_ = (float)pow(Common::decay_factor, (float)(time_now - this->time_updated_)) * this->density_; /* (5) - research paper */
	this->time_updated_ = time_now;
}

DSTREAM_PUBLIC dstream_char_vect_t * dstream_char_vect_new(uint64_t time_now)
{
	return AS_TYPE(dstream_char_vect_s, new CharacteristicVector(time_now));
}

DSTREAM_PUBLIC void dstream_char_vect_free(void ** vect)
{
	if (!*vect)
		return;
	delete AS_TYPE(CharacteristicVector, *vect);
}

DSTREAM_PUBLIC void dstream_char_vect_add_record(dstream_char_vect_t * vect, uint64_t time_now)
{
	AS_TYPE(CharacteristicVector, vect)->AddRecord(time_now);
}

DSTREAM_PUBLIC void dstream_char_vect_merge(dstream_char_vect_t * vect, uint8_t * buffer)
{
	AS_TYPE(CharacteristicVector, vect)->Merge(buffer);
}

DSTREAM_PUBLIC void dstream_char_vect_print(dstream_char_vect_t * vect)
{
	AS_TYPE(CharacteristicVector, vect)->Print();
}

DSTREAM_PUBLIC void dstream_char_vect_serialize(dstream_char_vect_t * vect, uint8_t * buffer)
{
	AS_TYPE(CharacteristicVector, vect)->Serialize(buffer);
}

DSTREAM_PUBLIC void dstream_char_vect_update_density(dstream_char_vect_t * vect, uint64_t time_now)
{
	AS_TYPE(CharacteristicVector, vect)->UpdateDensity(time_now);
}