#ifndef _DSTREAM_CHARACTERISTIC_VECTOR_H_
#define _DSTREAM_CHARACTERISTIC_VECTOR_H_

#include <stdint.h>
#include <tuple>

#include "Common.h"

class CharacteristicVector
{	
public:
	CharacteristicVector();

	CharacteristicVector(uint64_t time_now);

	CharacteristicVector(uint8_t * buffer);

	void AddRecord(uint64_t time_now);
	
	void Serialize(uint8_t * buffer);

	void SetStatus(uint8_t status);

	void UpdateDensity(uint64_t time_now);

	// Setters

	void set_label(unsigned int label) {
		this->label_ = label;
	}

	void set_changed() { this->changed_ = true; }

	void set_unchanged() { this->changed_ = false; }
	
	// Getters

	uint64_t get_time_updated() {
		return time_updated_;
	}
	
	//Returns density
	float get_density() { return this->density_; }
	
	//Returns label of the grid (ID of the cluster)
	uint32_t get_label() { return label_; }

	//Returns status, SPORADIC, SPARSE, TRANSITIONAL etc; 
	uint8_t get_status() { return status_; }

	//Returns, if status has changed
	bool is_changed() {
		return changed_;
	}
private:
	uint64_t time_updated_;
	float density_;
	uint32_t label_;
	uint8_t status_; 
	bool changed_; 
};

#endif // !_DSTREAM_CHARACTERISTIC_VECTOR_H_