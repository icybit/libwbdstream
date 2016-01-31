#pragma once

#define NO_CLASS 0
#define DECAY_FACTOR 0.998

class CharacteristicVector
{
private:
	unsigned __int64 time_updated_;
	unsigned __int64 time_created_;
	float density_;
	unsigned int label_; //cluster label
	bool status_; //is sporadic?	
public:
	CharacteristicVector();

	CharacteristicVector(unsigned __int64 time_created);

	void UpdateDensity(unsigned __int64 time_now);

	void AddRecord(unsigned __int64 time_now);

	// Setters

	void set_label(unsigned int label) {
		this->label_ = label;
	}

	void set_status(bool sporadic) { 
		this->status_ = sporadic; 
	}
	
	// Getters

	unsigned __int64 get_time_updated() {
		return time_updated_;
	}

	unsigned __int64 get_time_removed() {
		return time_created_;
	}
	
	float get_density() { return density_; }
	
	//Returns label of the grid (ID of the cluster)
	unsigned int get_label() { return label_; }

	//Returns status, is grid sporadic 
	bool is_sporadic() { return status_; }

	

};

