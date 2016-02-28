#ifndef STM_D_STREAM_CHARACTERISTIC_VECTOR_H_
#define STM_D_STREAM_CHARACTERISTIC_VECTOR_H_

#define NO_CLASS 0

#define SPORADIC 0

#define SPARSE_FROM_TRANSITIONAL 1
#define SPARSE 2
#define SPARSE_FROM_DENSE 3

#define TRANSITIONAL_FROM_SPARSE 4
#define TRANSITIONAL 5
#define TRANSITIONAL_FROM_DENSE 6

#define DENSE_FROM_SPARSE 7
#define DENSE 8
#define DENSE_FROM_TRANSITIONAL 9

#define INITIAL 10

#define DECAY_FACTOR 0.998f

class CharacteristicVector
{	
public:
	CharacteristicVector();

	CharacteristicVector(unsigned __int64 time_updated, float density, unsigned int label,
		unsigned __int8 status, bool changed);

	void AddRecord(unsigned __int64 time_now);

	void SetStatus(unsigned __int8 status);

	void UpdateDensity(unsigned __int64 time_now);

	// Setters

	void set_label(unsigned int label) {
		this->label_ = label;
	}

	void set_changed() { this->changed_ = true; }

	void set_unchanged() { this->changed_ = false; }
	
	// Getters

	unsigned __int64 get_time_updated() {
		return time_updated_;
	}
	
	float get_density() { return density_; }
	
	//Returns label of the grid (ID of the cluster)
	unsigned int get_label() { return label_; }

	//Returns status, is grid sporadic 
	unsigned __int8 get_status() { return status_; }

	//Returns, if status has changed
	bool is_changed() {
		return changed_;
	}

private:
	unsigned __int64 time_updated_;
	float density_;
	unsigned int label_; //cluster label
	unsigned __int8 status_; // sporadic/sparse/transitional/dense
	bool changed_; 
};

#endif // !STM_D_STREAM_CHARACTERISTIC_VECTOR_H_