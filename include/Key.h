#ifndef STM_D_STREAM_KEY_H_
#define STM_D_STREAM_KEY_H_

#define EPSILON 0.00001f

#include <functional>  

struct Key
{
	float x_;
	float y_;
	Key() : x_(-1.0f), y_(-1.0f) {}
	Key(float x, float y) : x_(x), y_(y) {}

	bool empty()
	{
		return (abs(this->x_ + 1.0f) < EPSILON) && (abs(this->y_ + 1.0f) < EPSILON);
	}

	bool operator==(const Key & k) const
	{
		return (abs(this->x_ - k.x_) < EPSILON) && (abs(this->y_ - k.y_) < EPSILON);
	}

	bool operator<(const Key & k) const
	{
		return (this->x_ + EPSILON < k.x_) || ((abs(this->x_ - k.x_) < EPSILON) && (this->y_ + EPSILON < k.y_));
	}

	bool operator>(const Key & k) const
	{
		return (this->x_ - EPSILON > k.x_) || ((abs(this->x_ - k.x_) < EPSILON) && (this->y_ - EPSILON > k.y_));
	}
};

namespace std
{
	template<>
	struct hash<Key>
	{
		size_t operator()(const Key & key) const
		{
			size_t const h1(std::hash<float>()(key.x_));
			size_t const h2(std::hash<float>()(key.y_));
			return h1 ^ (h2 << 1);
		}
	};
}

#endif // !STM_D_STREAM_KEY_H_
