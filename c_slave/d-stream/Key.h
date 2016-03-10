#ifndef STM_D_STREAM_KEY_H_
#define STM_D_STREAM_KEY_H_

#include <functional>  

struct Key
{
	float x_;
	float y_;
	Key() {}
	Key(float x, float y) : x_(x), y_(y) {}

	bool operator==(const Key & k) const
	{
		return this->x_ == k.x_ && this->y_ == k.y_;
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
