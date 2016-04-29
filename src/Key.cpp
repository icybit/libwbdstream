#include "Common.h"
#include "Key.h"

bool Key::operator==(const Key & k) const
{
	return (abs(this->x_ - k.x_) < EPSILON) && (abs(this->y_ - k.y_) < EPSILON);
}

bool Key::operator<(const Key & k) const
{
	return (this->x_ + EPSILON < k.x_) || ((abs(this->x_ - k.x_) < EPSILON) && (this->y_ + EPSILON < k.y_));
}

bool Key::operator>(const Key & k) const
{
	return (this->x_ - EPSILON > k.x_) || ((abs(this->x_ - k.x_) < EPSILON) && (this->y_ - EPSILON > k.y_));
}