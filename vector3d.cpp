#include "vector3d.h"

void Vector3d::set(double _x, double _y, double _z)
{
    x = _x;
    y = _y;
    z = _z;
}

Vector3d Vector3d::operator-(const Vector3d &b)
{
    return Vector3d(x - b.x, y - b.y, z - b.z);
}

Vector3d Vector3d::operator+(const Vector3d &b)
{
    return Vector3d(x + b.x, y + b.y, z + b.z);
}



