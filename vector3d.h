#ifndef VECTOR3D_H
#define VECTOR3D_H


class Vector3d
{
public:
    double x, y, z;
    void set(double _x, double _y, double _z);
    Vector3d(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
    Vector3d(){x = y = z = 0;}
    Vector3d operator-(const Vector3d& b);
    Vector3d operator+(const Vector3d& b);
};


#endif // VECTOR3D_H
