#ifndef PARTICLE_H
#define PARTICLE_H

#include "vector3d.h"

class Particle
{
private:
    Vector3d positionLC; // in lattice constant
    int type;
public:
    int getType(){return type;}
    void setType(int _type){type = _type;}
    Particle(double x, double y, double z, int type);
    Vector3d getPosition(double, int, double, int a = 0, int b = 0, int c = 0);
    Vector3d getPositionLC();
    void print();
};

#endif // PARTICLE_H
