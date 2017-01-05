#include "particle.h"
#include <iostream>
#include <cmath>

Particle::Particle(double x, double y, double z, int _type) {
    type = _type;
    positionLC.set(x, y, z);
}

Vector3d Particle::getPosition(double a, int typeTransform, double alpha, int dx, int dy, int dz)
{
    if (typeTransform == 0){
        return Vector3d((positionLC.x + dx) * a * (1.0 + alpha),
                        (positionLC.y + dy) * a * (1.0 + alpha),
                        (positionLC.z + dz) * a * (1.0 + alpha));
    }
    else if (typeTransform == 11){
         return Vector3d((positionLC.x + dx) * a * (1.0 + alpha),
                         (positionLC.y + dy) * a * (1.0 + alpha),
                         (positionLC.z + dz) * a);
    }
    else if(typeTransform == 12){
        return Vector3d((positionLC.x + dx) * a * (1.0 + alpha),
                        (positionLC.y + dy) * a * (1.0 - alpha),
                        (positionLC.z + dz) * a);
    }
    else if(typeTransform == 44){
        return Vector3d((positionLC.x + dx) * a + (positionLC.y + dy) * a * alpha,
                        (positionLC.x + dx) * a * alpha + (positionLC.y + dy) * a,
                        (positionLC.z + dz) * a / (1.0 - pow(alpha, 2)));
    }
}

Vector3d Particle::getPositionLC()
{
    return positionLC;
}

void Particle::print()
{
    std::cerr << positionLC.x << " " << positionLC.y << " " <<
                 positionLC.z << std::endl;
}


