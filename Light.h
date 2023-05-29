#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "LSystems.h"

class Light {
public:
    Color ambientLight;
    Color diffuseLight;
    Color specularLight;
};

class InfLight: public Light{
public:
    Vector3D ldVector;
};

class Pointlight: public Light{
public:
    Vector3D location;
    double spotAngle;
};

typedef std::list<Light> Lights3D;
#endif //ENGINE_LIGHT_H
