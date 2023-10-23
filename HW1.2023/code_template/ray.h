#ifndef __RAY_H__
#define __RAY_H__

#include "parser.h"

namespace ray
{
    struct Ray
    {
        parser::Vec3f origin;
        parser::Vec3f direction;
    };

    struct Hit
    {
        bool hit;
        parser::Vec3f intersectionPoint;
        parser::Vec3f normal;
        int materialID;
        int type; //spheres are 0, triangles 1 and meshes 2
        int objectID;
        float t;
    };

    Ray GenerateRay(const parser::Camera &camera, int i, int j);
    void NormalizeVec3f(parser::Vec3f& vec);
    parser::Vec3f CrossProduct(const parser::Vec3f &v1, const parser::Vec3f &v2);
    inline float DotProduct(const parser::Vec3f &v1, const parser::Vec3f &v2) {return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;};
    parser::Vec3f AddVec3f(const parser::Vec3f &v1, const parser::Vec3f &v2);
    parser::Vec3f SubtractVec3f(const parser::Vec3f &v1, const parser::Vec3f &v2);
    parser::Vec3f MultiplyVec3f(const parser::Vec3f &v, const float& value);
    parser::Vec3f NegateVec3f(const parser::Vec3f &v);
    float Determinant(const parser::Vec3f &v0, const parser::Vec3f &v1, const parser::Vec3f &v2);
}
#endif