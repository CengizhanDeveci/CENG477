#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <limits>
#include <cmath>

using namespace parser;
using namespace std;

#ifndef RAY_H
#define RAY_H

class Ray{
    public:
        Vec3f origin;
        Vec3f direction;
};

class Hit{
        public:
        bool hit_happened;
        Vec3f intersection_point;
        Vec3f normal;
        int material_id;
        float t;
        int type; //spheres are 0, triangles 1 and meshes 2
        int object_id;
};

Vec3f subtract(const Vec3f &obj1, const Vec3f &obj2);

Vec3f add(const Vec3f &obj1, const Vec3f &obj2);

Vec3f multiply(const Vec3f &obj1, float x);

float dot(const Vec3f &obj1, const Vec3f &obj2);

Vec3f cross(const Vec3f &obj1, const Vec3f &obj2);

float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2);

float findLength(const Vec3f &obj);

float findDistance(const Vec3f &obj1, const Vec3f &obj2);

Vec3f normalize(const Vec3f &obj);

Vec3f findIntersectionPoint(const Ray &ray, float t);

Ray generateRay(const Camera &camera, int i, int j);

Vec3f findIrradiance(const PointLight &point_light, const Vec3f &intersection_point);
#endif