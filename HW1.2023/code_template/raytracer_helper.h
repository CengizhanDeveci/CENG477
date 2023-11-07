#ifndef __HELPER__
#define __HELPER__

#include "parser.h"
#include "ray.h"

using namespace ray;

parser::Vec3f FindIntersectionPoint(const Ray &ray, float t);
bool HitSphere(const Ray &ray, parser::Vec3f &center, float radius, int materialID, int objectID, Hit& Hit);


#endif