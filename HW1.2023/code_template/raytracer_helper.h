#ifndef __HELPER__
#define __HELPER__

#include "parser.h"
#include "ray.h"

using namespace ray;

parser::Vec3f FindIntersectionPoint(const Ray &ray, float t);
bool HitSphere(const Ray &ray, parser::Vec3f &center, float radius, int materialID, int objectID, Hit& Hit);
bool HitTriangle(const Ray &ray, parser::Vec3f& a, parser::Vec3f& b, parser::Vec3f& c, int obj_id, Hit& hit);


#endif