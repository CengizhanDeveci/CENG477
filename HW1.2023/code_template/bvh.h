#ifndef __BVH__
#define __BVH__

#include "ray.h"
#include "parser.h"
#include <vector>

using namespace ray;
using namespace parser;

class BoundingBox
{
    public:
        parser::Vec3f min;
        parser::Vec3f max;

        BoundingBox()
        {
            min = parser::Vec3f();
            max = parser::Vec3f();
        }

        bool Intersect(ray::Ray& ray, ray::Hit& hit);

};



struct BVHNode
{
    BVHNode* left;
    BVHNode* right;
    BoundingBox* boundingBox;
    Object* object;
};



class BVH
{
    public:
        BVHNode* root;
        BVHNode* BuildBVH(std::vector<Object*>& objects, int start, int end, int depth);
        BVH(std::vector<Object*>& objects);
        bool Intersect(ray::Ray& ray, BVHNode* node, ray::Hit& hit);

};

#endif