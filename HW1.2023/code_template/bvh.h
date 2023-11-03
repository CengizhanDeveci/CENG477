#ifndef __BVH__
#define __BFH__

#include "ray.h"
#include "parser.h"

using namespace ray;
using namespace parser;

struct Object
{

};

struct BVHNode
{
    BoundingBox boundingBox;
    BVHNode* left;
    BVHNode* right;
    Object* object;
};

class BoundingBox
{
    public:
        // variables
        parser::Vec3f minCorner;
        parser::Vec3f maxCorner;
        parser::Vec3f center;

        // methods
        bool Intersect(const Ray& ray, Hit& hit) const;

        // static methods
        static BoundingBox CalculateBoundingBox(const Object& object);

};

class BoundingVolume
{
    public:
        // methods
        BVHNode* BuildBVH(std::vector<Object>& objects, int start, int end);
        

        // variables
        BVHNode* parent;

};

int ChooseSplitAxis();
bool CompareOnAxis(int splitAxis, const Object& obj1, const Object& obj2);
#endif