#include "bvh.h"

bool BoundingBox::Intersect(const Ray& ray, Hit& hit) const
{

}

BVHNode* BoundingVolume::BuildBVH(std::vector<Object>& objects, int start, int end)
{
    if(start == end)
    {
        BVHNode* leaf = new BVHNode();
        leaf->boundingBox = BoundingBox::CalculateBoundingBox(objects[start]);
        leaf->object = &objects[start];
        return leaf;
    }

    int splitAxis = ChooseSplitAxis();

}

int ChooseSplitAxis() // TODO choose split axis
{
    return 0;
}

bool CompareOnAxis(int splitAxis, const Object& obj1, const Object& obj2)
{
    
}