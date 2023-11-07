#include <algorithm>
#include <iostream>
#include "parser.h"
#include "ray.h"
#include "bvh.h"
#include <limits>
#include "raytracer_helper.h"

#define epsilon 10e-3

using namespace parser;

BVH::BVH(std::vector<Object*> &objects)
{
    root = BuildBVH(objects, 0, objects.size() - 1, 0);
}

BVHNode *BVH::BuildBVH(std::vector<Object*> &objects, int start, int end, int depth)
{
    if(start == end)
    {
        BVHNode* leaf = new BVHNode();
        leaf->boundingBox = objects[start]->boundingBox;
        leaf->object = objects[start];
        return leaf;
    }

    int splitAxis = depth % 3;
    int mid = (start + end) / 2;

    std::sort(objects.begin() + start, objects.begin() + end, [splitAxis](const Object* a, const Object* b) -> bool
        {
            if(splitAxis == 0)
            {
                return a->boundingBox->min.x < b->boundingBox->min.x;
            }
            else if(splitAxis == 1)
            {
                return a->boundingBox->min.y < b->boundingBox->min.y;
            }
            else
            {
                return a->boundingBox->min.z < b->boundingBox->min.z;
            }
        });

    BVHNode* left = BuildBVH(objects, start, mid, depth + 1);
    BVHNode* right = BuildBVH(objects, mid + 1, end, depth + 1);

    BVHNode* parent = new BVHNode();
    parent->left = left;
    parent->right = right;

    parent->boundingBox = new BoundingBox();
    parent->boundingBox->min.x = std::min(left->boundingBox->min.x, right->boundingBox->min.x);
    parent->boundingBox->min.y = std::min(left->boundingBox->min.y, right->boundingBox->min.y);
    parent->boundingBox->min.z = std::min(left->boundingBox->min.z, right->boundingBox->min.z);
    parent->boundingBox->max.x = std::max(left->boundingBox->max.x, right->boundingBox->max.x);
    parent->boundingBox->max.y = std::max(left->boundingBox->max.y, right->boundingBox->max.y);
    parent->boundingBox->max.z = std::max(left->boundingBox->max.z, right->boundingBox->max.z);

    return parent;
}



bool BVH::Intersect(Ray &ray, BVHNode* node, Hit& hit)
{
    hit.hit = false;
    if(node == nullptr)
    {
        return false;
    }

    Hit tmpHit;
    if(node->boundingBox->Intersect(ray, tmpHit))
    {
        if(node->left == nullptr && node->right == nullptr && node->object != nullptr)
        {
            
            if(node->object->isSphere)
            {
                hit.hit = false;
                bool tmp = HitSphere(ray, node->object->center, node->object->radius, node->object->material_id, node->object->object_id, hit);
                if(hit.hit) std::cout << hit.t << std::endl;
            }
            else 
            {
                bool checkHit = node->object->boundingBox->Intersect(ray, hit);
                if(checkHit)
                {
                    hit.hit = true;
                    hit.intersectionPoint = AddVec3f(ray.origin, MultiplyVec3f(ray.direction, hit.t));
                    hit.materialID = node->object->material_id;
                    hit.normal = node->object->normal; // TODO calculate normal
                    hit.objectID = node->object->object_id;
                    hit.type = node->object->type;
                }
                else
                {
                    hit.hit = false;
                }
            }
            return hit.hit;
            
            return node->object->boundingBox->Intersect(ray, hit);
        }

        Hit hit1, hit2;
        hit.t = std::numeric_limits<float>::max();

        bool leftHitCheck = Intersect(ray, node->left, hit1);
        bool rightHitCheck = Intersect(ray, node->right, hit2);

        if(leftHitCheck)
        {
            hit = hit1;
        }

        if(rightHitCheck && hit2.t < hit.t)
        {
            hit = hit2;
        }

        return leftHitCheck || rightHitCheck;

        // if(node->left == nullptr && node->right == nullptr && node->object != nullptr)
        // {
        //     return node->object->boundingBox->Intersect(ray, hit);
        // }
        // else
        // {
        //     bool leftHitCheck = Intersect(ray, node->left, hit);
        //     bool rightHitCheck = Intersect(ray, node->right, hit);
        //     return leftHitCheck || rightHitCheck;
        // }
    }
    else
    {
        return false;
    }
}

bool BoundingBox::Intersect(Ray& ray, Hit& hit)
{
    float tmin = (min.x - ray.origin.x) / ray.direction.x;
    float tmax = (max.x - ray.origin.x) / ray.direction.x;

    float tymin = (min.y - ray.origin.y) / ray.direction.y;
    float tymax = (max.y - ray.origin.y) / ray.direction.y;

    if(tmin <= epsilon && tmax <= epsilon)
    {
        return false;
    }

    if(tymin <= epsilon && tymax <= epsilon)
    {
        return false;
    }

    if(tmin > tmax)
    {
        std::swap(tmin, tmax);
    }

    if(tymin > tymax)
    {
        std::swap(tymin, tymax);
    }

    if((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }

    if(tymin > tmin)
    {
        tmin = tymin;
    }

    if(tymax < tmax)
    {
        tmax = tymax;
    }

    float tzmin = (min.z - ray.origin.z) / ray.direction.z;
    float tzmax = (max.z - ray.origin.z) / ray.direction.z;

    if(tzmin <= epsilon && tzmax <= epsilon)
    {
        return false;
    }

    if(tzmin > tzmax)
    {
        std::swap(tzmin, tzmax);
    }

    if((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }
    
    if(tzmin > tmin)
    {
        tmin = tzmin;
    }

    if(tzmax < tmax)
    {
        tmax = tzmax;
    }

    
    hit.t = tmin;
    hit.hit = true;

    return true;	
}


