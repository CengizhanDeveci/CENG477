#include "ray.h"
#include <cmath>

using namespace ray;

Ray ray::GenerateRay(const parser::Camera & camera, int i, int j)
{
    Ray result;
    parser::Vec3f camPos = camera.position;

    parser::Vec3f u, v, w;
    parser::Vec3f m, q, s;

    float left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float bottom = camera.near_plane.z;
    float top = camera.near_plane.w;

    w = ray::NegateVec3f(camera.gaze);
    v = camera.up;
    u = ray::CrossProduct(v, w);
    u = ray::NormalizeVec3f(u);
    v = ray::NormalizeVec3f(v);
    w = ray::NormalizeVec3f(w);

    float su = (i + 0.5f) * ((right - left) / camera.image_width);
    float sv = (j + 0.5f) * ((top - bottom) / camera.image_height);


    m = ray::SubtractVec3f(camPos, ray::MultiplyVec3f(w, camera.near_distance)); 
    q = ray::AddVec3f(m, ray::AddVec3f(ray::MultiplyVec3f(u, left), ray::MultiplyVec3f(v, top)));
    s = ray::AddVec3f(q, ray::SubtractVec3f(ray::MultiplyVec3f(u, su), ray::MultiplyVec3f(v, sv)));

    result.origin = camPos;
    result.direction = ray::SubtractVec3f(s, camPos);
    result.direction = ray::NormalizeVec3f(result.direction);
    
    return result;
}

parser::Vec3f ray::NormalizeVec3f(parser::Vec3f& vec)
{
    parser::Vec3f result;
    float sumSqrt = sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    result.x = vec.x / sumSqrt;
    result.y = vec.y / sumSqrt;
    result.z = vec.z / sumSqrt;

    return result;
}

parser::Vec3f ray::CrossProduct(const parser::Vec3f &v1, const parser::Vec3f &v2)
{
    parser::Vec3f result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;

    return result;
}

parser::Vec3f ray::AddVec3f(const parser::Vec3f &v1, const parser::Vec3f &v2)
{
    parser::Vec3f result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;

    return result;
}

parser::Vec3f ray::SubtractVec3f(const parser::Vec3f &v1, const parser::Vec3f &v2)
{
    parser::Vec3f result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;

    return result;
}

parser::Vec3f ray::MultiplyVec3f(const parser::Vec3f &v, const float& value)
{
    parser::Vec3f result; 
    result.x = v.x * value;
    result.y = v.y * value;
    result.z = v.z * value;

    return result;
}

parser::Vec3f ray::NegateVec3f(const parser::Vec3f &v)
{
    parser::Vec3f result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

float ray::Determinant(const parser::Vec3f &v0, const parser::Vec3f &v1, const parser::Vec3f &v2)
{
	return v0.x * (v1.y*v2.z - v2.y*v1.z) + v0.y * (v2.x*v1.z - v1.x*v2.z) + v0.z * (v1.x*v2.y - v1.y*v2.x);
}