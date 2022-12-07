#include "ray.h"



Vec3f subtract(const Vec3f &obj1, const Vec3f &obj2){
    Vec3f result;
    result.x = obj1.x - obj2.x;
    result.y = obj1.y - obj2.y;
    result.z = obj1.z - obj2.z;
    return result;
}

Vec3f add(const Vec3f &obj1, const Vec3f &obj2){
    Vec3f result;
    result.x = obj1.x + obj2.x;
    result.y = obj1.y + obj2.y;
    result.z = obj1.z + obj2.z;
    return result;
}

Vec3f multiply(const Vec3f &obj1, float number){
    Vec3f result;
    result.x = obj1.x * number;
    result.y = obj1.y * number;
    result.z = obj1.z * number;
    return result;
}

float dot(const Vec3f &obj1, const Vec3f &obj2){
    float result = 0;
    result += obj1.x * obj2.x + obj1.y * obj2.y + obj1.z * obj2.z;
    return result;
}

Vec3f cross(const Vec3f &obj1, const Vec3f &obj2){
    Vec3f result;
    result.x = obj1.y * obj2.z - obj1.z * obj2.y;
    result.y = obj1.z * obj2.x - obj1.x * obj2.z;
    result.z = obj1.x * obj2.y - obj1.y * obj2.x;
    return result;
}

float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2){
	return v0.x * (v1.y*v2.z - v2.y*v1.z) + v0.y * (v2.x*v1.z - v1.x*v2.z) + v0.z * (v1.x*v2.y - v1.y*v2.x);
}

float findLength(const Vec3f &obj){
    float result = 0;
    result = sqrt(obj.x * obj.x + obj.y * obj.y + obj.z * obj.z);
    return result;
}

float findDistance(const Vec3f &obj1, const Vec3f &obj2){
    return sqrt(pow(obj1.x - obj2.x, 2) + pow(obj1.y - obj2.y, 2) + pow(obj1.z - obj2.z, 2));
}

Vec3f normalize(const Vec3f &obj){
    Vec3f result;
    result.x = obj.x / findLength(obj);
    result.y = obj.y / findLength(obj);
    result.z = obj.z / findLength(obj);
    return result;
}

Vec3f findIntersectionPoint(const Ray &ray, float t){
	Vec3f result;
	result.x = ray.origin.x + t*ray.direction.x;
	result.y = ray.origin.y + t*ray.direction.y;
	result.z = ray.origin.z + t*ray.direction.z;
	return result;
}

Ray generateRay(const Camera &camera, int i, int j){
    Ray result;
    float su, sv;
    Vec3f m, q, s; 
    Vec3f gaze = camera.gaze;
    gaze = normalize(gaze);

    Vec3f position = camera.position;
    float top = camera.near_plane.w;
    float bottom = camera.near_plane.z;
    float right = camera.near_plane.y;
    float left = camera.near_plane.x;

    Vec3f u, w, v;

    su = (i + 0.5) * (right - left) / camera.image_width;
    sv = (j + 0.5) * (top - bottom) / camera.image_height;
    
    m = add(position, multiply(gaze, camera.near_distance));
    u = cross(gaze, camera.up);
    u = normalize(u);

    v = cross(u, gaze);

    q = add(m, add(multiply(u, left), multiply(v, top)));
    s = add(q, subtract(multiply(u, su), multiply(v, sv)));

    result.origin = position;
    result.direction = subtract(s, position);
    result.direction = normalize(result.direction);

    return result;
}

Vec3f findIrradiance(const PointLight &point_light, const Vec3f &intersection_point){
    Vec3f irradiance;
    Vec3f intensity;
    Vec3f d;
    intensity = point_light.intensity;
    d = subtract(point_light.position, intersection_point);
    float d_dot_d = dot(d, d);

    if(d_dot_d == 0.0) return irradiance;
    
    irradiance.x = intensity.x/d_dot_d;
    irradiance.y = intensity.y/d_dot_d;
    irradiance.z = intensity.z/d_dot_d;
    
    return irradiance;
}