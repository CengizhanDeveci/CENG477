#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <limits>
#include <cmath>

using namespace std;
using namespace parser;

typedef unsigned char RGB[3];

class Ray{
    public:
        Vec3f origin;
        Vec3f direction;
        bool shadowRay;
}

class Hit{
    bool hit_happened;
    Vec3f intersection_point;
    Vec3f normal;
    int material_id;
    float t;
    int type;
    int object_id;
}

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

float findLength(const Vec3f &obj){
    float result = 0;
    result = sqrt(obj.x * obj.x + obj.y * obj.y + obj.z * obj.z);
    return result;
}

float findDistance(const Vec3f &obj1, const Vec3f &obj2){
    return sqrt(pow(obj1.x - obj2.x, 2) + pow(obj1.y - obj2.y, 2) + pow(obj1.z - obj2.z), 2);
}

Vec3f normalize(const Vec3f &obj){
    Vec3f result;
    result.x = obj1.x / findLength(obj);
    result.y = obj1.y / findLength(obj);
    result.z = obj1.z / findLength(obj);
    return result;
}

Vec3f findIntersectionPoint(const Ray &ray, float t){
	Vec3f result;
	result.x = ray.origin.x + t*ray.direction.x;
	result.y = ray.origin.y + t*ray.direction.y;
	result.z = ray.origin.z + t*ray.direction.z;
	return result;
}

Hit findHit(const Scene &scene, const Ray &ray){
    Hit hitted;
    hitted.hit_happened = false;
    float t = numeric_limits<float>::max();

    
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    
    int number_of_cameras = scene.cameras.size();

    for(int camera_no = 0; camera_no < number_of_cameras; camera_no++){
        int width = scene.cameras[camera_no].image_width;
        int height = scene.cameras[camera_no].image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                Ray ray = GenerateRay(scene.cameras[camera_no], j, i);
                Vec3f color;

                color = computeColor(ray, scene, scene.cameras[camera_no], scene.max_recursion_depth);

                // if the color values greater than 255 we set 255 otherwise round the nearest integer value
                if(color.x > 255) image[3* i * width + 3 * j] = 255;
                else image[3* i * width + 3 * j] = (int) (color.x + 0.5);

                if(color.y > 255) image[3 * i * width + 3 * j + 1] = 255;
                else image[3* i * width + 3 * j + 1] = (int) (color.y + 0.5);

                if(color.z > 255) image[3 * i * width + 3 * j + 2] = 255;
                else image[3 * i * width + 3 * j + 2] = (int) (color.z + 0.5);
                
            }
        }
        
        write_ppm(scene.cameras[camera_no].image_name.c_str(), image, width, height);

    }
}

