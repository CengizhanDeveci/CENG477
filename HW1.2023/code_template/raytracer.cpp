#include <iostream>
#include <cmath>
#include "parser.h"
#include "ppm.h"
#include "ray.h"
#include <limits>

using namespace ray;

typedef unsigned char RGB[3];

parser::Vec3f FindIntersectionPoint(const Ray &ray, float t)
{
    parser::Vec3f result;
    result.x = ray.origin.x + t * ray.direction.x;
    result.y = ray.origin.y + t * ray.direction.y;
    result.z = ray.origin.z + t * ray.direction.z;

    return result;
}

Hit HitSphere(const Ray &ray, parser::Vec3f &center, float radius, int materialID, int objectID)
{
    Hit hit;

    float A = DotProduct(ray.direction, ray.direction);
	float B = 2 * DotProduct(ray.direction, SubtractVec3f(ray.origin, center));
	float C = DotProduct(SubtractVec3f(ray.origin, center), SubtractVec3f(ray.origin, center)) - radius * radius;

	float delta = B * B - 4 * A * C;

    if(delta < 0) // no intersection
    {
        hit.hit = false;
    }
    else // intersects
    {
        float t1 = (-1 * B + sqrtf(delta))/2*A;
		float t2 = (-1 * B - sqrtf(delta))/2*A;

        hit.hit = true;

		hit.materialID = materialID;
		hit.type = 0;
		hit.objectID = objectID;

		float t = (t1 > t2) ? t2 : t1;
        hit.t = t;
		hit.intersectionPoint = FindIntersectionPoint(ray, t);

		hit.normal = SubtractVec3f(hit.intersectionPoint, center);
		hit.normal.x /= radius;
		hit.normal.y /= radius;
		hit.normal.z /= radius;
        NormalizeVec3f(hit.normal);
    }

    return hit;
}

Hit HitTriangle(const Ray &ray, const parser::Scene &scene, parser::Triangle &triangle, int obj_id){
	Hit hit;
	hit.hit = false;

    parser::Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
    parser::Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
    parser::Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];

    parser::Vec3f normal;
    normal = CrossProduct(SubtractVec3f(c, b), SubtractVec3f(a, b));
    NormalizeVec3f(normal);

    // some calculations for barycentric coordinates
    float detA = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), ray.direction);
    float beta = Determinant(SubtractVec3f(a, ray.origin), SubtractVec3f(a, c), ray.direction) / detA;
    float gamma = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, ray.origin), ray.direction) / detA;
    float t = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), SubtractVec3f(a, ray.origin)) / detA;

    // if the values are negative it does not intersects or we cannot see the object
	if(detA == 0.0 || t < 0 || gamma < 0 || gamma > 1 - beta || beta < 0 || beta > 1 - gamma){
		return hit;
	}

	hit.hit = true;
	hit.type = 1;
	hit.objectID = obj_id;
	hit.materialID = triangle.material_id;
	hit.t = t;
	hit.intersectionPoint = FindIntersectionPoint(ray, t);
	hit.normal = CrossProduct(SubtractVec3f(b, a), SubtractVec3f(c, a));
	NormalizeVec3f(hit.normal);

	return hit;
}

Hit HitMesh(const Ray &ray, const parser::Mesh &mesh, const parser::Scene &scene, int obj_id){
    float t_min = std::numeric_limits<float>::max();
    Hit hit;
    hit.hit = false;
    for(int i = 0; i < mesh.faces.size(); i++){
        // every mesh faces is an triangle therefore we will use same algorithm which we used in triangle intersection.
        // we will get smallest positive t value.
        parser::Vec3f a = scene.vertex_data[mesh.faces[i].v0_id - 1];
        parser::Vec3f b = scene.vertex_data[mesh.faces[i].v1_id - 1];
        parser::Vec3f c = scene.vertex_data[mesh.faces[i].v2_id - 1];

        parser::Vec3f normal;
        normal = CrossProduct(SubtractVec3f(c, b), SubtractVec3f(a, b));
        NormalizeVec3f(normal);

        // some calculations for barycentric coordinates
        float detA = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), ray.direction);
        float beta = Determinant(SubtractVec3f(a, ray.origin), SubtractVec3f(a, c), ray.direction) / detA;
        float gamma = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, ray.origin), ray.direction) / detA;
        float t = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), SubtractVec3f(a, ray.origin)) / detA;

        // if the values are negative it does not intersects
        if(beta < 0 || beta > 1 - gamma) continue;
        if(gamma < 0 || gamma > 1 - beta) continue;
        // we cannot see the object
        if(t < 0) continue;

        if(t < t_min){
            hit.hit = true;
            
            hit.materialID = mesh.material_id;
            hit.normal = normal;
            t_min = t;
        }
    }

    if(hit.hit)
    {
        hit.t = t_min;
        hit.intersectionPoint = FindIntersectionPoint(ray, t_min);
    }

    hit.objectID = obj_id;
    hit.type = 2;

    return hit;
}


parser::Vec3f AmbientShading(const parser::Scene &scene, int materialID){
    parser::Vec3f result;
    result.x = scene.materials[materialID - 1].ambient.x * scene.ambient_light.x;
    result.y = scene.materials[materialID - 1].ambient.y * scene.ambient_light.y;
    result.z = scene.materials[materialID - 1].ambient.z * scene.ambient_light.z;
    return result;
}

Hit FindHit(const parser::Scene &scene, const Ray &ray)
{
    Hit hit;
    hit.hit = false;
    hit.t = std::numeric_limits<float>::max();

    // loop for sphere hit
    int spheresSize = scene.spheres.size();
    for (int sphereNumber = 0; sphereNumber < spheresSize; sphereNumber++)
    {
        parser::Sphere currentSphere = scene.spheres[sphereNumber];
        parser::Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
        float radius = currentSphere.radius;

        Hit hitCheck = HitSphere(ray, center, radius, currentSphere.material_id, sphereNumber);

        if(hitCheck.hit && hitCheck.t > 0 && hitCheck.t < hit.t)
        {
            hit.hit = true;
            hit.intersectionPoint = hitCheck.intersectionPoint;
            hit.materialID = hitCheck.materialID;
            hit.normal = hitCheck.normal;
            hit.objectID = hitCheck.objectID;
            hit.t = hitCheck.t;
            hit.type = hitCheck.type;
        }
    }
    
    // loop for triangle hit
    int trianglesSize = scene.triangles.size();
    for (int triangleNumber = 0; triangleNumber < trianglesSize; triangleNumber++)
    {
        parser::Triangle currentTriangle = scene.triangles[triangleNumber];
        Hit hitCheck = HitTriangle(ray, scene, currentTriangle, triangleNumber);

        if(hitCheck.hit && hitCheck.t && hitCheck.t < hit.t)
        {
            hit.hit = true;
            hit.intersectionPoint = hitCheck.intersectionPoint;
            hit.materialID = hitCheck.materialID;
            hit.normal = hitCheck.normal;
            hit.objectID = hitCheck.objectID;
            hit.t = hitCheck.t;
            hit.type = hitCheck.type;
        }
    }

    // loop for meshes hit
    int meshesSize = scene.meshes.size();
    for (int meshNumber = 0; meshNumber < meshesSize; meshNumber++)
    {
        parser::Mesh currentMesh = scene.meshes[meshNumber];
        Hit hitCheck = HitMesh(ray, currentMesh, scene, meshNumber);
        
        if(hitCheck.hit && hitCheck.t && hitCheck.t < hit.t)
        {
            hit.hit = true;
            hit.intersectionPoint = hitCheck.intersectionPoint;
            hit.materialID = hitCheck.materialID;
            hit.normal = hitCheck.normal;
            hit.objectID = hitCheck.objectID;
            hit.t = hitCheck.t;
            hit.type = hitCheck.type;
        }
    }
    
    return hit;
}

parser::Vec3f ComputeColor(const parser::Scene &scene, const parser::Camera &camera, const Hit &hit, const Ray &ray, int recursionCount)
{
    parser::Vec3f color;
    
    if(hit.hit)
    {
        color = AmbientShading(scene, hit.materialID);
    }
    else
    {
        color.x = 0;
        color.y = 0;
        color.z = 0;
    }
    return color;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int numberOfCameras = scene.cameras.size();

    for(int cameraNo = 0; cameraNo < numberOfCameras; cameraNo++)
    {
        int width = scene.cameras[cameraNo].image_width;
        int height = scene.cameras[cameraNo].image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        for(int j = 0; j < height; j++)
        {
            for(int i = 0; i < width; i++)
            {
                Ray ray = GenerateRay(scene.cameras[cameraNo], i, j);
                Hit hit = FindHit(scene, ray);
                parser::Vec3f color = ComputeColor(scene, scene.cameras[cameraNo], hit, ray, scene.max_recursion_depth);
                image[3 * j * width + 3 * i] = color.x;
                image[3 * j * width + 3 * i + 1] = color.y;
                image[3 * j * width + 3 * i + 2] = color.z;
            }
        }

        write_ppm(scene.cameras[cameraNo].image_name.c_str(), image, width, height);
    }
}
