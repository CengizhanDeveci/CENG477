#include <iostream>
#include <cmath>
#include "parser.h"
#include "ppm.h"
#include "ray.h"
#include <limits>
#include <thread>
#include "bvh.h"
#include "raytracer_helper.h"

using namespace ray;

typedef unsigned char RGB[3];

std::vector<Object*> objects;

parser::Vec3f FindIntersectionPoint(const Ray &ray, float t)
{
    parser::Vec3f result;
    result.x = ray.origin.x + t * ray.direction.x;
    result.y = ray.origin.y + t * ray.direction.y;
    result.z = ray.origin.z + t * ray.direction.z;

    return result;
}

bool HitSphere(const Ray &ray, parser::Vec3f &center, float radius, int materialID, int objectID, Hit& hit)
{
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
        if(t < 0)
        {
            hit.hit = false;
            return false;
        }
        hit.t = t;
		hit.intersectionPoint = FindIntersectionPoint(ray, t);

		hit.normal = SubtractVec3f(hit.intersectionPoint, center);
		hit.normal.x /= radius;
		hit.normal.y /= radius;
		hit.normal.z /= radius;
        hit.normal = ray::NormalizeVec3f(hit.normal);

    }
    return hit.hit;
}

bool HitTriangle(const Ray &ray, parser::Vec3f& a, parser::Vec3f& b, parser::Vec3f& c, int obj_id, Hit& hit){
	hit.hit = false;

    parser::Vec3f normal;
    normal = CrossProduct(SubtractVec3f(c, b), SubtractVec3f(a, b));
    normal = ray::NormalizeVec3f(normal);

    // some calculations for barycentric coordinates
    float detA = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), ray.direction);
    float beta = Determinant(SubtractVec3f(a, ray.origin), SubtractVec3f(a, c), ray.direction) / detA;
    float gamma = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, ray.origin), ray.direction) / detA;
    float t = Determinant(SubtractVec3f(a, b), SubtractVec3f(a, c), SubtractVec3f(a, ray.origin)) / detA;

    // if the values are negative it does not intersects or we cannot see the object
	if(detA == 0.0 || t < 0 || gamma < 0 || gamma > 1 - beta || beta < 0 || beta > 1 - gamma){
		return false;
	}

	hit.hit = true;
	hit.type = 1;
	hit.objectID = obj_id;
	hit.t = t;
	hit.intersectionPoint = FindIntersectionPoint(ray, t);

	return true;
}

Hit HitTriangle(const Ray &ray, const parser::Scene &scene, parser::Triangle &triangle, int obj_id){
	Hit hit;
	hit.hit = false;

    parser::Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
    parser::Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
    parser::Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];

    parser::Vec3f normal;
    normal = CrossProduct(SubtractVec3f(c, b), SubtractVec3f(a, b));
    normal = ray::NormalizeVec3f(normal);

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
	hit.normal = ray::NormalizeVec3f(hit.normal);

	return hit;
}

parser::Vec3f FindIrradiance(const parser::PointLight &pointLight, const parser::Vec3f &intersectionPoint)
{
    parser::Vec3f result;
    parser::Vec3f d = SubtractVec3f(pointLight.position, intersectionPoint);
    float dDotD = DotProduct(d, d);

    if(dDotD <= 0) result;

    result.x = pointLight.intensity.x / dDotD;
    result.y = pointLight.intensity.y / dDotD;
    result.z = pointLight.intensity.z / dDotD;
    
    return result;
}


parser::Vec3f AmbientShading(const parser::Scene &scene, int materialID){
    parser::Vec3f result;
    result.x = scene.materials[materialID - 1].ambient.x * scene.ambient_light.x;
    result.y = scene.materials[materialID - 1].ambient.y * scene.ambient_light.y;
    result.z = scene.materials[materialID - 1].ambient.z * scene.ambient_light.z;
    return result;
}

parser::Vec3f DiffuseShading(const parser::PointLight &pointLight, const parser::Scene &scene, const Hit &hit)
{
    parser::Vec3f color;
    parser::Vec3f irradiance = FindIrradiance(pointLight, hit.intersectionPoint);
    parser::Vec3f wi = SubtractVec3f(pointLight.position, hit.intersectionPoint);
    wi = ray::NormalizeVec3f(wi);

    float cosThetaPrime = DotProduct(wi, hit.normal);
    if(cosThetaPrime <= 0) {color.x = 0; color.y=0; color.z=0;return color;}

    color.x = scene.materials[hit.materialID - 1].diffuse.x * cosThetaPrime * irradiance.x;
    color.y = scene.materials[hit.materialID - 1].diffuse.y * cosThetaPrime * irradiance.y;
    color.z = scene.materials[hit.materialID - 1].diffuse.z * cosThetaPrime * irradiance.z;

    return color;
}

parser::Vec3f SpecularShading(const parser::PointLight &pointLight, const parser::Scene &scene, const Ray &ray, const Hit &hit)
{
    parser::Vec3f color;
    parser::Vec3f irradiance = FindIrradiance(pointLight, hit.intersectionPoint);

    parser::Vec3f wi = SubtractVec3f(pointLight.position, hit.intersectionPoint);
	wi = ray::NormalizeVec3f(wi);

	parser::Vec3f h = SubtractVec3f(wi, ray.direction);
	h = ray::NormalizeVec3f(h);

	float cos_alpha_prime = DotProduct(hit.normal, h);
	if(cos_alpha_prime <= 0) {color.x = 0; color.y=0; color.z=0;return color;}

    color.x = scene.materials[hit.materialID - 1].specular.x * pow(cos_alpha_prime, scene.materials[hit.materialID - 1].phong_exponent) * irradiance.x;
	color.y = scene.materials[hit.materialID - 1].specular.y * pow(cos_alpha_prime, scene.materials[hit.materialID - 1].phong_exponent) * irradiance.y;
	color.z = scene.materials[hit.materialID - 1].specular.z * pow(cos_alpha_prime, scene.materials[hit.materialID - 1].phong_exponent) * irradiance.z;

    return color;
}

bool ShadowCheck(const parser::Scene &scene, Ray shadowRay, float t, BVH& bvh)
{
    Hit hit;
    hit.hit = bvh.Intersect(shadowRay, bvh.root, hit);
    if(hit.hit && hit.t < t)
    {
        return true;
    }
    return false;
}

parser::Vec3f ComputeColor(const parser::Scene &scene, const parser::Camera &camera, const Hit &hit, const Ray &ray, int recursionCount, BVH& bvh)
{
    parser::Vec3f color;
    
    if(hit.hit)
    {
        color = AmbientShading(scene, hit.materialID);
        for(int lightNumber = 0; lightNumber < scene.point_lights.size(); lightNumber++)
        {
            parser::Vec3f wiEpsilon;
            wiEpsilon.x = hit.normal.x * scene.shadow_ray_epsilon;
            wiEpsilon.y = hit.normal.y * scene.shadow_ray_epsilon;
            wiEpsilon.z = hit.normal.z * scene.shadow_ray_epsilon;

            Ray shadowRay;
            shadowRay.origin = AddVec3f(wiEpsilon, hit.intersectionPoint);
            parser::Vec3f tmp1 = hit.intersectionPoint;
            parser::Vec3f tmp2 = scene.point_lights[lightNumber].position;
            shadowRay.direction = SubtractVec3f(tmp2, tmp1);
            shadowRay.direction = ray::NormalizeVec3f(shadowRay.direction);
            float tLight = (tmp2.x - shadowRay.origin.x) / shadowRay.direction.x;
            bool isShadow = ShadowCheck(scene, shadowRay, tLight, bvh);

            if(!isShadow)
            {
                parser::Vec3f diffuse = DiffuseShading(scene.point_lights[lightNumber], scene, hit);
                parser::Vec3f specular = SpecularShading(scene.point_lights[lightNumber], scene, ray, hit);

                color.x += diffuse.x + specular.x;
                color.y += diffuse.y + specular.y;
                color.z += diffuse.z + specular.z;
            }
        }

        if(scene.materials[hit.materialID - 1].is_mirror && recursionCount > 0)
        {
            parser::Vec3f mirrorColor;

            float cosTheta = DotProduct(ray.direction, hit.normal);
            parser::Vec3f wr;
            wr.x = -2 * hit.normal.x * cosTheta + ray.direction.x;
            wr.y = -2 * hit.normal.y * cosTheta + ray.direction.y;
            wr.z = -2 * hit.normal.z * cosTheta + ray.direction.z;
            wr = ray::NormalizeVec3f(wr);
            
            parser::Vec3f wiEpsilon;
            wiEpsilon.x = wr.x * scene.shadow_ray_epsilon;
            wiEpsilon.y = wr.y * scene.shadow_ray_epsilon;
            wiEpsilon.z = wr.z * scene.shadow_ray_epsilon;

            Ray reflectionRay;
            reflectionRay.origin = AddVec3f(hit.intersectionPoint, wiEpsilon);
            reflectionRay.direction = wr;

            Hit newHit;
            newHit.hit = false;
            bool hitHappenedMirror = bvh.Intersect(reflectionRay, bvh.root, newHit);
            if(newHit.hit)
            {
                mirrorColor = ComputeColor(scene, camera, newHit, reflectionRay, recursionCount - 1, bvh);
                color.x += mirrorColor.x * scene.materials[hit.materialID - 1].mirror.x;
                color.y += mirrorColor.y * scene.materials[hit.materialID - 1].mirror.y;
                color.z += mirrorColor.z * scene.materials[hit.materialID - 1].mirror.z;
            }
            
        }
    }
    else
    {
        color.x = scene.background_color.x;
        color.y = scene.background_color.y;
        color.z = scene.background_color.z;
    }
    return color;
}

void multiThread(parser::Scene scene, int cameraNo, unsigned char* &image, int height, int width, int start, BVH& bvh)
{
    for(int j = start; j < height; j+=4)
    {
        for(int i = 0; i < width; i++)
        {   
            Ray ray = GenerateRay(scene.cameras[cameraNo], i, j);
            Hit hit;
            hit.hit = false;
            bool hitHappened = bvh.Intersect(ray, bvh.root, hit);

            parser::Vec3f color = ComputeColor(scene, scene.cameras[cameraNo], hit, ray, scene.max_recursion_depth, bvh);
            image[3 * j * width + 3 * i] = color.x < 255 ? color.x : 255;
            image[3 * j * width + 3 * i + 1] = color.y < 255 ? color.y : 255;
            image[3 * j * width + 3 * i + 2] = color.z < 255 ? color.z : 255;
        }
    }
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for(parser::Triangle& triangle : scene.triangles)
    {
        triangle.indices.boundingBox = new BoundingBox();
        triangle.indices.boundingBox->min.x = std::min(scene.vertex_data[triangle.indices.v0_id - 1].x, std::min(scene.vertex_data[triangle.indices.v1_id - 1].x, scene.vertex_data[triangle.indices.v2_id - 1].x));
        triangle.indices.boundingBox->min.y = std::min(scene.vertex_data[triangle.indices.v0_id - 1].y, std::min(scene.vertex_data[triangle.indices.v1_id - 1].y, scene.vertex_data[triangle.indices.v2_id - 1].y));
        triangle.indices.boundingBox->min.z = std::min(scene.vertex_data[triangle.indices.v0_id - 1].z, std::min(scene.vertex_data[triangle.indices.v1_id - 1].z, scene.vertex_data[triangle.indices.v2_id - 1].z));
        triangle.indices.boundingBox->max.x = std::max(scene.vertex_data[triangle.indices.v0_id - 1].x, std::max(scene.vertex_data[triangle.indices.v1_id - 1].x, scene.vertex_data[triangle.indices.v2_id - 1].x));
        triangle.indices.boundingBox->max.y = std::max(scene.vertex_data[triangle.indices.v0_id - 1].y, std::max(scene.vertex_data[triangle.indices.v1_id - 1].y, scene.vertex_data[triangle.indices.v2_id - 1].y));
        triangle.indices.boundingBox->max.z = std::max(scene.vertex_data[triangle.indices.v0_id - 1].z, std::max(scene.vertex_data[triangle.indices.v1_id - 1].z, scene.vertex_data[triangle.indices.v2_id - 1].z));

        triangle.indices.a = scene.vertex_data[triangle.indices.v0_id - 1];
        triangle.indices.b = scene.vertex_data[triangle.indices.v1_id - 1];
        triangle.indices.c = scene.vertex_data[triangle.indices.v2_id - 1];

        parser::Vec3f normal;
        normal = CrossProduct(SubtractVec3f(triangle.indices.c, triangle.indices.b), SubtractVec3f(triangle.indices.a, triangle.indices.b));
        normal = ray::NormalizeVec3f(normal);   
        triangle.indices.normal = normal;
        triangle.indices.material_id = triangle.material_id;
        triangle.indices.isSphere = false;

        objects.push_back(&triangle.indices);   
    }
    for (parser::Sphere& sphere : scene.spheres)
    {
        sphere.boundingBox = new BoundingBox();
        sphere.boundingBox->min = SubtractVec3f(scene.vertex_data[sphere.center_vertex_id - 1], parser::Vec3f(sphere.radius, sphere.radius, sphere.radius));
        sphere.boundingBox->max = AddVec3f(scene.vertex_data[sphere.center_vertex_id - 1], parser::Vec3f(sphere.radius, sphere.radius, sphere.radius));
        sphere.normal = parser::Vec3f(0, 0, 0);
        sphere.isSphere = true;
        sphere.center = scene.vertex_data[sphere.center_vertex_id - 1];
        sphere.materialID = sphere.material_id;
        sphere.tmp = sphere.radius;
        objects.push_back(&sphere);
    }
    for (parser::Mesh& mesh : scene.meshes)
    {
        for(parser::Face& face : mesh.faces)
        {
            face.boundingBox = new BoundingBox();
            face.boundingBox->min.x = std::min(scene.vertex_data[face.v0_id - 1].x, std::min(scene.vertex_data[face.v1_id - 1].x, scene.vertex_data[face.v2_id - 1].x));
            face.boundingBox->min.y = std::min(scene.vertex_data[face.v0_id - 1].y, std::min(scene.vertex_data[face.v1_id - 1].y, scene.vertex_data[face.v2_id - 1].y));
            face.boundingBox->min.z = std::min(scene.vertex_data[face.v0_id - 1].z, std::min(scene.vertex_data[face.v1_id - 1].z, scene.vertex_data[face.v2_id - 1].z));
            face.boundingBox->max.x = std::max(scene.vertex_data[face.v0_id - 1].x, std::max(scene.vertex_data[face.v1_id - 1].x, scene.vertex_data[face.v2_id - 1].x));
            face.boundingBox->max.y = std::max(scene.vertex_data[face.v0_id - 1].y, std::max(scene.vertex_data[face.v1_id - 1].y, scene.vertex_data[face.v2_id - 1].y));
            face.boundingBox->max.z = std::max(scene.vertex_data[face.v0_id - 1].z, std::max(scene.vertex_data[face.v1_id - 1].z, scene.vertex_data[face.v2_id - 1].z));

            face.a = scene.vertex_data[face.v0_id - 1];
            face.b = scene.vertex_data[face.v1_id - 1];
            face.c = scene.vertex_data[face.v2_id - 1];

            parser::Vec3f normal;
            normal = CrossProduct(SubtractVec3f(face.c, face.b), SubtractVec3f(face.a, face.b));
            normal = ray::NormalizeVec3f(normal);
            face.normal = normal;
            face.material_id = mesh.material_id;
            face.isSphere = false;

            objects.push_back(&face);
        }
    }
    BVH bvh(objects);

    int numberOfCameras = scene.cameras.size();

    for(int cameraNo = 0; cameraNo < numberOfCameras; cameraNo++)
    {

        int width = scene.cameras[cameraNo].image_width;
        int height = scene.cameras[cameraNo].image_height;
        unsigned char* image = new unsigned char [width * height * 3];
        std::thread t1(multiThread, scene, cameraNo, std::ref(image),height, width, 0, std::ref(bvh));
        std::thread t2(multiThread, scene, cameraNo, std::ref(image),height, width, 1, std::ref(bvh));
        std::thread t3(multiThread, scene, cameraNo, std::ref(image),height, width, 2, std::ref(bvh));
        std::thread t4(multiThread, scene, cameraNo, std::ref(image),height, width, 3, std::ref(bvh));

        t1.join();
        t2.join();
        t3.join();
        t4.join();

        write_ppm(scene.cameras[cameraNo].image_name.c_str(), image, width, height);
    }
}