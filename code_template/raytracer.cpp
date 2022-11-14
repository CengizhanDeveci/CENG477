#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <limits>
#include <cmath>
#include "ray.h"
#include <thread>

using namespace std;
using namespace parser;

typedef unsigned char RGB[3];

Hit hitSphere(const Ray &ray, Vec3f &center, float radius, int material_id, int object_id){
	Hit hit;

	float A = dot(ray.direction, ray.direction);
	float B = 2 * dot(ray.direction, subtract(ray.origin, center));
	float C = dot(subtract(ray.origin, center), subtract(ray.origin, center)) - radius * radius;

	float delta = B * B - 4 * A * C;

	if(delta < 0){ // no intersection
		hit.hit_happened = false;
	}
	else{ // intersects
		float t1 = (-1 * B + sqrtf(delta))/2*A;
		float t2 = (-1 * B - sqrtf(delta))/2*A;

        hit.hit_happened = true;

		hit.material_id = material_id;
		hit.type = 0;
		hit.object_id = object_id;

		float t = (t1 > t2) ? t2 : t1;
        hit.t = t;
		hit.intersection_point = findIntersectionPoint(ray, t);

		hit.normal = subtract(hit.intersection_point, center);
		hit.normal.x /= radius;
		hit.normal.y /= radius;
		hit.normal.z /= radius;
        hit.normal = normalize(hit.normal);
	}

	return hit;
}

Hit hitTriangle(const Ray &ray, const Scene &scene, Triangle &triangle, int obj_id){
	Hit hit;
	hit.hit_happened = false;

    Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
    Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
    Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];

    Vec3f normal;
    normal = cross(subtract(c, b),subtract(a, b));
    normal = normalize(normal);

    // some calculations for barycentric coordinates
    float detA = determinant(subtract(a, b), subtract(a, c), ray.direction);
    float beta = determinant(subtract(a, ray.origin), subtract(a, c), ray.direction) / detA;
    float gamma = determinant(subtract(a, b), subtract(a, ray.origin), ray.direction) / detA;
    float t = determinant(subtract(a, b), subtract(a, c), subtract(a, ray.origin)) / detA;

    // if the values are negative it does not intersects or we cannot see the object
	if(detA == 0.0 || t < 0 || gamma < 0 || gamma > 1 - beta || beta < 0 || beta > 1 - gamma){
		return hit;
	}

	hit.hit_happened = true;
	hit.type = 1;
	hit.object_id = obj_id;
	hit.material_id = triangle.material_id;
	hit.t = t;
	hit.intersection_point = findIntersectionPoint(ray, t);
	hit.normal = cross(subtract(b, a), subtract(c, a));
	hit.normal = normalize(hit.normal);

	return hit;
}

Hit hitMesh(const Ray &ray, const Mesh &mesh, const Scene &scene, int obj_id){
    float t_min = numeric_limits<float>::max();
    Hit hit;
    hit.hit_happened = false;
    for(int i = 0; i < mesh.faces.size(); i++){
        // every mesh faces is an triangle therefore we will use same algorithm which we used in triangle intersection.
        // we will get smallest positive t value.
        Vec3f a = scene.vertex_data[mesh.faces[i].v0_id - 1];
        Vec3f b = scene.vertex_data[mesh.faces[i].v1_id - 1];
        Vec3f c = scene.vertex_data[mesh.faces[i].v2_id - 1];

        Vec3f normal;
        normal = cross(subtract(c, b), subtract(a, b));
        normal = normalize(normal);

        // some calculations for barycentric coordinates
        float detA = determinant(subtract(a, b), subtract(a, c), ray.direction);
        float beta = determinant(subtract(a, ray.origin), subtract(a, c), ray.direction) / detA;
        float gamma = determinant(subtract(a, b), subtract(a, ray.origin), ray.direction) / detA;
        float t = determinant(subtract(a, b), subtract(a, c), subtract(a, ray.origin)) / detA;

        // if the values are negative it does not intersects
        if(beta < 0 || beta > 1 - gamma) continue;
        if(gamma < 0 || gamma > 1 - beta) continue;
        // we cannot see the object
        if(t < 0) continue;

        if(t < t_min){
            hit.hit_happened = true;
            hit.intersection_point = findIntersectionPoint(ray, t);
            hit.material_id = mesh.material_id;
            hit.normal = normal;
            hit.object_id = obj_id;
            hit.type = 2;
            hit.t = t;
        }
    }

    return hit;
}

Vec3f ambientShading(Scene scene, int material_id){
    Vec3f result;
    result.x = scene.materials[material_id - 1].ambient.x * scene.ambient_light.x;
    result.y = scene.materials[material_id - 1].ambient.y * scene.ambient_light.y;
    result.z = scene.materials[material_id - 1].ambient.z * scene.ambient_light.z;
    return result;
}

const Vec3f diffuseShading(const PointLight &current_light, const Scene &scene, Hit hit){
	Vec3f color;
	Vec3f irradiance = findIrradiance(current_light, hit.intersection_point);
	Vec3f w_i = subtract(current_light.position, hit.intersection_point);
	w_i = normalize(w_i);

	float dotPro = dot(w_i, hit.normal);
    dotPro = (dotPro < 0) ? 0 : dotPro;

	color.x = scene.materials[hit.material_id - 1].diffuse.x * dotPro * irradiance.x;
	color.y = scene.materials[hit.material_id - 1].diffuse.y * dotPro * irradiance.y;
	color.z = scene.materials[hit.material_id - 1].diffuse.z * dotPro * irradiance.z;
	
	return color;
}

Vec3f specularShading(const PointLight &current_light, const Scene &scene, const Ray &ray, Hit hit){
	Vec3f color;
	Vec3f irradiance = findIrradiance(current_light, hit.intersection_point);

	Vec3f wi = subtract(current_light.position, hit.intersection_point);
	wi = normalize(wi);

	Vec3f h = subtract(wi, ray.direction);
	h = normalize(h);

	float cos_alpha_prime = dot(hit.normal, h);
	cos_alpha_prime = (cos_alpha_prime < 0) ? 0 : cos_alpha_prime;

	color.x = scene.materials[hit.material_id - 1].specular.x * pow(cos_alpha_prime, scene.materials[hit.material_id - 1].phong_exponent) * irradiance.x;
	color.y = scene.materials[hit.material_id - 1].specular.y * pow(cos_alpha_prime, scene.materials[hit.material_id - 1].phong_exponent) * irradiance.y;
	color.z = scene.materials[hit.material_id - 1].specular.z * pow(cos_alpha_prime, scene.materials[hit.material_id - 1].phong_exponent) * irradiance.z;

	return color;
}

// it will find object which intersects with ray first
Hit findHit(const Scene &scene, const Ray &ray){
    Hit hitted;
    hitted.hit_happened = false;
    hitted.t = numeric_limits<float>::max();
    // loop for sphere hit
    for(int sphere_number = 0; sphere_number < scene.spheres.size(); sphere_number++){
        Sphere current_sphere = scene.spheres[sphere_number];
        Vec3f center = scene.vertex_data[current_sphere.center_vertex_id - 1];
        float radius = current_sphere.radius;

        Hit hit = hitSphere(ray, center, radius, current_sphere.material_id, sphere_number);

        if(hit.hit_happened && hit.t >= 0)
        {
            if(hit.t < hitted.t){
                hitted.hit_happened = true;
                hitted.intersection_point = hit.intersection_point;
                hitted.material_id = hit.material_id;
                hitted.normal = hit.normal;
                hitted.object_id = hit.object_id;
                hitted.t = hit.t;
                hitted.type = hit.type;
            }
        }
    }
    // loop for triangle hit
    for(int triangle_number = 0; triangle_number < scene.triangles.size(); triangle_number++){
        Triangle current_triangle = scene.triangles[triangle_number];
        Vec3f v0 = scene.vertex_data[current_triangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[current_triangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[current_triangle.indices.v2_id - 1];

        Hit hit = hitTriangle(ray, scene, current_triangle, triangle_number);
        if(hit.hit_happened && hit.t >= 0)
        {
            if(hit.t < hitted.t){
                hitted.hit_happened = true;
                hitted.intersection_point = hit.intersection_point;
                hitted.material_id = hit.material_id;
                hitted.normal = hit.normal;
                hitted.object_id = hit.object_id;
                hitted.t = hit.t;
                hitted.type = hit.type;
            }
        }
    }
    // loop for mesh hit
    for(int mesh_number = 0; mesh_number < scene.meshes.size(); mesh_number++){
        Mesh current_mesh = scene.meshes[mesh_number];
        Hit hit = hitMesh(ray, current_mesh, scene, mesh_number);
        if(hit.hit_happened && hit.t >= 0)
        {
            if(hit.t < hitted.t){
                hitted.hit_happened = true;
                hitted.intersection_point = hit.intersection_point;
                hitted.material_id = hit.material_id;
                hitted.normal = hit.normal;
                hitted.object_id = hit.object_id;
                hitted.t = hit.t;
                hitted.type = hit.type;
            }
        }
    }

    return hitted;
}

bool shadow(Ray shadow_ray, Scene scene, float t){
    // loop for sphere hit
    for(int sphereNumber = 0; sphereNumber < scene.spheres.size(); sphereNumber++){
        Sphere currentSphere = scene.spheres[sphereNumber];
        Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
        float radius = currentSphere.radius;
        Hit shadow_hit = hitSphere(shadow_ray, center, radius, currentSphere.material_id, sphereNumber);

        if(shadow_hit.hit_happened && t > shadow_hit.t && shadow_hit.t >= 0) return true;
    }

    // loop for triangle hit
    for(int triangleNumber = 0; triangleNumber < scene.triangles.size(); triangleNumber++)
    {
        Triangle currentTriangle = scene.triangles[triangleNumber];
        Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];
        Hit shadow_hit = hitTriangle(shadow_ray, scene, currentTriangle, triangleNumber);

        if(shadow_hit.hit_happened && t > shadow_hit.t && shadow_hit.t >= 0) return true;
    }

    // loop for mesh hit
    for(int meshNumber = 0; meshNumber < scene.meshes.size(); meshNumber++)
    {
        Mesh currentMesh = scene.meshes[meshNumber];
        Hit shadow_hit = hitMesh(shadow_ray, currentMesh, scene, meshNumber);

        if(shadow_hit.hit_happened && t > shadow_hit.t && shadow_hit.t >= 0) return true;
    }
}

Vec3f computeColor(const Scene &scene, Hit &hit, const Camera &camera, Ray &ray, int recursion_count){
    Vec3f color;

    if(hit.hit_happened){
        int material_id = hit.material_id;

        color = ambientShading(scene, material_id);


        for(int light_number = 0; light_number < scene.point_lights.size(); light_number++){
            bool is_shadow = false;

            PointLight current_light = scene.point_lights[light_number];
            float light_to_cam = findDistance(current_light.position, hit.intersection_point);

            Vec3f wi = subtract(current_light.position, hit.intersection_point);
			wi = normalize(wi);

            Vec3f wi_epsilon;
			wi_epsilon.x = wi.x * scene.shadow_ray_epsilon;
			wi_epsilon.y = wi.y * scene.shadow_ray_epsilon;
			wi_epsilon.z = wi.z * scene.shadow_ray_epsilon;

            Ray shadow_ray;
            shadow_ray.direction = wi;
            shadow_ray.origin = add(hit.intersection_point,wi_epsilon);

            float t_light = subtract(current_light.position, shadow_ray.origin).x / shadow_ray.direction.x;
            is_shadow = shadow(shadow_ray, scene, t_light);

            // if no shadowRay intersection but there is object-ray intersection
	        if(!is_shadow || (is_shadow && light_to_cam == 0))
	        {

		        Vec3f diffuse = diffuseShading(current_light, scene, hit);
				
		        Vec3f specular = specularShading(current_light, scene, ray, hit);
		                    
              	color.x += diffuse.x + specular.x;
  		        color.y += diffuse.y + specular.y;
  		        color.z += diffuse.z + specular.z;
            
	        }
        }

        // checks the mirror image
        bool is_mirror = scene.materials[hit.material_id - 1].is_mirror;
        Vec3f mirror_color;

        if(is_mirror && recursion_count > 0){
            float w_i = -2 * dot(ray.direction, hit.normal);
            Vec3f normal_wi;
            hit.normal = normalize(hit.normal);
            ray.direction = normalize(ray.direction);
            normal_wi.x = hit.normal.x * w_i + ray.direction.x;
            normal_wi.y = hit.normal.y * w_i + ray.direction.y;
            normal_wi.z = hit.normal.z * w_i + ray.direction.z;
            
            normal_wi = normalize(normal_wi);

            Vec3f wi_epsilon;
            wi_epsilon.x = normal_wi.x * scene.shadow_ray_epsilon;
            wi_epsilon.y = normal_wi.y * scene.shadow_ray_epsilon;
            wi_epsilon.z = normal_wi.z * scene.shadow_ray_epsilon;

            Ray reflection_ray;
            reflection_ray.origin = add(hit.intersection_point, wi_epsilon);
            reflection_ray.direction = normal_wi;

            Hit new_hit = findHit(scene, reflection_ray);
            if(!(new_hit.type == hit.type && new_hit.object_id == hit.object_id)){
                mirror_color = computeColor(scene, new_hit, camera, reflection_ray, recursion_count - 1);
                color.x += mirror_color.x * scene.materials[hit.material_id - 1].mirror.x;
                color.y += mirror_color.y * scene.materials[hit.material_id - 1].mirror.y;
                color.z += mirror_color.z * scene.materials[hit.material_id - 1].mirror.z;
            }
        }
    }else{
        color.x = scene.background_color.x;
        color.y = scene.background_color.y;
        color.z = scene.background_color.z;
    }

    return color;
}

// void multiThread(Scene scene, int camera_no, unsigned char* &image, int height, int width){
//     for(int i = 0; i < height; i++){
//         for(int j = 0; j < width; j++){
//             Ray ray = generateRay(scene.cameras[camera_no], j, i);
//             Vec3f color;

//             Hit hit = findHit(scene, ray);
//             color = computeColor(scene, hit, scene.cameras[camera_no],ray, scene.max_recursion_depth);

//             // if the color values greater than 255 we set 255 otherwise round the nearest integer value
//             if(color.x > 255) image[3* i * width + 3 * j] = 255;
//             else image[3* i * width + 3 * j] = (int) (color.x + 0.5);

//             if(color.y > 255) image[3 * i * width + 3 * j + 1] = 255;
//             else image[3* i * width + 3 * j + 1] = (int) (color.y + 0.5);

//             if(color.z > 255) image[3 * i * width + 3 * j + 2] = 255;
//             else image[3 * i * width + 3 * j + 2] = (int) (color.z + 0.5);
//         }
//     }
// }

int main(int argc, char* argv[]){
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    
    int number_of_cameras = scene.cameras.size();

    for(int camera_no = 0; camera_no < number_of_cameras; camera_no++){
        int width = scene.cameras[camera_no].image_width;
        int height = scene.cameras[camera_no].image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        //thread(multiThread, scene, camera_no, image, height, width);

        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                Ray ray = generateRay(scene.cameras[camera_no], j, i);
                Vec3f color;

                Hit hit = findHit(scene, ray);
                color = computeColor(scene, hit, scene.cameras[camera_no],ray, scene.max_recursion_depth);

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

