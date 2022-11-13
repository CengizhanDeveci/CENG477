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
};

class Hit{
        public:
        bool hit_happened;
        Vec3f intersection_point;
        Vec3f normal;
        int material_id;
        float t;
        int type;
        int object_id;
};

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
    float left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float bottom = camera.near_plane.z;
    float top = camera.near_plane.w;

    Vec3f gaze = normalize(camera.gaze);

    float su = (right - left)*(j + 0.5)/camera.image_width;
    float sv = (top - bottom)*(i + 0.5)/camera.image_height;

    Vec3f m, q, u, v;
    m.x = camera.position.x + (gaze.x * camera.near_distance);
    m.y = camera.position.y + (gaze.y * camera.near_distance);
    m.z = camera.position.z + (gaze.z * camera.near_distance);

    u = cross(gaze, camera.up);
    u = normalize(u);

    v = cross(u, gaze);

    q.x = m.x + (u.x*left) + (v.x*top);
    q.y = m.y + (u.y*left) + (v.y*top);
    q.z = m.z + (u.z*left) + (v.z*top);

    Vec3f s;
    s.x = q.x + (u.x*su) - (v.x * sv);
    s.y = q.y + (u.y*su) - (v.y * sv);
    s.z = q.z + (u.z*su) - (v.z * sv);

    Ray ray;
    ray.origin = camera.position;
    ray.direction = subtract(s, camera.position);
    ray.direction = normalize(ray.direction);
    ray.shadowRay = false;

    return ray;
}

Hit sphereIntersection(const Ray &ray, const Vec3f &center, float radius, int material_id, int obj_id)
{
	Hit hit;

	const float A = dot(ray.direction, ray.direction);
	Vec3f e_minus_c = subtract(ray.origin, center);
	const float B = 2 * dot(ray.direction, e_minus_c);
	const float C = dot(e_minus_c, e_minus_c) - radius * radius;

	const float discriminant = B*B - 4*A*C;

	if(discriminant < 0)			// no intersection
	{
		hit.hit_happened = false;
	}
	else							// intersection at 1 or 2 points 
	{	
		const float t1 = (-1 * B + sqrtf(discriminant))/2*A;
		const float t2 = (-1 * B - sqrtf(discriminant))/2*A;

		hit.material_id = material_id;
		hit.hit_happened = true;
		hit.type = 0;
		hit.object_id = obj_id;

		const float t = fmin(t1, t2);
		hit.intersection_point = findIntersectionPoint(ray, t);
		hit.normal = subtract(hit.intersection_point, center);
		hit.normal.x /= radius;
		hit.normal.y /= radius;
		hit.normal.z /= radius;
        hit.normal = normalize(hit.normal);
		hit.t = t;
	}

	return hit;
}

Hit triangleIntersection(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, int material_id, int obj_id)
{
	Hit hit;
	hit.hit_happened = false;

	Vec3f o = ray.origin;
	Vec3f d = ray.direction;

	Vec3f a_minus_b = subtract(a, b);
	Vec3f a_minus_c = subtract(a, c);
	Vec3f a_minus_o = subtract(a, o);

	float detA = determinant(a_minus_b, a_minus_c, d);
	if(detA == 0.0)
	{
		return hit;
	}

	float t = (determinant(a_minus_b, a_minus_c, a_minus_o))/detA;
	if(t <= 0.0) {
		return hit;
	}

	float gamma = (determinant(a_minus_b,a_minus_o, d))/detA;
	if(gamma < 0 || gamma > 1) {
		return hit;
	}

	float beta = (determinant(a_minus_o, a_minus_c, d))/detA;
	if(beta < 0 || beta > (1 - gamma)) {
		return hit;
	}

	hit.hit_happened = true;
	hit.type = 1;
	hit.object_id = obj_id;
	hit.material_id = material_id;
	hit.t = t;
	hit.intersection_point = findIntersectionPoint(ray, t);
	hit.normal = cross(subtract(b, a), subtract(c, a));
	hit.normal = normalize(hit.normal);

	return hit;
}

Hit meshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int material_id, int obj_id)
{
    float t_min = numeric_limits<float>::max();
    Hit hit;
    hit.hit_happened = false;
    for(int i = 0; i < mesh.faces.size(); i++){
        Vec3f a = scene.vertex_data[mesh.faces[i].v0_id - 1];
        Vec3f b = scene.vertex_data[mesh.faces[i].v1_id - 1];
        Vec3f c = scene.vertex_data[mesh.faces[i].v2_id - 1];

        Vec3f normal;
        normal = cross(subtract(c, b), subtract(a, b));
        normal = normalize(normal);

        float detA = determinant(subtract(a, b), subtract(a, c), ray.direction);
        float beta = determinant(subtract(a, ray.origin), subtract(a, c), ray.direction) / detA;
        float gamma = determinant(subtract(a, b), subtract(a, ray.origin), ray.direction) / detA;
        float t = determinant(subtract(a, b), subtract(a, c), subtract(a, ray.origin)) / detA;

        if(beta < 0 || beta > 1 - gamma) continue;
        if(gamma < 0 || gamma > 1 - beta) continue;

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

Vec3f findIrradiance(const PointLight &point_light, const Vec3f &intersection_point){
    Vec3f irradiance;
    Vec3f d = subtract(point_light.position, intersection_point);
    float d_square = dot(d, d);

    if(d_square != 0.0){
	    irradiance.x = point_light.intensity.x/d_square;
	    irradiance.y = point_light.intensity.y/d_square;
	    irradiance.z = point_light.intensity.z/d_square;
    }
    return irradiance;
}

const Vec3f findDiffuse(const PointLight &current_light, const Scene &scene, int material_id, const Vec3f &normal, const Vec3f &intersection_point)
{
	Vec3f diffuse;

	Vec3f irradiance = findIrradiance(current_light, intersection_point);

	Vec3f l = subtract(current_light.position, intersection_point);
	l = normalize(l);


	float dotPro = dot(l, normal);
	if(dotPro < 0)
	{
		dotPro = 0;
	}

	diffuse.x = scene.materials[material_id - 1].diffuse.x * dotPro * irradiance.x;
	diffuse.y = scene.materials[material_id - 1].diffuse.y * dotPro * irradiance.y;
	diffuse.z = scene.materials[material_id - 1].diffuse.z * dotPro * irradiance.z;
	
	return diffuse;
}

Vec3f findSpecular(const PointLight &current_light, const Scene &scene, const Ray &ray, int material_id, const Vec3f &normal, const Vec3f &intersection_point){
	Vec3f specular;

	Material material = scene.materials[material_id - 1];

	Vec3f irradiance = findIrradiance(current_light, intersection_point);

	Vec3f wi = subtract(current_light.position, intersection_point);
	wi = normalize(wi);

	Vec3f h = subtract(wi, ray.direction);
	h = normalize(h);

	float dotPro = dot(normal, h);
	if(dotPro < 0){
		dotPro = 0;
	}

	specular.x = material.specular.x * pow(dotPro, material.phong_exponent) * irradiance.x;
	specular.y = material.specular.y * pow(dotPro, material.phong_exponent) * irradiance.y;
	specular.z = material.specular.z * pow(dotPro, material.phong_exponent) * irradiance.z;

	return specular;
}


Hit findHit(const Scene &scene, const Ray &ray){
    Hit hitted;
    hitted.hit_happened = false;
    hitted.t = numeric_limits<float>::max();

    for(int sphere_number = 0; sphere_number < scene.spheres.size(); sphere_number++){
        Sphere current_sphere = scene.spheres[sphere_number];
        Vec3f center = scene.vertex_data[current_sphere.center_vertex_id - 1];
        float radius = current_sphere.radius;

        Hit hit = sphereIntersection(ray, center, radius, current_sphere.material_id, sphere_number);

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

    for(int triangle_number = 0; triangle_number < scene.triangles.size(); triangle_number++){
        Triangle current_triangle = scene.triangles[triangle_number];
        Vec3f v0 = scene.vertex_data[current_triangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[current_triangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[current_triangle.indices.v2_id - 1];

        Hit hit = triangleIntersection(ray, v0, v1, v2, current_triangle.material_id, triangle_number);
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

    for(int mesh_number = 0; mesh_number < scene.meshes.size(); mesh_number++){
        Mesh current_mesh = scene.meshes[mesh_number];
        Hit hit = meshIntersection(ray, current_mesh, scene, current_mesh.material_id, mesh_number);
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
    for(int sphereNumber = 0; sphereNumber < scene.spheres.size(); sphereNumber++){
        Sphere currentSphere = scene.spheres[sphereNumber];
        Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
        float radius = currentSphere.radius;

        Hit shadow_hit = sphereIntersection(shadow_ray, center, radius, currentSphere.material_id, sphereNumber);

        if(shadow_hit.hit_happened)
        {
            if(t > shadow_hit.t && shadow_hit.t >= 0)
            {
                return true;
            }
        }
    }

    for(int triangleNumber = 0; triangleNumber < scene.triangles.size(); triangleNumber++)
    {
        Triangle currentTriangle = scene.triangles[triangleNumber];
        Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];

        Hit shadow_hit = triangleIntersection(shadow_ray, v0, v1, v2, currentTriangle.material_id, triangleNumber);

        if(shadow_hit.hit_happened)
        {
            if(t > shadow_hit.t && shadow_hit.t >= 0)
            {
                return true;
            }
        }
    }

    for(int meshNumber = 0; meshNumber < scene.meshes.size(); meshNumber++)
    {
        Mesh currentMesh = scene.meshes[meshNumber];

        Hit shadow_hit = meshIntersection(shadow_ray, currentMesh, scene, currentMesh.material_id, meshNumber);

        if(shadow_hit.hit_happened)
        {
            if(t > shadow_hit.t && shadow_hit.t >= 0)
            {
                return true;
            }
        }
    }
}

Vec3f computeColor(const Scene &scene, Hit &hit, const Camera &camera, Ray &ray, int recursion_count){
    Vec3f color;
    float r = 0;
    float g = 0;
    float b = 0;

    if(hit.hit_happened){
        int material_id = hit.material_id;
        r = scene.materials[material_id - 1].ambient.x * scene.ambient_light.x;
        g = scene.materials[material_id - 1].ambient.y * scene.ambient_light.y;
        b = scene.materials[material_id - 1].ambient.z * scene.ambient_light.z;

        for(int light_number = 0; light_number < scene.point_lights.size(); light_number++){
            bool shadow_flag = false;

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
            shadow_ray.shadowRay = true;

            Hit shadow_hit;

            float t_light = subtract(current_light.position, shadow_ray.origin).x / shadow_ray.direction.x;
            shadow_flag = shadow(shadow_ray, scene, t_light);

            // if no shadowRay intersection but there is object-ray intersection
	        if(!shadow_flag || (shadow_flag && light_to_cam == 0))
	        {
		        int material_id = hit.material_id;

		        Vec3f diffuse = findDiffuse(current_light, scene, material_id, hit.normal, hit.intersection_point);
				
		        Vec3f specular = findSpecular(current_light, scene, ray, material_id, hit.normal, hit.intersection_point);
		                    
              	r += diffuse.x + specular.x;
  		        g += diffuse.y + specular.y;
  		        b += diffuse.z + specular.z;
            
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
            reflection_ray.shadowRay = false;

            Hit new_hit = findHit(scene, reflection_ray);
            if(!(new_hit.type == hit.type && new_hit.object_id == hit.object_id)){
                mirror_color = computeColor(scene, new_hit, camera, reflection_ray, recursion_count - 1);
                r += mirror_color.x * scene.materials[hit.material_id - 1].mirror.x;
                g += mirror_color.y * scene.materials[hit.material_id - 1].mirror.y;
                b += mirror_color.z * scene.materials[hit.material_id - 1].mirror.z;
            }
        }
    }else{
        r = scene.background_color.x;
        g = scene.background_color.y;
        b = scene.background_color.z;
    }

    color.x = r;
    color.y = g;
    color.z = b;

    return color;
}

int main(int argc, char* argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    
    int number_of_cameras = scene.cameras.size();

    for(int camera_no = 0; camera_no < number_of_cameras; camera_no++){
        int width = scene.cameras[camera_no].image_width;
        int height = scene.cameras[camera_no].image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                Ray ray = generateRay(scene.cameras[camera_no], i, j);
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

