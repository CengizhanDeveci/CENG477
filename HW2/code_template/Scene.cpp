#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

// it takes the amount of translation and returns translation matrix
Matrix4 translation(Translation t){
	Matrix4 result;

	result.val[0][0] = 1;
	result.val[1][1] = 1;
	result.val[2][2] = 1;
	result.val[3][3] = 1;

	result.val[0][3] = t.tx;
	result.val[1][3] = t.ty;
	result.val[2][3] = t.tz;

	return result;
}

// it takes an angle as an argument and return the rotation matrix
Matrix4 rotation(Rotation r){
	Matrix4 result;

	// initializing u
	Vec3 u;
	u.x = r.ux;
	u.y = r.uy;
	u.z = r.uz;
	normalizeVec3(u);

	// finding v
	// firstly we need to find smallest absolute value
	double smallest = abs(u.x);
	char smallestIndex = 'x';
	if(smallest > abs(u.y)){
		smallest = abs(u.y);
		smallestIndex = 'y';
	}
	if(smallest > abs(u.z)){
		smallest = abs(u.z);
		smallestIndex = 'z';
	}

	
	// swapping the other two by making one of them negating
	Vec3 v;
	if(smallestIndex == 'x'){
		v.x = 0;
		v.y = -u.z;
		v.z = u.y;
	}else if(smallestIndex == 'y'){
		v.x = -u.z;
		v.y = 0;
		v.z = u.x;
	}else{
		v.x = -u.y;
		v.y = u.x;
		v.z = 0;
	}

	// finding the w by making cross product
	Vec3 w;
	w = crossProductVec3(u, v);
	normalizeVec3(w);
	normalizeVec3(v);
	normalizeVec3(u);

	Matrix4 matrixInverse;
	matrixInverse.val[0][0] = u.x;
	matrixInverse.val[0][1] = v.x;
	matrixInverse.val[0][2] = w.x;
	matrixInverse.val[1][0] = u.y;
	matrixInverse.val[1][1] = v.y;
	matrixInverse.val[1][2] = w.y;
	matrixInverse.val[2][0] = u.z;
	matrixInverse.val[2][1] = v.z;
	matrixInverse.val[2][2] = w.z;
	matrixInverse.val[3][3] = 1;

	Matrix4 matrixM;
	matrixM.val[0][0] = u.x;
	matrixM.val[0][1] = u.y;
	matrixM.val[0][2] = u.z;
	matrixM.val[1][0] = v.x;
	matrixM.val[1][1] = v.y;
	matrixM.val[1][2] = v.z;
	matrixM.val[2][0] = w.x;
	matrixM.val[2][1] = w.y;
	matrixM.val[2][2] = w.z;
	matrixM.val[3][3] = 1;

	// Matrix4 translationMatrix = translation(Translation(-1, center.x, center.y, center.z));

	// Matrix4 inverseTranslationMatrix = translation(Translation(-1, -center.x, -center.y, -center.z));
	
	// converts to radian from degree
	double radian = r.angle * 3.14 / 180;
	double cosTheta = cos(radian);
	double sinTheta = sin(radian);

	// we will rotate from the axis x now
	Matrix4 rotationMatrix;
	rotationMatrix.val[0][0] = 1;
	rotationMatrix.val[1][1] = cosTheta;
	rotationMatrix.val[1][2] = -sinTheta;
	rotationMatrix.val[2][1] = sinTheta;
	rotationMatrix.val[2][2] = cosTheta;
	rotationMatrix.val[3][3] = 1;
	
	result = multiplyMatrixWithMatrix(matrixInverse, multiplyMatrixWithMatrix(rotationMatrix, matrixM));

	return result;
}

Matrix4 scaling(Scaling scale){
	Matrix4 result;

	result.val[0][0] = scale.sx;
	result.val[1][1] = scale.sy;
	result.val[2][2] = scale.sz;
	result.val[3][3] = 1;

	return result;
}


Matrix4 cameraTransform(Vec3 center, Vec3 u, Vec3 v, Vec3 w){
	Matrix4 result;

	result.val[0][0] = u.x;
	result.val[0][1] = u.y;
	result.val[0][2] = u.z;
	result.val[1][0] = v.x;
	result.val[1][1] = v.y;
	result.val[1][2] = v.z;
	result.val[2][0] = w.x;
	result.val[2][1] = w.y;
	result.val[2][2] = w.z;
	result.val[3][3] = 1;

	result = multiplyMatrixWithMatrix(result, translation(Translation(-1, -center.x, -center.y, -center.z)));

	return result;
}

Matrix4 orthographicProjection(Camera* camera){
	Matrix4 result;

	result.val[0][0] = 2 / (camera->right - camera->left);
	result.val[0][3] = -(camera->right + camera->left) / (camera->right - camera->left);
	result.val[1][1] = 2 / (camera->top - camera->bottom);
	result.val[1][3] = - (camera->top + camera->bottom) / (camera->top - camera->bottom);
	result.val[2][2] = -2 / (camera->far - camera->near);
	result.val[2][3] = -(camera->far + camera->near) / (camera->far - camera->near); 
	result.val[3][3] = 1;

	return result;
}

Matrix4 pers_p2o(Camera* camera){
	Matrix4 result;

	result.val[0][0] = camera->near;
	result.val[1][1] = camera->near;
	result.val[2][2] = camera->near + camera->far;
	result.val[2][3] = camera->near * camera->far;
	result.val[3][2] = -1;

	return result;
}


// finds the perspective projection matrix by multiplying orthogonal matrix by p2o matrix.
Matrix4 perspectiveProjection(Camera* camera){
	Matrix4 result;
	
	Matrix4 orthMatrix = orthographicProjection(camera);
	Matrix4 p2oMatrix = pers_p2o(camera);

	result = multiplyMatrixWithMatrix(orthMatrix, p2oMatrix);

	return result;
}

Matrix4 viewPortTransformation(double nx, double ny){
	Matrix4 result;

	result.val[0][0] = nx / 2;
	result.val[0][3] = (nx - 1) / 2;
	result.val[1][1] = ny / 2;
	result.val[1][3] = (ny - 1) / 2;
	result.val[2][2] = 0.5;
	result.val[2][3] = 0.5;
	result.val[3][3] = 1.0;

	return result;
}

// for clipping it will return whether it is visible
bool isVisible(double den, double num, double &tEnter, double &tLeave){
	double t;
	if(den > 0){ // potentially entering
		t = num / den;
		if(t > tLeave) return false;
		if(t > tEnter) tEnter = t;
	}else if(den < 0){ // potentially leaving
		t = num / den;
		if(t < tEnter) return false;
		if(t < tLeave) tLeave = t;
	}else if(num > 0){ // line parallel to edge
		return false;
	}
	return true;
}


// I will use Liang-Barsky Algorithm for lines clipping
void clipping(Vec3& point1, Vec3& point2){
	double dx = point2.x - point1.x;
	double dy = point2.y - point1.y;
	double dz = point2.z - point1.z;
	double tE = 0;
	double tL = 1;
	bool visible = false;
	if(isVisible(dx, -1 - point1.x, tE, tL))
		if(isVisible(-dx, point1.x - 1, tE, tL))
			if(isVisible(dy, -1 - point1.y, tE, tL))
				if(isVisible(-dy, point1.y - 1, tE, tL))
					if(isVisible(dz, -1 - point1.z, tE, tL))
						if(isVisible(-dz, point1.z - 1, tE, tL)){
							visible = true;
							if(tL < 1){ // if it is visible change accordingly
								point2.x = point1.x + dx * tL;
								point2.y = point1.y + dy * tL;
								point2.z = point1.z + dz * tL;
							}
							if(tE > 0){
								point1.x = point1.x + dx * tE;
								point1.y = point1.y + dy * tE;
								point1.z = point1.z + dz * tE;
							}
						}
}


// push all transformations matrixes
void Scene::transformations(){

	for(int i = 0; i < this->translations.size(); i++){
		this->translationsMatrix.push_back(translation(*translations[i]));
	}
	
	for(int i = 0; i < this->scalings.size(); i++){
		this->scalingsMatrix.push_back(scaling(*scalings[i]));
	}
	
	for(int i = 0; i < this->rotations.size(); i++){
		this->rotationsMatrix.push_back(rotation(*rotations[i]));
	}
}

 // loop for the transform all meshes
void Scene::meshesTransformations(){

	// calculate the all transformations given in scene
	transformations();
	// for all meshes it finds the required transformations matrix
	for(int i = 0; i < this->meshes.size(); i++){
		Matrix4 result = getIdentityMatrix();
		for(int j = 0; j < this->meshes[i]->numberOfTransformations; j++){
			if(this->meshes[i]->transformationTypes[j] == 't'){
				result = multiplyMatrixWithMatrix(this->translationsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}else if(this->meshes[i]->transformationTypes[j] == 's'){
				result = multiplyMatrixWithMatrix(this->scalingsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}else if(this->meshes[i]->transformationTypes[j] == 'r'){
				result = multiplyMatrixWithMatrix(this->rotationsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}
		}
		this->transformationsResult.push_back(result);
	}
}


// it will return the camera transformation matrix
Matrix4 Scene::cameraTransformation(Camera* camera){
	Matrix4 result;
	if(camera->projectionType == 0){
		result = orthographicProjection(camera);
	}else if(camera->projectionType == 1){
		result = perspectiveProjection(camera);
	}

	result = multiplyMatrixWithMatrix(result, cameraTransform(camera->pos, camera->u, camera->v, camera->w));

	return result;
}


// for every mesh it will find the new vertex position according to transformations and camera matrix
void Scene::transformVertices(Matrix4 cameraMatrix, Matrix4 viewPortMatrix, int meshNumber){
	Matrix4 result = multiplyMatrixWithMatrix(cameraMatrix, transformationsResult[meshNumber]);

	for(int i = 0; i < this->vertices.size(); i++){
			
		Vec4 vertex;
		vertex.x = this->vertices[i]->x;
		vertex.y = this->vertices[i]->y;
		vertex.z = this->vertices[i]->z;
		vertex.t = 1;
		vertex.colorId = this->vertices[i]->colorId;

		vertex = multiplyMatrixWithVec4(result, vertex);

		vertex = multiplyMatrixWithVec4(viewPortMatrix, vertex);

		Vec3 resultVertex;
		double w = vertex.t;

		resultVertex.x = vertex.x / w;
		resultVertex.y = vertex.y / w;
		resultVertex.z = vertex.z / w;
		resultVertex.colorId = vertex.colorId;

		transformedVertices[meshNumber].push_back(resultVertex);
	}
}


Vec3 findNormal(Vec3 a, Vec3 b, Vec3 c){
	Vec3 normal;

	normal = crossProductVec3(subtractVec3(c, b), subtractVec3(a, b));
	normalizeVec3(normal);

	return normal;
}


bool Scene::backfaceCulling(Vec3 inverseW, Triangle triangle, int meshNumber){
	Vec3 a = this->transformedVertices[meshNumber][triangle.getFirstVertexId() - 1];
	Vec3 b = this->transformedVertices[meshNumber][triangle.getSecondVertexId() - 1];
	Vec3 c = this->transformedVertices[meshNumber][triangle.getThirdVertexId() - 1];
	Vec3 n = findNormal(a, b, c);

	double dot = dotProductVec3(inverseW, n);
	return true ? dot > 0 : false;
}

void Scene::draw(int x, int y, Color c) {
	this->image[x][y] = c;
}

void Scene::midpointWithInterpolation(int x0, int y0, int x1, int y1, Color c0, Color c1){
	
	double m = (double)(y1 - y0) / ((double)(x1 - x0));

	if(m > 1.0){
		if(y0 > y1){
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;
		}
		int x = x0;
		int d = 2 * (x0 - x1) + (y1 - y0);

		Color c = c0;
		double r, g, b;
		r = (c1.r - c0.r) / (y1 - y0);
		g = (c1.g - c0.g) / (y1 - y0);
		b = (c1.b - c0.b) / (y1 - y0);
		Color dc(r, g, b);

		for(int y = y0; y <= y1; y++){
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));
			draw(x, y, c);

			if(d < 0){
				x += 1;
				d += 2 * ((x0 - x1) + (y1 - y0));
			}else {
				d += 2 * (x0 - x1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

		}

	}else if(m <= 1.0 && m > 0.0){
		if(x0 > x1){
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;
		}
		int y = y0;
		int d = 2 * (y0 - y1) + (x1 - x0);

		Color c = c0;
		double r, g, b; // color change ratio for r g b
		r = (c1.r - c0.r) / (x1 - x0);
		g = (c1.g - c0.g) / (x1 - x0);
		b = (c1.b - c0.b) / (x1 - x0);
		Color dc(r, g, b);

		for(int x = x0; x <= x1; x++){
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));
			draw(x, y, c);

			if(d < 0){
				y += 1;
				d += 2 * ((y0 - y1) + (x1 - x0));
			} else {
				d += 2 * (y0 - y1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;
		}
	}else if(m <= 0.0 && m > -1.0){
		if(x0 > x1){
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;
		}
		int y = y0;
		int d = 2 * (y1 - y0) + (x1 - x0);

		Color c = c0;
		double r, g, b; // color change ratio for r g b
		r = (c1.r - c0.r) / (x1 - x0);
		g = (c1.g - c0.g) / (x1 - x0);
		b = (c1.b - c0.b) / (x1 - x0);
		Color dc(r, g, b);

		for(int x = x0; x <= x1; x++){
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));
			draw(x, y, c);

			if(d < 0){
				y -= 1;
				d += 2 * ((y1 - y0) + (x1 - x0));
			} else {
				d += 2 * (y1 - y0);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;
		}
	}else if(m <= -1.0){
		if(y0 < y1){
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;
		}
		int x = x0;
		int d = 2 * (x0 - x1) + (y0 - y1);

		Color c = c0;
		double r, g, b;
		r = (c1.r - c0.r) / (y0 - y1);
		g = (c1.g - c0.g) / (y0 - y1);
		b = (c1.b - c0.b) / (y0 - y1);
		Color dc(r, g, b);

		for(int y = y0; y >= y1; y--){
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));
			draw(x, y, c);

			if(d < 0){
				x += 1;
				d += 2 * ((x0 - x1) + (y0 - y1));
			}else {
				d += 2 * (x0 - x1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

		}
	}
}

float Scene::f01(int x, int y, int x0, int y0, int x1, int y1) {
	return x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1;
}

float Scene::f12(int x, int y, int x1, int y1, int x2, int y2) {
	return x * (y1 - y2) + y * (x2 - x1) + x1 * y2 - y1 * x2;
}

float Scene::f20(int x, int y, int x2, int y2, int x0, int y0) {
	return x * (y2 - y0) + y * (x0 - x2) + x2 * y0 - y2 * x0;
}

void Scene::triangleRasterization(int x0, int y0, int x1, int y1, int x2, int y2, Color c0, Color c1, Color c2) {
	float alpha, beta, gamma;
	int ymin, xmin, ymax, xmax;
	Color c;
	ymin = min(min(y0, y1), y2);
	xmin = min(min(x0, x1), x2);
	ymax = max(max(y0, y1), y2);
	xmax = max(max(x0, x1), x2);
	for (int y = ymin; y <= ymax; y++)
	{
		for (int x = xmin; x <= xmax; x++)
		{
			alpha = f12(x, y, x1, y1, x2, y2) / f12(x0, y0, x1, y1, x2, y2);
			beta = f20(x, y, x2, y2, x0, y0) / f20(x1, y1, x2, y2, x0, y0);
			gamma = f01(x, y, x0, y0, x1, y1) / f01(x2, y2, x0, y0, x1, y1);
			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				c.r = alpha * c0.r + beta * c1.r + gamma * c2.r;
				c.g = alpha * c0.g + beta * c1.g + gamma * c2.g;
				c.b = alpha * c0.b + beta * c1.b + gamma * c2.b;
				
				c.r = round(c.r);
				c.g = round(c.g);
				c.b = round(c.b);
				
				draw(x, y, c);
			}
		}
	}
}

void Scene::drawMeshes(Camera* camera){
	for(int meshNumber = 0; meshNumber < this->meshes.size(); meshNumber++){
		for(int i = 0; i < this->meshes[meshNumber]->numberOfTriangles; i++){

			if(this->meshes[meshNumber]->type == 0){ // wireframe
					int id1 = this->meshes[meshNumber]->triangles[i].getFirstVertexId();
					int id2 = this->meshes[meshNumber]->triangles[i].getSecondVertexId();
					int id3 = this->meshes[meshNumber]->triangles[i].getThirdVertexId();
					// for id1 and id2
					midpointWithInterpolation(this->transformedVertices[meshNumber][id1-1].x,
					this->transformedVertices[meshNumber][id1 - 1].y,
					this->transformedVertices[meshNumber][id2 - 1].x,
					this->transformedVertices[meshNumber][id2 - 1].y,
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1]);

					// for id1 and id3 
					midpointWithInterpolation(this->transformedVertices[meshNumber][id1-1].x,
					this->transformedVertices[meshNumber][id1 - 1].y,
					this->transformedVertices[meshNumber][id3 - 1].x,
					this->transformedVertices[meshNumber][id3 - 1].y,
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1]);

					//for id2 and id3
					midpointWithInterpolation(this->transformedVertices[meshNumber][id2-1].x,
					this->transformedVertices[meshNumber][id2 - 1].y,
					this->transformedVertices[meshNumber][id3 - 1].x,
					this->transformedVertices[meshNumber][id3 - 1].y,
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1],
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1]);

				}
			else if(!cullingEnabled || (cullingEnabled && backfaceCulling(inverseVec3(camera->w), this->meshes[meshNumber]->triangles[i],meshNumber))){
				 // solid
				int id1 = this->meshes[meshNumber]->triangles[i].getFirstVertexId();
				int id2 = this->meshes[meshNumber]->triangles[i].getSecondVertexId();
				int id3 = this->meshes[meshNumber]->triangles[i].getThirdVertexId();
				triangleRasterization(this->transformedVertices[meshNumber][id1-1].x, 
				this->transformedVertices[meshNumber][id1 - 1].y,
				this->transformedVertices[meshNumber][id2 - 1].x,
				this->transformedVertices[meshNumber][id2 - 1].y,
				this->transformedVertices[meshNumber][id3 - 1].x,
				this->transformedVertices[meshNumber][id3 - 1].y,
				*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
				*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1],
				*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1]);
			}
		}
	}
}


void Scene::clearTransformed(){
	this->transformedVertices.clear();
	this->transformationsResult.clear();
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function.
	meshesTransformations();
	Matrix4 viewPortMatrix = viewPortTransformation(camera->horRes, camera->verRes);

	Matrix4 cameraMatrix = cameraTransformation(camera);

	for(int meshNumber = 0; meshNumber < this->meshes.size(); meshNumber++){

		vector<Vec3> tmp;
		transformedVertices.push_back(tmp);
	
		transformVertices(cameraMatrix, viewPortMatrix, meshNumber);

	}
	drawMeshes(camera);

	clearTransformed();
}

/*	
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();
		
		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);
			
			if (result != EOF) {
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}