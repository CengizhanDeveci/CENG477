#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	void transformations();
	void meshesTransformations();

	vector<Matrix4> translationsMatrix;
	vector<Matrix4> scalingsMatrix;
	vector<Matrix4> rotationsMatrix;
	vector<Matrix4> transformationsResult; // each index equals the transformations matrix result for that mesh
	Matrix4 cameraTransformation(Camera* camera);
	vector<vector<Vec3>> transformedVertices;
	void transformVertices(Matrix4 cameraMatrix, Matrix4 viewPortMatrix, int meshNumber);
	bool backfaceCulling(Vec3 inverseW, Triangle triangle, int meshNumber);

	void midpointWithInterpolation(int x0, int y0, int x1, int y1, Color c0, Color c1);
	float f01(int x, int y, int x0, int y0, int x1, int y1);
	float f12(int x, int y, int x1, int y1, int x2, int y2);
	float f20(int x, int y, int x2, int y2, int x0, int y0);
	void triangleRasterization(int x0, int y0, int x1, int y1, int x2, int y2, Color c0, Color c1, Color c2);
	void draw(int x, int y, Color c);

};

#endif
