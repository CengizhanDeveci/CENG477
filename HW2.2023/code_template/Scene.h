#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);

	void transformations();
	void meshesTransformations();
	void transformVertices(Matrix4 cameraMatrix, Matrix4 viewPortMatrix, int meshNumber);
	bool backfaceCulling(Vec3 inverseW, Triangle triangle, int meshNumber, Matrix4 viewportMatrix);
	void midpointWithInterpolation(Vec3 p0, Vec3 p1, Color c0, Color c1, Camera* camera);
	float f01(int x, int y, int x0, int y0, int x1, int y1);
	float f12(int x, int y, int x1, int y1, int x2, int y2);
	float f20(int x, int y, int x2, int y2, int x0, int y0);
	void triangleRasterization(Vec3 p0, Vec3 p1, Vec3 p2, Color c0, Color c1, Color c2, Camera* camera);
	void draw(int x, int y, double z, Color c);
	void drawWireframe(int x, int y, Color c);
	void drawMeshes(Camera* camera, Matrix4 viewPortMatrix);
	void clearTransformed();
	std::vector<Vec3> clipping(Vec3& point1, Vec3& point2, bool& visible);
	bool isVisible(double den, double num, double &tEnter, double &tLeave);
	Matrix4 cameraTransformation(Camera* camera);

	std::vector<Matrix4> translationsMatrix;
	std::vector<Matrix4> scalingsMatrix;
	std::vector<Matrix4> rotationsMatrix;
	std::vector<Matrix4> transformationsResult; // each index equals the transformations matrix result for that mesh
	std::vector<std::vector<Vec3>> transformedVertices;

	
};

#endif
