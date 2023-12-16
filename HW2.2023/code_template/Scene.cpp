#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
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
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
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

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
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
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	// command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	// system(command.c_str());
	
	command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	meshesTransformations();
	Matrix4 cameraMatrix = cameraTransformation(camera);
	Matrix4 viewPortMatrix = getViewPortTransformation(camera->horRes, camera->verRes);

	for(int meshNumber = 0; meshNumber < this->meshes.size(); meshNumber++)
	{
		vector<Vec3> tmp;
		transformedVertices.push_back(tmp);

		transformVertices(cameraMatrix, viewPortMatrix, meshNumber);
	}
	drawMeshes(camera, viewPortMatrix);

	clearTransformed();
}

// for every mesh it will find the new vertex position according to transformations and camera matrix
void Scene::transformVertices(Matrix4 cameraMatrix, Matrix4 viewPortMatrix, int meshNumber)
{
	Matrix4 result = multiplyMatrixWithMatrix(cameraMatrix, transformationsResult[meshNumber]);

	for(int i = 0; i < this->vertices.size(); i++)
	{
			
		Vec4 vertex;
		vertex.x = this->vertices[i]->x;
		vertex.y = this->vertices[i]->y;
		vertex.z = this->vertices[i]->z;
		vertex.t = 1;
		vertex.colorId = this->vertices[i]->colorId;

		vertex = multiplyMatrixWithVec4(result, vertex);

		Vec3 resultVertex;
		double w = vertex.t;

		resultVertex.x = vertex.x / w;
		resultVertex.y = vertex.y / w;
		resultVertex.z = vertex.z / w;
		resultVertex.colorId = vertex.colorId;

		transformedVertices[meshNumber].push_back(resultVertex);
	}
}

void Scene::transformations()
{
	for(int i = 0; i < this->translations.size(); i++)
	{
		this->translationsMatrix.push_back(getTranslationMatrix(*translations[i]));
	}

	for(int i = 0; i < this->scalings.size(); i++)
	{
		this->scalingsMatrix.push_back(getScalingMatrix(*scalings[i]));
	}

	for(int i = 0; i < this->rotations.size(); i++)
	{
		this->rotationsMatrix.push_back(getRotationMatrix(*rotations[i]));
	}
}

void Scene::meshesTransformations()
{
	transformations();
	for(int i = 0; i < this->meshes.size(); i++)
	{
		Matrix4 result = getIdentityMatrix();
		for(int j = 0; j < this->meshes[i]->numberOfTransformations; j++)
		{
			if(this->meshes[i]->transformationTypes[j] == 't')
			{
				result = multiplyMatrixWithMatrix(this->translationsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}
			else if(this->meshes[i]->transformationTypes[j] == 's')
			{
				result = multiplyMatrixWithMatrix(this->scalingsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}
			else if(this->meshes[i]->transformationTypes[j] == 'r')
			{
				result = multiplyMatrixWithMatrix(this->rotationsMatrix[this->meshes[i]->transformationIds[j] - 1], result);
			}
		}

		this->transformationsResult.push_back(result);
	}
}

float Scene::f01(int x, int y, int x0, int y0, int x1, int y1) 
{
	return x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1;
}

float Scene::f12(int x, int y, int x1, int y1, int x2, int y2) 
{
	return x * (y1 - y2) + y * (x2 - x1) + x1 * y2 - y1 * x2;
}

float Scene::f20(int x, int y, int x2, int y2, int x0, int y0) 
{
	return x * (y2 - y0) + y * (x0 - x2) + x2 * y0 - y2 * x0;
}

bool Scene::backfaceCulling(Vec3 inverseW, Triangle triangle, int meshNumber, Matrix4 viewportMatrix)
{
	Vec3 a = this->transformedVertices[meshNumber][triangle.vertexIds[0] - 1];
	Vec3 b = this->transformedVertices[meshNumber][triangle.vertexIds[1] - 1];
	Vec3 c = this->transformedVertices[meshNumber][triangle.vertexIds[2] - 1];
	Vec4 a4 = multiplyMatrixWithVec4(viewportMatrix, Vec4(a.x, a.y, a.z, 1, a.colorId));
	Vec4 b4 = multiplyMatrixWithVec4(viewportMatrix, Vec4(b.x, b.y, b.z, 1, b.colorId));
	Vec4 c4 = multiplyMatrixWithVec4(viewportMatrix, Vec4(c.x, c.y, c.z, 1, c.colorId));

	a = perspectiveDivide(a4);
	b = perspectiveDivide(b4);
	c = perspectiveDivide(c4);

	Vec3 n = findNormal(a, b, c);
	n = normalizeVec3(n);
	inverseW = normalizeVec3(inverseW);
	double dot = dotProductVec3(inverseW, n);

	return (dot > -10e-3);
}

void Scene::draw(int x, int y, double z, Color c) 
{
	if(this->depth[x][y] > z && z >= -1)
	{
		this->depth[x][y] = z;
		this->image[x][y] = c;
	}
}

void Scene::drawWireframe(int x, int y, Color c)
{
	this->image[x][y] = c;
}

void Scene::midpointWithInterpolation(Vec3 p0, Vec3 p1, Color c0, Color c1, Camera* camera)
{
	int x0 = p0.x;
	int y0 = p0.y;
	int x1 = p1.x;
	int y1 = p1.y;
	double z0 = p0.z;
	double z1 = p1.z;

	double m = (double)(y1 - y0) / ((double)(x1 - x0));

	if(m > 1.0)
	{
		if(y0 > y1)
		{
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;

			double tmpz;
			tmpz = z1;
			z1 = z0;
			z0 = tmp;
		}
		int x = x0;
		int d = 2 * (x0 - x1) + (y1 - y0);

		Color c = c0;
		double r, g, b;
		r = (c1.r - c0.r) / (y1 - y0);
		g = (c1.g - c0.g) / (y1 - y0);
		b = (c1.b - c0.b) / (y1 - y0);
		Color dc(r, g, b);

		double dz = (p1.z - p0.z) / (y1 - y0);

		for(int y = y0; y <= y1; y++)
		{
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));

			drawWireframe(x, y, c);
			z0 += dz;

			if(d < 0)
			{
				x += 1;
				d += 2 * ((x0 - x1) + (y1 - y0));
			}
			else 
			{
				d += 2 * (x0 - x1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

		}

	}else if(m <= 1.0 && m > 0.0)
	{
		if(x0 > x1)
		{
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;

			double tmpz;
			tmpz = z1;
			z1 = z0;
			z0 = tmpz;
		}
		int y = y0;
		int d = 2 * (y0 - y1) + (x1 - x0);

		Color c = c0;
		double r, g, b; // color change ratio for r g b
		r = (c1.r - c0.r) / (x1 - x0);
		g = (c1.g - c0.g) / (x1 - x0);
		b = (c1.b - c0.b) / (x1 - x0);
		Color dc(r, g, b);

		double dz = (p1.z - p0.z) / (x1 - x0);

		for(int x = x0; x <= x1; x++)
		{
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));

			drawWireframe(x, y, c);

			if(d < 0)
			{
				y += 1;
				d += 2 * ((y0 - y1) + (x1 - x0));
			} 
			else 
			{
				d += 2 * (y0 - y1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

			z0 += dz;
		}
	}
	else if(m <= 0.0 && m > -1.0)
	{
		if(x0 > x1)
		{
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;

			double tmpz;
			tmpz = z1;
			z1 = z0;
			z0 = tmp;
		}
		int y = y0;
		int d = 2 * (y1 - y0) + (x1 - x0);

		Color c = c0;
		double r, g, b; // color change ratio for r g b
		r = (c1.r - c0.r) / (x1 - x0);
		g = (c1.g - c0.g) / (x1 - x0);
		b = (c1.b - c0.b) / (x1 - x0);
		Color dc(r, g, b);

		double dz = (p1.z - p0.z) / (x1 - x0);

		for(int x = x0; x <= x1; x++)
		{
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));

			drawWireframe(x, y, c);
			if(d < 0)
			{
				y -= 1;
				d += 2 * ((y1 - y0) + (x1 - x0));
			} 
			else 
			{
				d += 2 * (y1 - y0);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

			z0 += dz;
		}
	}
	else if(m <= -1.0)
	{
		if(y0 < y1)
		{
			int tmp = x1;
			x1 = x0;
			x0 = tmp;
			tmp = y1;
			y1 = y0;
			y0 = tmp;
			Color colorTmp = c1;
			c1 = c0;
			c0 = colorTmp;

			double tmpz;
			tmpz = z1;
			z1 = z0;
			z0 = tmp;
		}
		int x = x0;
		int d = 2 * (x0 - x1) + (y0 - y1);

		Color c = c0;
		double r, g, b;
		r = (c1.r - c0.r) / (y0 - y1);
		g = (c1.g - c0.g) / (y0 - y1);
		b = (c1.b - c0.b) / (y0 - y1);
		Color dc(r, g, b);

		double dz = (p1.z - p0.z) / (y0 - y1);

		for(int y = y0; y >= y1; y--)
		{
			c.r = makeBetweenZeroAnd255(round(c.r));
			c.g = makeBetweenZeroAnd255(round(c.g));
			c.b = makeBetweenZeroAnd255(round(c.b));

			drawWireframe(x, y, c);

			if(d < 0)
			{
				x += 1;
				d += 2 * ((x0 - x1) + (y0 - y1));
			}
			else 
			{
				d += 2 * (x0 - x1);
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;

			z0 += dz;
		}
	}
}

void Scene::triangleRasterization(Vec3 p0, Vec3 p1, Vec3 p2, Color c0, Color c1, Color c2, Camera* camera) 
{
	int x0 = p0.x;
	int y0 = p0.y;
	int x1 = p1.x;
	int y1 = p1.y;
	int x2 = p2.x;
	int y2 = p2.y;

	double z0 = p0.z;
	double z1 = p1.z;
	double z2 = p2.z;



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
				
				double z = alpha * z0 + beta * z1 + gamma * z2;

				c.r = makeBetweenZeroAnd255(round(c.r));
				c.g = makeBetweenZeroAnd255(round(c.g));
				c.b = makeBetweenZeroAnd255(round(c.b));
				if(x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes)
				{
					draw(x, y, z, c);
				}
			}
		}
	}
}	

// for clipping it will return whether it is visible
bool Scene::isVisible(double den, double num, double &tEnter, double &tLeave)
{
	double t;
	if(den > 0)
	{ // potentially entering
		t = num / den;
		if(t > tLeave) return false;
		if(t > tEnter) tEnter = t;
	}
	else if(den < 0)
	{ // potentially leaving
		t = num / den;
		if(t < tEnter) return false;
		if(t < tLeave) tLeave = t;
	}
	else if(num > 0)
	{ // line parallel to edge
		return false;
	}
	return true;
}

// I will use Liang-Barsky Algorithm for lines clipping
std::vector<Vec3> Scene::clipping(Vec3& point1, Vec3& point2, bool& visible)
{
	double dx = point2.x - point1.x;
	double dy = point2.y - point1.y;
	double dz = point2.z - point1.z;
	double tE = 0;
	double tL = 1;
	visible = false;
	vector<Vec3> result;
	Vec3 p1 = point1;
	Vec3 p2 = point2;
	if(isVisible(dx, -1 - point1.x, tE, tL))
		if(isVisible(-dx, point1.x - 1, tE, tL))
			if(isVisible(dy, -1 - point1.y, tE, tL))
				if(isVisible(-dy, point1.y - 1, tE, tL))
					if(isVisible(dz, -1 - point1.z, tE, tL))
						if(isVisible(-dz, point1.z - 1, tE, tL)){
							visible = true;
							if(tL < 1){ // if it is visible change accordingly
								p2.x = p1.x + dx * tL;
								p2.y = p1.y + dy * tL;
								p2.z = p1.z + dz * tL;
							}
							if(tE > 0){
								p1.x = p1.x + dx * tE;
								p1.y = p1.y + dy * tE;
								p1.z = p1.z + dz * tE;
							}
						}
	result.push_back(p1);
	result.push_back(p2);
	return result;
}

void Scene::drawMeshes(Camera* camera, Matrix4 viewPortMatrix){
	for(int meshNumber = 0; meshNumber < this->meshes.size(); meshNumber++){
		for(int i = 0; i < this->meshes[meshNumber]->numberOfTriangles; i++){
			if(!cullingEnabled || (cullingEnabled && backfaceCulling(camera->w, this->meshes[meshNumber]->triangles[i],meshNumber, viewPortMatrix))){
				if(this->meshes[meshNumber]->type == 0){ // wireframe
					int id1 = this->meshes[meshNumber]->triangles[i].vertexIds[0];
					int id2 = this->meshes[meshNumber]->triangles[i].vertexIds[1];
					int id3 = this->meshes[meshNumber]->triangles[i].vertexIds[2];
					// for id1 and id2
					Vec3 p1, p2;
					Vec4 pv1, pv2;
					bool visible;
					visible = false;
					vector<Vec3> tmp = clipping(this->transformedVertices[meshNumber][id1 - 1], this->transformedVertices[meshNumber][id2 - 1], visible);
					if(visible)
					{
						p1 = tmp[0];
						p2 = tmp[1];

						pv1 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p1.x, p1.y, p1.z, 1, this->transformedVertices[meshNumber][id1 - 1].colorId));
						pv2 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p2.x, p2.y, p2.z, 1, this->transformedVertices[meshNumber][id2 - 1].colorId));
											
						p1 = perspectiveDivide(pv1);
						p2 = perspectiveDivide(pv2);						

						midpointWithInterpolation(p1, p2,
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1], camera);
					}
					

					// for id1 and id3 
					vector<Vec3> tmp2 = clipping(this->transformedVertices[meshNumber][id1 - 1], this->transformedVertices[meshNumber][id3 - 1], visible);
					if(visible)
					{
						p1 = tmp2[0];
						p2 = tmp2[1];

						pv1 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p1.x, p1.y, p1.z, 1, this->transformedVertices[meshNumber][id1 - 1].colorId));
						pv2 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p2.x, p2.y, p2.z, 1, this->transformedVertices[meshNumber][id3 - 1].colorId));
											
						p1 = perspectiveDivide(pv1);
						p2 = perspectiveDivide(pv2);

						midpointWithInterpolation(p1, p2,
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1], camera);
					}
					

					//for id2 and id3
					vector<Vec3> tmp3 = clipping(this->transformedVertices[meshNumber][id2 - 1], this->transformedVertices[meshNumber][id3 - 1], visible);
					if(visible)
					{
						p1 = tmp3[0];
						p2 = tmp3[1];

						pv1 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p1.x, p1.y, p1.z, 1, this->transformedVertices[meshNumber][id2 - 1].colorId));
						pv2 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p2.x, p2.y, p2.z, 1, this->transformedVertices[meshNumber][id3 - 1].colorId));
											
						p1 = perspectiveDivide(pv1);
						p2 = perspectiveDivide(pv2);
						
						midpointWithInterpolation(p1, p2,
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1],
						*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1], camera);
					}			
				}
				else
				{ // solid
					int id1 = this->meshes[meshNumber]->triangles[i].vertexIds[0];
					int id2 = this->meshes[meshNumber]->triangles[i].vertexIds[1];
					int id3 = this->meshes[meshNumber]->triangles[i].vertexIds[2];
					
					Vec3 p1, p2, p3;
					p1 = this->transformedVertices[meshNumber][id1 - 1];
					p2 = this->transformedVertices[meshNumber][id2 - 1];
					p3 = this->transformedVertices[meshNumber][id3 - 1];

					Vec4 pv1, pv2, pv3;
					pv1 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p1.x, p1.y, p1.z, 1, this->transformedVertices[meshNumber][id1 - 1].colorId));
					pv2 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p2.x, p2.y, p2.z, 1, this->transformedVertices[meshNumber][id2 - 1].colorId));
					pv3 = multiplyMatrixWithVec4(viewPortMatrix, Vec4(p3.x, p3.y, p3.z, 1, this->transformedVertices[meshNumber][id3 - 1].colorId));

					p1 = perspectiveDivide(pv1);
					p2 = perspectiveDivide(pv2);
					p3 = perspectiveDivide(pv3);

					triangleRasterization(p1, p2, p3,
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id1 - 1].colorId - 1],
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id2 - 1].colorId - 1],
					*this->colorsOfVertices[this->transformedVertices[meshNumber][id3 - 1].colorId - 1], camera);
				}
			}
		}
	}
}

Matrix4 Scene::cameraTransformation(Camera* camera)
{
	Matrix4 result;
	if(camera->projectionType == 0)
	{
		result = getOrthographicProjectionMatrix(camera);
	}
	else if(camera->projectionType == 1)
	{
		result = getPerspectiveTransformationMatrix(camera);
	}

	result = multiplyMatrixWithMatrix(result, getCameraTransformMatrix(camera->position, camera->u, camera->v, camera->w));

	return result;
}

void Scene::clearTransformed()
{
	this->transformedVertices.clear();
	this->transformationsResult.clear();
}



