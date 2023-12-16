#include <iostream>
#include <cmath>
#include "Helpers.h"

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

/*
 * It takes the amount of translation and returns translation matrix.
 */
Matrix4 getTranslationMatrix(Translation t)
{
    Matrix4 result;

    result.values[0][0] = 1.0;
    result.values[1][1] = 1.0;
    result.values[2][2] = 1.0;
    result.values[3][3] = 1.0;
    result.values[0][3] = t.tx;
    result.values[1][3] = t.ty;
    result.values[2][3] = t.tz;


    return result;
}

/*
* It takes the amount of rotation and returns rotation matrix.
*/
Matrix4 getRotationMatrix(Rotation r)
{
    Matrix4 result;
    
    Vec3 u;
    u.x = r.ux;
    u.y = r.uy;
    u.z = r.uz;
    u = normalizeVec3(u);

    // finding v
    double smallest = abs(u.x);
    char smallestIndex = 'x';
    if (abs(u.y) < smallest)
    {
        smallest = abs(u.y);
        smallestIndex = 'y';
    }
    if (abs(u.z) < smallest)
    {
        smallest = abs(u.z);
        smallestIndex = 'z';
    }

    // make smallest 0 and swapping by negating one
    Vec3 v;
    if (smallestIndex == 'x')
    {
        v.x = 0.0;
        v.y = -u.z;
        v.z = u.y;
    }
    else if (smallestIndex == 'y')
    {
        v.x = -u.z;
        v.y = 0.0;
        v.z = u.x;
    }
    else
    {
        v.x = -u.y;
        v.y = u.x;
        v.z = 0.0;
    }

    Vec3 w = crossProductVec3(u, v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    Matrix4 alignMatrix;
    alignMatrix.values[0][0] = u.x;
    alignMatrix.values[0][1] = u.y;
    alignMatrix.values[0][2] = u.z;
    alignMatrix.values[0][3] = 0.0;
    alignMatrix.values[1][0] = v.x;
    alignMatrix.values[1][1] = v.y;
    alignMatrix.values[1][2] = v.z;
    alignMatrix.values[1][3] = 0.0;
    alignMatrix.values[2][0] = w.x;
    alignMatrix.values[2][1] = w.y;
    alignMatrix.values[2][2] = w.z;
    alignMatrix.values[2][3] = 0.0;
    alignMatrix.values[3][0] = 0.0;
    alignMatrix.values[3][1] = 0.0;
    alignMatrix.values[3][2] = 0.0;
    alignMatrix.values[3][3] = 1.0;

    Matrix4 alignMatrixInverse;
    alignMatrixInverse.values[0][0] = u.x;
    alignMatrixInverse.values[0][1] = v.x;
    alignMatrixInverse.values[0][2] = w.x;
    alignMatrixInverse.values[0][3] = 0.0;
    alignMatrixInverse.values[1][0] = u.y;
    alignMatrixInverse.values[1][1] = v.y;
    alignMatrixInverse.values[1][2] = w.y;
    alignMatrixInverse.values[1][3] = 0.0;
    alignMatrixInverse.values[2][0] = u.z;
    alignMatrixInverse.values[2][1] = v.z;
    alignMatrixInverse.values[2][2] = w.z;
    alignMatrixInverse.values[2][3] = 0.0;
    alignMatrixInverse.values[3][0] = 0.0;
    alignMatrixInverse.values[3][1] = 0.0;
    alignMatrixInverse.values[3][2] = 0.0;
    alignMatrixInverse.values[3][3] = 1.0;

    // convert angle to radian
    double angle = r.angle * 3.14 / 180.0;
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);

    Matrix4 rotationMatrix;
    rotationMatrix.values[0][0] = 1.0;
    rotationMatrix.values[1][1] = cosTheta;
    rotationMatrix.values[1][2] = -sinTheta;
    rotationMatrix.values[2][1] = sinTheta;
    rotationMatrix.values[2][2] = cosTheta;
    rotationMatrix.values[3][3] = 1.0;

    result = multiplyMatrixWithMatrix(alignMatrixInverse, rotationMatrix);
    result = multiplyMatrixWithMatrix(result, alignMatrix);

    return result;  
}

/* 
 * It takes the amount of scaling and returns scaling matrix.
 */
Matrix4 getScalingMatrix(Scaling s)
{
    Matrix4 result;

    result.values[0][0] = s.sx;
    result.values[1][1] = s.sy;
    result.values[2][2] = s.sz;
    result.values[3][3] = 1.0;

    return result;
}

Matrix4 getCameraTransformMatrix(Vec3 center, Vec3 u, Vec3 v, Vec3 w)
{

    Matrix4 result;

    result.values[0][0] = u.x;
    result.values[0][1] = u.y;
    result.values[0][2] = u.z;
    result.values[0][3] = 0.0;

    result.values[1][0] = v.x;
    result.values[1][1] = v.y;
    result.values[1][2] = v.z;
    result.values[1][3] = 0.0;

    result.values[2][0] = w.x;
    result.values[2][1] = w.y;
    result.values[2][2] = w.z;
    result.values[2][3] = 0.0;

    result.values[3][0] = 0.0;
    result.values[3][1] = 0.0;
    result.values[3][2] = 0.0;
    result.values[3][3] = 1.0;

    Matrix4 translationMatrix;
    translationMatrix.values[0][0] = 1.0;
    translationMatrix.values[1][1] = 1.0;
    translationMatrix.values[2][2] = 1.0;
    translationMatrix.values[3][3] = 1.0;
    translationMatrix.values[0][3] = -center.x;
    translationMatrix.values[1][3] = -center.y;
    translationMatrix.values[2][3] = -center.z;

    result = multiplyMatrixWithMatrix(result, translationMatrix);

    return result;
}

Matrix4 getOrthographicProjectionMatrix(Camera* camera)
{
    Matrix4 result;

    result.values[0][0] = 2.0 / (camera->right - camera->left);
    result.values[0][1] = 0.0;
    result.values[0][2] = 0.0;
    result.values[0][3] = -(camera->right + camera->left) / (camera->right - camera->left);

    result.values[1][0] = 0.0;
    result.values[1][1] = 2.0 / (camera->top - camera->bottom);
    result.values[1][2] = 0.0;
    result.values[1][3] = -(camera->top + camera->bottom) / (camera->top - camera->bottom);

    result.values[2][0] = 0.0;
    result.values[2][1] = 0.0;
    result.values[2][2] = -2.0 / (camera->far - camera->near);
    result.values[2][3] = -(camera->far + camera->near) / (camera->far - camera->near);

    result.values[3][0] = 0.0;
    result.values[3][1] = 0.0;
    result.values[3][2] = 0.0;
    result.values[3][3] = 1.0;

    return result;
}

Matrix4 getPersp2o(Camera* camera)
{
    Matrix4 result;

    result.values[0][0] = camera->near;
    result.values[1][1] = camera->near;
    result.values[2][2] = camera->near + camera->far;
    result.values[2][3] = camera->near * camera->far;
    result.values[3][2] = -1.0;

    return result;
}

Matrix4 getPerspectiveTransformationMatrix(Camera* camera)
{
    Matrix4 result;

    result = multiplyMatrixWithMatrix(getOrthographicProjectionMatrix(camera), getPersp2o(camera));

    return result;
}

Vec3 perspectiveDivide(Vec4 point)
{
    Vec3 result;

    result.x = point.x / point.t;
    result.y = point.y / point.t;
    result.z = point.z / point.t;

    return result;
}

Vec3 findNormal(Vec3 a, Vec3 b, Vec3 c)
{
	Vec3 normal;

	normal = crossProductVec3(subtractVec3(c, b), subtractVec3(a, b));
	normalizeVec3(normal);

	return normal;
}

Matrix4 getViewPortTransformation(double nx, double ny)
{
	Matrix4 result;

	result.values[0][0] = nx / 2;
	result.values[0][3] = (nx - 1) / 2;
	result.values[1][1] = ny / 2;
	result.values[1][3] = (ny - 1) / 2;
	result.values[2][2] = 0.5;
	result.values[2][3] = 0.5;
	result.values[3][3] = 1.0;

	return result;
}

