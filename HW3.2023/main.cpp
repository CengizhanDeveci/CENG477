#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <GL/glew.h>    // The GL Header File
#include <GL/gl.h>      // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <ft2build.h>
#include FT_FREETYPE_H
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

using namespace std;

GLuint gProgram[4];
GLint gIntensityLoc;
float gIntensity = 1000;
int gWidth = 800, gHeight = 800;
int score = 0;
float speed = 2.0f;

bool objects[3];

float xPosition = 0.0f;

struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) {}
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) {}
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) {}
    GLfloat x, y, z;
};

struct Face
{
    Face(int v[], int t[], int n[])
    {
        vIndex[0] = v[0];
        vIndex[1] = v[1];
        vIndex[2] = v[2];
        tIndex[0] = t[0];
        tIndex[1] = t[1];
        tIndex[2] = t[2];
        nIndex[0] = n[0];
        nIndex[1] = n[1];
        nIndex[2] = n[2];
    }
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
int gVertexDataByIndex[4];
int gNormalDataByIndex[4];
int faceSize[3];

bool getHit = false;
int moveSide = 0;
bool endGame = false;
double currentTime;
bool restarted = false;

/// Holds all state information relevant to a character as loaded using FreeType
struct Character
{
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing; // Offset from baseline to left/top of glyph
    GLuint Advance;     // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

void restart();

void randomize()
{
    int randomNo = rand() % 3;
    if (randomNo == 0)
    {
        objects[0] = true;
        objects[1] = false;
        objects[2] = false;
    }
    else if (randomNo == 1)
    {
        objects[0] = false;
        objects[1] = true;
        objects[2] = false;
    }
    else if (randomNo == 2)
    {
        objects[0] = false;
        objects[1] = false;
        objects[2] = true;
    }
}

bool ParseObj(const string &fileName)
{
    fstream myfile;
    static int i = 1;

    // Open the input
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;
        gVertexDataByIndex[i] = gVertexDataByIndex[i - 1];
        gNormalDataByIndex[i] = gNormalDataByIndex[i - 1];

        faceSize[i - 1] = 0;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                        gNormalDataByIndex[i] += 3 * sizeof(GLfloat);
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                        gVertexDataByIndex[i] += 3 * sizeof(GLfloat);
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
                    char c;
                    int vIndex[3], nIndex[3], tIndex[3];
                    str >> vIndex[0];
                    str >> c >> c; // consume "//"
                    str >> nIndex[0];
                    str >> vIndex[1];
                    str >> c >> c; // consume "//"
                    str >> nIndex[1];
                    str >> vIndex[2];
                    str >> c >> c; // consume "//"
                    str >> nIndex[2];

                    assert(vIndex[0] == nIndex[0] &&
                           vIndex[1] == nIndex[1] &&
                           vIndex[2] == nIndex[2]); // a limitation for now

                    // make indices start from 0
                    for (int c = 0; c < 3; ++c)
                    {
                        vIndex[c] -= 1;
                        nIndex[c] -= 1;
                        tIndex[c] -= 1;
                    }

                    faceSize[i - 1]++;
                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            // data += curLine;
            if (!myfile.eof())
            {
                // data += "\n";
            }
        }

        i++;

        myfile.close();
    }
    else
    {
        return false;
    }

    assert(gVertices.size() == gNormals.size());

    return true;
}

bool ReadDataFromFile(
    const string &fileName, ///< [in]  Name of the shader file
    string &data)           ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

void createVS(GLuint &program, const string &filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar *shader = (const GLchar *)shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(program, vs);
}

void createFS(GLuint &program, const string &filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar *shader = (const GLchar *)shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(program, fs);
}

void initShaders()
{
    gProgram[0] = glCreateProgram();
    gProgram[1] = glCreateProgram();
    gProgram[2] = glCreateProgram();
    gProgram[3] = glCreateProgram();

    createVS(gProgram[0], "vert0.glsl");
    createFS(gProgram[0], "frag0.glsl");

    createVS(gProgram[1], "vert1.glsl");
    createFS(gProgram[1], "frag1.glsl");

    createVS(gProgram[2], "vert2.glsl");
    createFS(gProgram[2], "frag2.glsl");

    createVS(gProgram[3], "vert_text.glsl");
    createFS(gProgram[3], "frag_text.glsl");

    glBindAttribLocation(gProgram[0], 0, "inVertex");
    glBindAttribLocation(gProgram[0], 1, "inNormal");
    glBindAttribLocation(gProgram[1], 0, "inVertex");
    glBindAttribLocation(gProgram[1], 1, "inNormal");
    glBindAttribLocation(gProgram[2], 0, "inVertex");
    glBindAttribLocation(gProgram[2], 1, "inNormal");
    glBindAttribLocation(gProgram[3], 2, "vertex");

    glLinkProgram(gProgram[0]);
    glLinkProgram(gProgram[1]);
    glLinkProgram(gProgram[2]);
    glLinkProgram(gProgram[3]);
    glUseProgram(gProgram[0]);

    gIntensityLoc = glGetUniformLocation(gProgram[0], "intensity");
    glUniform1f(gIntensityLoc, gIntensity);
}

void initVBO()
{
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &gVertexAttribBuffer);
    glGenBuffers(1, &gIndexBuffer);

    assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
    gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
    GLfloat *vertexData = new GLfloat[gVertices.size() * 3];
    GLfloat *normalData = new GLfloat[gNormals.size() * 3];
    GLuint *indexData = new GLuint[gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

    for (int i = 0; i < gVertices.size(); ++i)
    {
        vertexData[3 * i] = gVertices[i].x;
        vertexData[3 * i + 1] = gVertices[i].y;
        vertexData[3 * i + 2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
    }

    for (int i = 0; i < gNormals.size(); ++i)
    {
        normalData[3 * i] = gNormals[i].x;
        normalData[3 * i + 1] = gNormals[i].y;
        normalData[3 * i + 2] = gNormals[i].z;
    }

    for (int i = 0; i < gFaces.size(); ++i)
    {
        indexData[3 * i] = gFaces[i].vIndex[0];
        indexData[3 * i + 1] = gFaces[i].vIndex[1];
        indexData[3 * i + 2] = gFaces[i].vIndex[2];
    }

    glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // done copying; can free now
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    // glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[3]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/liberation/LiberationSerif-Italic.ttf", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer);
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x};
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void init()
{
    gVertexDataByIndex[0] = 0;
    gNormalDataByIndex[0] = 0;

    //

    ParseObj("quad.obj");

    ParseObj("bunny.obj");

    ParseObj("cube.obj");

    glEnable(GL_DEPTH_TEST);
    initShaders();
    initFonts(gWidth, gHeight);
    initVBO();
}

void drawModel(int i = 0)
{
    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataByIndex[i]));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes + gNormalDataByIndex[i]));

    // glDrawElements(GL_TRIANGLES, faceSize[i] * 3, GL_UNSIGNED_INT, 0);
    if (i == 2)
    {
        glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0); // i dont know why but it works
    }
    else
    {
        glDrawElements(GL_TRIANGLES, faceSize[i] * 3, GL_UNSIGNED_INT, 0);
    }
}

void renderText(const std::string &text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    float static height = y;
    // Activate corresponding render state
    glUseProgram(gProgram[3]);
    glUniform3f(glGetUniformLocation(gProgram[3], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++)
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = height - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            {xpos, ypos + h, 0.0, 0.0},
            {xpos, ypos, 0.0, 1.0},
            {xpos + w, ypos, 1.0, 1.0},

            {xpos, ypos + h, 0.0, 0.0},
            {xpos + w, ypos, 1.0, 1.0},
            {xpos + w, ypos + h, 1.0, 0.0}};

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        // glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}

void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    static float jumpY = 0;
    static float offset = -9.0f;
    static float length = -100.0f;
    static float rotate = 0.0f;

    // quad


    if(restarted)
    {
        jumpY = 0;
        offset = -9.0f;
        length = -100.0f;
        rotate = 0.0f;
        restarted = false;
    }

    glUseProgram(gProgram[1]);
    glLoadIdentity();

    glm::mat4 T = glm::translate(glm::mat4(1.f), glm::vec3(0.f, -8.f, -11.f + offset));
    glm::mat4 R = glm::rotate(glm::mat4(1.f), glm::radians(90.0f), glm::vec3(1, 0, 0));
    glm::mat4 S = glm::scale(glm::mat4(1.f), glm::vec3(20.f, 20.f, 1000.0f));
    glm::mat4 modelMat = T * S * R;
    glm::mat4 modelMatInv = glm::transpose(glm::inverse(modelMat));
    glm::mat4 perspMat = glm::perspective(glm::radians(90.0f), 1.f, 0.1f, 200.0f);

    GLuint offsetLocation = glGetUniformLocation(gProgram[1], "offset");
    glUniform1f(offsetLocation, offset);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    drawModel(0);

    assert(glGetError() == GL_NO_ERROR);

    // bunny

    glUseProgram(gProgram[0]);
    glLoadIdentity();

    T = glm::translate(glm::mat4(1.f), glm::vec3(xPosition, 0.f - 5.0f + jumpY, -10.f));
    if (endGame)
    {
        R = glm::rotate(glm::mat4(1.f), glm::radians(-90.0f), glm::vec3(0, 0, 1)) * glm::rotate(glm::mat4(1.f), glm::radians(-90.0f), glm::vec3(0, 1, 0));
    }
    else
    {
        R = glm::rotate(glm::mat4(1.f), glm::radians(-90.0f + rotate), glm::vec3(0, 1, 0));
    }
    modelMat = T * R;
    modelMatInv = glm::transpose(glm::inverse(modelMat));
    perspMat = glm::perspective(glm::radians(90.0f), 1.f, 0.1f, 200.0f);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    drawModel(1);

    // cubes

    glUseProgram(gProgram[2]);
    glLoadIdentity();

    T = glm::translate(glm::mat4(1.f), glm::vec3(-6.5f, -2.f, length));
    S = glm::scale(glm::mat4(1.f), glm::vec3(1.f, 1.5f, 1.0f));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    GLuint collectable = glGetUniformLocation(gProgram[2], "collectable");
    float temp;

    if (objects[0])
        temp = 0.0f;
    else
        temp = 1.0f;

    glUniform1f(collectable, temp);

    drawModel(2);

    glUseProgram(gProgram[2]);
    glLoadIdentity();

    T = glm::translate(glm::mat4(1.f), glm::vec3(0.f, -2.f, length));
    S = glm::scale(glm::mat4(1.f), glm::vec3(1.f, 1.5f, 1.0f));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    collectable = glGetUniformLocation(gProgram[2], "collectable");
    if (objects[1])
        temp = 0.0f;
    else
        temp = 1.0f;

    glUniform1f(collectable, temp);

    drawModel(2);

    glUseProgram(gProgram[2]);
    glLoadIdentity();

    T = glm::translate(glm::mat4(1.f), glm::vec3(6.5f, -2.f, length));
    S = glm::scale(glm::mat4(1.f), glm::vec3(1.f, 1.5f, 1.0f));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    collectable = glGetUniformLocation(gProgram[2], "collectable");

    if (objects[2])
        temp = 0.0f;
    else
        temp = 1.0f;

    glUniform1f(collectable, temp);

    drawModel(2);

    assert(glGetError() == GL_NO_ERROR);

    renderText("Score: " + to_string(score), 0, gHeight - 50, 1, glm::vec3(0, 1, 1));

    assert(glGetError() == GL_NO_ERROR);

    float time = glfwGetTime() * speed * 5.0f; // Assuming you use GLFW for time
    jumpY = 1.5f * sin(time * 2.0f);           // Adjust the frequency to control the jump speed

    offset += speed / 10;
    offset = offset > 11.0f ? -9.0f : offset;
    speed += 0.001;

    length += speed / 3;

    if (endGame)
    {
        rotate = 90.0f;
        moveSide = 0;
        speed = 0.0;
        jumpY = 0;
        double restartTime = glfwGetTime();
        if (restartTime - currentTime > 2.0f)
        {
            jumpY = 0;
            offset = -9.0f;
            length = -100.0f;
            rotate = 0.0f;

            restart();
        }
    }
    else
    {
        if (length > -5.0f) // align with the cube
        {
            if (xPosition < -5.0f)
            {
                if (objects[0])
                {
                    getHit = true;
                    score += 1000;
                }
                else
                {
                    currentTime = glfwGetTime();
                    endGame = true;
                }
            }
            else if (xPosition > -1.5f && xPosition < 1.5f)
            {
                if (objects[1])
                {
                    getHit = true;
                    score += 1000;
                }
                else
                {
                    currentTime = glfwGetTime();
                    endGame = true;
                }
            }
            else if (xPosition > 5.0f)
            {
                if (objects[2])
                {
                    getHit = true;
                    score += 1000;
                }
                else
                {
                    currentTime = glfwGetTime();
                    endGame = true;
                }
            }
            length = -100.0f;
            randomize();
        }

        if (getHit)
        {
            rotate += 3.f * speed;
            if (rotate > 360.f)
            {
                rotate = 0.f;
                getHit = false;
            }
        }
        score += 1 * speed;

        xPosition += moveSide * 0.075f * speed;
        if (xPosition > 9.f)
        {
            xPosition = 9.f;
        }
        else if (xPosition < -9.f)
        {
            xPosition = -9.f;
        }
    }
}

void reshape(GLFWwindow *window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);
}

void restart()
{
    getHit = false;
    moveSide = 0;
    score = 0;
    speed = 2.0f;
    xPosition = 0.0f;
    endGame = false;
    restarted = true;
    randomize();
}

void keyboard(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }

    if (key == GLFW_KEY_R && (action == GLFW_PRESS  || action == GLFW_REPEAT)) // RESTART
    {
        restart();
    }

    if (key == GLFW_KEY_A && (action == GLFW_PRESS)) // LEFT
    {
        if (action == GLFW_PRESS || action == GLFW_REPEAT)
        {
            moveSide = -1;
        }
        else
        {
            moveSide = 0;
        }
    }
    else if (key == GLFW_KEY_D) // RIGHT
    {
        if (action == GLFW_PRESS || action == GLFW_REPEAT)
        {
            moveSide = 1;
        }
        else
        {
            moveSide = 0;
        }
    }
    else
    {
        moveSide = 0;
    }
}

void mainLoop(GLFWwindow *window)
{
    randomize();
    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        display();
        glfwSwapBuffers(window);
    }
}

int main(int argc, char **argv) // Create Main Function For Bringing It All Together
{
    GLFWwindow *window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(gWidth, gHeight, "Bunny Run", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    glfwSetWindowTitle(window, "Bunny Run");

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, gWidth, gHeight); // need to call this once ourselves
    mainLoop(window);                 // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
