#version 330 core

layout(location = 0) in vec3 inVertex;
layout(location = 1) in vec3 inNormal;

out vec4 fragPos;
out vec3 N;
out vec3 FragVertex; // Add this line to declare an output variable for inVertex

uniform mat4 modelingMat;
uniform mat4 modelingMatInvTr;
uniform mat4 perspectiveMat;


void main(void)
{
    vec4 p = modelingMat * vec4(inVertex, 1); // translate to world coordinates
    vec3 Nw = vec3(modelingMatInvTr * vec4(inNormal, 0)); // provided by the programmer

    N = normalize(Nw);
    FragVertex = inVertex; // Pass inVertex to FragVertex
    fragPos = p;

    gl_Position = perspectiveMat * modelingMat * vec4(inVertex, 1);
}
