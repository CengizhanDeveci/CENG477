#version 330 core

in vec4 fragPos;
in vec3 N;
in vec3 FragVertex;

uniform vec3 lightPos;
uniform vec3 eyePos;

uniform float offset;

out vec4 fragColor;

vec3 Iamb = vec3(0.8, 0.8, 0.8);

vec3 kaWhite = vec3(1, 1, 1);
vec3 kaBlue = vec3(0.27, 0.51, 0.7);

void main(void)
{
    // Apply the boolean logic from your provided code

	bool x = (FragVertex.x > 0.5) || ((FragVertex.x < 0.0) && (FragVertex.x) > -0.5) ;

	float tempY = FragVertex.y * 100.0 - 6.0;
	while(tempY < -2.0)
	{
		tempY = tempY + 2.0;
	}
	bool y = (tempY) > -1.0;
	

    bool xorXY = x == y;

    // Use the boolean logic to determine the color
    vec3 ambientColor;
    if (xorXY) 
	{
        ambientColor = kaWhite * Iamb;
    } 
	else 
	{
        ambientColor = kaBlue * Iamb;
    }

    fragColor = vec4(ambientColor, 1.0);
}
