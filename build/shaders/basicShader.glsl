#version 330

in vec2 texCoord0;
in vec3 normal0;
in vec3 color0;
in vec3 position0;

uniform vec4 lightColor;
uniform sampler2D sampler1;
uniform vec4 lightDirection;

out vec4 Color;


bool inCircle(vec2 circleCenter, vec2 pos, float radius){
    return distance(circleCenter, pos) <= radius;
}

void main()
{
    vec3 color = vec3(0.0, 0.0, 0.0);
    float minCorner = 0.05;
    float maxCorner = 0.95;
    if(inCircle(vec2(0.5, 0.5), texCoord0, 0.6) && texCoord0.x < maxCorner && texCoord0.x > minCorner && texCoord0.y < maxCorner && texCoord0.y > minCorner){
        color = abs(normal0) + max(normal0.zxy, 0.0);
    }
    Color = vec4(color, 1.0);
}
