#version 330

attribute vec3 position;
attribute vec3 color;
attribute vec3 normal;
attribute vec2 texCoords;

out vec2 texCoord0;
out vec3 normal0;
out vec3 color0;
out vec3 position0;

uniform mat4 Proj;
uniform mat4 View;
uniform mat4 Model;
uniform mat4 LookAt;

void main()
{	
	texCoord0 = texCoords;
	color0 = color;
	normal0 = vec3(LookAt * vec4(normal, 0.0));
	position0 = vec3(LookAt * vec4(mat3(Proj) * position, 1.0));
	gl_Position = Proj *View * Model* vec4(position, 1.0); //you must have gl_Position
}
