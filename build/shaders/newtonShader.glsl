#version 330

in vec2 texCoord0;
in vec3 normal0;
in vec3 color0;
in vec3 position0;

uniform int iterationNum;
uniform vec4 coeffs;
uniform vec4 rootA;
uniform vec4 rootB;
uniform vec4 rootC;
uniform vec4 colorA;
uniform vec4 colorB;
uniform vec4 colorC;
uniform float translateX;
uniform float translateY;
uniform float zoomFactor;
out vec4 Color;

#define complex_mul(a, b) vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
#define complex_div(a, b) vec2((a.x*b.x + a.y*b.y)/(b.x * b.x + b.y * b.y), (a.y * b.x - a.x * b.y)/(b.x*b.x + b.y*b.y));

vec2 getCoordinatesFromScreen(vec2 texCoord, vec4 coordinateRange){
	float horizontalSize = coordinateRange.y-coordinateRange.x; //Size of horizontal part of screen
	float verticalSize = coordinateRange.w-coordinateRange.z; //Size of vertical part of screen
	vec2 toRet;
	toRet.x = texCoord.x*horizontalSize+coordinateRange.x; //Multiply by size and add initial offset position
	toRet.y = texCoord.y*verticalSize+coordinateRange.z;

	return toRet;
}

void main()
{
	vec4 coordRange = vec4(-2.0,2.0,-2.0,2.0)*zoomFactor + vec4(vec2(translateX), vec2(-translateY));
	vec2 pos = getCoordinatesFromScreen(texCoord0, coordRange); //vec2(position0.x, position0.y);
	float epsilon = 1e-4; // TODO: times the length of the interval
	vec2 z = pos;
	int i;
	for (i = 0; i < iterationNum; i++) {
		vec2 quadraticZ = complex_mul(z, z);
		vec2 cubicZ = complex_mul(quadraticZ, z);
	    vec2 f = coeffs[0] * cubicZ + coeffs[1] * quadraticZ + coeffs[2] * z + vec2(coeffs[3], 0);
		vec2 fd = coeffs[0] * 3 * quadraticZ + coeffs[1] * 2 * z + vec2(coeffs[2], 0);
		vec2 divided = complex_div(f, fd);
		z = z - divided;
		if (length(divided) < epsilon) {
			break;
		}
	}
	vec4[3] colors = vec4[3](colorA, colorB, colorC);
	vec4[3] roots = vec4[3](rootA, rootB, rootC);
	int currentMinIndex = 0;
	for (int j = 1; j < 3; j++) {
		if(distance(z, vec2(roots[j].x, roots[j].y)) < distance(z, vec2(roots[currentMinIndex].x, roots[currentMinIndex].y))) {
			currentMinIndex = j;
		}
	}
	vec4 color = colors[currentMinIndex];

//  Intensitiy based on distance
//	float d = distance(pos, vec2(roots[currentMinIndex].x, roots[currentMinIndex].y));
//	float intensitity = exp(-d / 2.5);

	// Intesity based on iterations count
	float intensity = (1.0 * (iterationNum - i)) / iterationNum;

	intensity = clamp(intensity, 0.25, 1);
	color = intensity * color;
	//vec4 color = texture(sampler1, texCoord0);
	Color = color; //you must have gl_FragColor
}
