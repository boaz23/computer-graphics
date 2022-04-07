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

out vec4 Color;

#define complex_mul(a, b) vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
#define complex_div(a, b) vec2((a.x*b.x + a.y*b.y)/(b.x * b.x + b.y * b.y), (a.y * b.x - a.x * b.y)/(b.x*b.x + b.y*b.y));


vec2 fx(vec2 z)
{
	vec2 cubicZ = vec2(pow(z.x,3.0)-3.0*z.x*pow(z.y,2.0), 3.0*pow(z.x,2.0)*z.y-pow(z.y,3.0));
	vec2 quadraticZ = complex_mul(z, z);
    vec2 xn=vec2(coeffs[0] * cubicZ.x + coeffs[1] * quadraticZ.x + coeffs[2] * z.x + coeffs[3], coeffs[0] * cubicZ.y + coeffs[1] * quadraticZ.y + coeffs[2] * z.y);
    return xn;
}
vec2 fpx(vec2 z)
{
	vec2 quadraticZ = complex_mul(z, z);
    vec2 xn=vec2(coeffs[0] * 3 * quadraticZ.x + coeffs[1] * 2 * z.x + coeffs[2], coeffs[0] * 3 * quadraticZ.y + coeffs[1] * 2 * z.y);
    return xn;
}
void main()
{
	vec2 z = vec2(position0.x, position0.y);
	for (int i = 0; i < iterationNum; i++) {
		vec2 f = fx(z);
		vec2 fd = fpx(z);
		vec2 divided = complex_div(f, fd);
		z = z - divided;
	}
	vec4[3] colors = {colorA, colorB, colorC};
	vec4[3] roots = {rootA, rootB, rootC};
	int currentMinIndex = 0;
	for (int j = 1; j < 3; j++) {
		if(distance(z, vec2(roots[j].x, roots[j].y)) < distance(z, vec2(roots[currentMinIndex].x, roots[currentMinIndex].y))) {
			currentMinIndex = j;
		}
	}
	vec4 color = colors[currentMinIndex];
	//vec4 color = texture(sampler1, texCoord0);
	Color = color; //you must have gl_FragColor
}
