#include "Assignment2.h"
#include <iostream>


static void printMat(const Eigen::Matrix4d& mat)
{
	std::cout<<" matrix:"<<std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout<< mat(j,i)<<" ";
		std::cout<<std::endl;
	}
}

Assignment2::Assignment2()
{
}

//Assignment2::Assignment2(float angle ,float relationWH, float near, float far) : Scene(angle,relationWH,near,far)
//{ 	
//}

void Assignment2::Init()
{		
	SceneParser("scene.txt", &sceneData);
	AddShader("shaders/ourRayTracingShader");
    //AddShader("shaders/rayTracingShader-roee");
	AddShape(Plane, -1, TRIANGLES, 0);
    //AddShape(Plane, -1, TRIANGLES, 1);
	xOld = 0;
	yOld = 0;
	xRel = 0;
	yRel = 0;

	SetShapeShader(0, 0);
    //SetShapeShader(1, 1);
    pickedShape = 0;
	SetShapeStatic(0);
    //SetShapeStatic(1);
	upDownAngle = 0.0;
	leftRightAngle = 0.0;
	scale = 1.0;
	ComputeAngleFromEye();
}

void Assignment2::Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((shapeIndx+1) & 0x000000FF) >>  0;
	int g = ((shapeIndx+1) & 0x0000FF00) >>  8;
	int b = ((shapeIndx+1) & 0x00FF0000) >> 16;

	//uniform vec4 eye;
	//uniform vec4 ambient;
	//uniform vec4[20] objects;
	//uniform vec4[20] objColors;
	//uniform vec4[10] lightsDirection;
	//uniform vec4[10] lightsIntensity;
	//uniform vec4[10] lightsPosition;
	//uniform ivec4 sizes;
	s->Bind();
	s->SetUniform4f("eye", sceneData.eye(0), sceneData.eye(1), sceneData.eye(2), sceneData.eye(3));	
	s->SetUniform4f("ambient", sceneData.ambient(0), sceneData.ambient(1), sceneData.ambient(2), sceneData.ambient(3));
	s->SetUniform4fv("objects", sceneData.objects.data(), sceneData.objects.size());
	s->SetUniform4fv("objColors", sceneData.colors.data(), sceneData.colors.size());
	s->SetUniform4fv("lightsDirection", sceneData.directions.data(), sceneData.directions.size());
	s->SetUniform4fv("lightsIntensity", sceneData.intensities.data(), sceneData.intensities.size());
	s->SetUniform4fv("lightsPosition", sceneData.lights.data(), sceneData.lights.size());
	s->SetUniform4i("sizes", sceneData.sizes(0), sceneData.sizes(1), sceneData.sizes(2), sceneData.sizes(3));

	s->SetUniformMat4f("Proj", Proj);
	s->SetUniformMat4f("View", View);
	s->SetUniformMat4f("Model", Model);
	s->Unbind();
}


void Assignment2::WhenRotate()
{
}

float Assignment2::UpdatePosition(float xpos, float ypos)
{
	xRel = xOld - xpos;
	yRel = yOld - ypos;
	xOld = xpos;
	yOld = ypos;
	return yRel;
}

void Assignment2::WhenTranslate()
{
	float movCoeff = 2.0f;
	float dx = -xRel / movCoeff;
	float dy = yRel / movCoeff;
	TranslateX(dx);
	TranslateY(-dy);
	//std::cout << sceneData.eye << std::endl;
	ComputeAngleFromEye();
}

void Assignment2::TranslateX(float dx)
{
	sceneData.eye.x() += TRANSLATION_SENSITIVITY * dx * radius;
}

void Assignment2::TranslateY(float dy)
{
	sceneData.eye.y() += TRANSLATION_SENSITIVITY * dy * radius;
}

void Assignment2::ChangeZoomBy(float dz)
{
	float nextZoom = radius / (1 + SCALE_SENSITIVITY * dz);
	if (nextZoom >= 0) {
		radius = nextZoom;
		ComputeEyeFromAngle();
	}
}


void Assignment2::RotateScene(int unitsUp, int unitsRight) {
	upDownAngle += unitsUp * ANGLE_STEP;
	leftRightAngle += unitsRight * ANGLE_STEP;
	if (upDownAngle <= 0) {
		upDownAngle += 2*EIGEN_PI;
	}
	else if (upDownAngle >= 2*EIGEN_PI) {
		upDownAngle -= 2*EIGEN_PI;
	}
	if (leftRightAngle <= 0) {
		leftRightAngle += 2*EIGEN_PI;
	}
	else if (leftRightAngle >= 2*EIGEN_PI) {
		leftRightAngle -= 2*EIGEN_PI;
	}
	ComputeEyeFromAngle();
}

void Assignment2::Animate() {
    if(isActive)
	{
		
	}
}

void Assignment2::ScaleAllShapes(float amt,int viewportIndx)
{
	for (int i = 1; i < data_list.size(); i++)
	{
		if (data_list[i]->Is2Render(viewportIndx))
		{
            data_list[i]->MyScale(Eigen::Vector3d(amt, amt, amt));
		}
	}
}

Assignment2::~Assignment2(void)
{

}

void Assignment2::ComputeEyeFromAngle() {
	float x = radius * sin(upDownAngle) * sin(leftRightAngle);
	float y = radius * cos(upDownAngle);
	float z = radius * sin(upDownAngle) * cos(leftRightAngle);
	sceneData.eye.x() = x;
	sceneData.eye.y() = y;
	sceneData.eye.z() = z;
	float dis = sqrt(x * x + y * y + z * z);
	//std::cout << "x: " << x << ", ";
	//std::cout << "y: " << y << ", ";
	//std::cout << "z: " << z << ", ";
	//std::cout << "dis: " << dis << ", ";
	//std::cout << "radius: " << radius << ", ";
	//std::cout << "updown: " << upDownAngle << ", ";
	//std::cout << "leftright: " << leftRightAngle << std::endl;
}

void Assignment2::ComputeAngleFromEye() {
	Eigen::Vector4f eye = sceneData.eye;
	radius = sqrt(eye.x() * eye.x() + eye.y() * eye.y() + eye.z() * eye.z());
	leftRightAngle = eye.z() != 0 ? atan(eye.x() / eye.z()) : 0.0;
	upDownAngle = acos(eye.y() / radius);
}
