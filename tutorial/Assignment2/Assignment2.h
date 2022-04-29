#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include <tutorial/Assignment2/sceneParser.h>
#include <GLFW/glfw3.h>

#define ANGLE_STEP EIGEN_PI / 50.0

class Assignment2 : public igl::opengl::glfw::Viewer
{
	
public:
	
	Assignment2();
//	Assignment2(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx);
	void WhenRotate();
	void WhenTranslate();
	void RotateScene(int unitsUp, int unitsRight);
	void Animate() override;
	void ScaleAllShapes(float amt, int viewportIndx);
	
	~Assignment2(void);
	void PrintParams();
private:
	SceneData sceneData;
	float upDownAngle;
	float leftRightAngle;
	float radius;
	float scale;
};


