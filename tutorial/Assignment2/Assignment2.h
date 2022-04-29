#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include <tutorial/Assignment2/sceneParser.h>
#include <GLFW/glfw3.h>

#define ANGLE_STEP EIGEN_PI / 200.0
#define SCALE_SENSITIVITY 0.01
#define TRANSLATION_SENSITIVITY 0.01

class Assignment2 : public igl::opengl::glfw::Viewer
{
	
public:
	
	Assignment2();
//	Assignment2(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx);
	void WhenRotate();
	void WhenTranslate();
	void ComputeAngleFromEye();
	void TransformObject();
	void RotateScene(int unitsUp, int unitsRight);
	void Animate() override;
	void ScaleAllShapes(float amt, int viewportIndx);

	float intersectionWithObject(Eigen::Vector4f objectPos, Eigen::Vector3f sourcePoint, Eigen::Vector3f dir);

	void intersection(Eigen::Vector3f cursorPoint);
	
	~Assignment2(void);
	void ComputeEyeFromAngle();
	float UpdatePosition(float xpos, float ypos);
	void TranslateX(float dx);
	void TranslateY(float dy);
	void ChangeZoomBy(float d);
private:
	SceneData sceneData;
	float upDownAngle;
	float leftRightAngle;
	float radius;
	float scale;
	float xOld, yOld, xRel, yRel;
	int pickedObjectIndex;
};


