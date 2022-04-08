#pragma once
#include "igl/opengl/glfw/Viewer.h"
#define SCALE_SENSITIVITY 0.01
#define TRANSLATION_SENSITIVITY 0.01
#define WIDTH 1200
#define HEIGHT 800
class Assignment1 : public igl::opengl::glfw::Viewer
{
private:
	float time;
	Eigen::Vector4f coeffs;
	Eigen::Vector3cf roots;
	int iterationsNum;
	int currentCoefIndex;
	float translateX, translateY, zoomNormalized;
	float xOld, yOld, xRel, yRel;

	Eigen::Vector3cf FindRootsOfReduceEquation(Eigen::Vector2cf reduceCoeffs);
	std::complex<float> Assignment1::NewtonCubicRoot(std::complex<float> num);

public:
	Assignment1();
	//	Assignment1(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx);
	void WhenRotate(const Eigen::Matrix4d& preMat, float dx, float dy);
	void WhenTranslate();
	float UpdatePosition(float xpos, float ypos);
	void Animate() override;
	void ScaleAllShapes(float amt, int viewportIndx);

	void SetCurrentCoefIndex(int index);
	void ChangeCurrentCoefBy(float d);
	void ChangeCurrentIterationsNumBy(int diff);
	void ResetZoom();
	void TranslateX(float dx);
	void TranslateY(float dy);
	void ChangeZoomBy(float d);
	void PrintCurrentCoefficient();
	float GetZoom() { return zoomNormalized; };

	Eigen::Vector3cf FindCubicRoots();

	~Assignment1(void);
};

