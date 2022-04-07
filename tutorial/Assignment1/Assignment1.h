#pragma once
#include "igl/opengl/glfw/Viewer.h"

class Assignment1 : public igl::opengl::glfw::Viewer
{
private:
	float time;
	Eigen::Vector4f coeffs;
	Eigen::Vector3cf roots;
	int iterationsNum;
	int currentCoefIndex;
	float translateX, translateY, zoomFactor;

	Eigen::Vector3cf FindRootsOfReduceEquation(Eigen::Vector2cf reduceCoeffs);
	unsigned char* CreateTexture();
	Eigen::Vector4f ComputePixelColor(std::vector<Eigen::Vector4f> colors, Eigen::Vector2f coordinates, Eigen::Vector4f coeef, Eigen::Vector3cf roots, int iterationNum);
	std::complex<float> Assignment1::NewtonCubicRoot(std::complex<float> num);

public:
	Assignment1();
	//	Assignment1(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx);
	void WhenRotate();
	void WhenTranslate();
	void Animate() override;
	void ScaleAllShapes(float amt, int viewportIndx);

	void SetCurrentCoefIndex(int index) { currentCoefIndex = index; }
	void ChangeCurrentCoefBy(float d)
	{
		coeffs(currentCoefIndex) += d;
		roots = FindCubicRoots();
	}
	void ChangeCurrentIterationsNumBy(int diff)
	{
		iterationsNum += diff;
		if (iterationsNum < 0)
		{
			iterationsNum = 0;
		}
	}
	void TranslateX(float dx) { translateX += dx; }
	void TranslateY(float dy) { translateY += dy; }
	void ChangeZoomBy(float d) { zoomFactor += d; }

	Eigen::Vector3cf FindCubicRoots();

	~Assignment1(void);
};

