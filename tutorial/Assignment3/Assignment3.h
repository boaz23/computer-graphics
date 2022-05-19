#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include "./CubeData.h";

#define CUBE_SIZE 1.0f

class Assignment3 : public igl::opengl::glfw::Viewer
{
	
public:
	
	Assignment3(int wallSize);
//	Assignment3(float angle,float relationWH,float near, float far);
	void Init();
	void GenerateCubes();
	void Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx);
	void WhenRotate();
	void WhenTranslate();
	void Animate() override;
	void ScaleAllShapes(float amt, int viewportIndx);
	
	~Assignment3(void);

private:
	int wallSize;
	std::vector<CubeData*> cubesData;
};


