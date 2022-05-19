#include "Assignment3.h"
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

Assignment3::Assignment3(int wallSize)
{
	this->wallSize = wallSize;
}

void Assignment3::Init()
{		
	AddShader("shaders/basicShader");
	GenerateCubes();
}

void Assignment3::GenerateCubes() {
	float offset = ((float)wallSize - CUBE_SIZE) / 2;
	bool flag = false;
	for (float x = -offset; x < offset+1; x++) {
		for (float y = -offset; y < offset+1; y++) {
			for (float z = -offset; z < offset+1; z++) {
				if (abs(x) != offset && abs(y) != offset && abs(z) != offset) {
					continue;
				}
				int newShapeIndex;
				if (flag) {
					newShapeIndex = AddShapeCopy(0, -1, TRIANGLES);
				}
				else {
					newShapeIndex = AddShape(Cube, -1, TRIANGLES);
					flag = true;
				}

				SetShapeShader(newShapeIndex, 0);
				Eigen::Vector3d toCenter = Eigen::Vector3d(x, y, z);
				data()->MyTranslate(toCenter, 1);
				data()->SetCenterOfRotation(-toCenter);
				cubesData.push_back(new CubeData(newShapeIndex, Eigen::Vector3d(x, y, z)));
				if (x == -offset+1) {
					Eigen::Matrix3d rotMat = Eigen::AngleAxisd(EIGEN_PI / 4.0, Eigen::Vector3d(1, 0, 0)).toRotationMatrix();
					data()->MyRotate(rotMat);
					cubesData.back()->RotateIndexes(rotMat.cast<int>());
				}
			}
		}
	}
}

void Assignment3::Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx)
{
	Shader *s = shaders[shaderIndx];
	s->Bind();
	s->SetUniformMat4f("Proj", Proj);
	s->SetUniformMat4f("View", View);
	s->SetUniformMat4f("Model", Model);
	s->Unbind();
}


void Assignment3::WhenRotate()
{
}

void Assignment3::WhenTranslate()
{
}

void Assignment3::Animate() {
    if(isActive)
	{
		
	}
}

void Assignment3::ScaleAllShapes(float amt,int viewportIndx)
{
	for (int i = 1; i < data_list.size(); i++)
	{
		if (data_list[i]->Is2Render(viewportIndx))
		{
            data_list[i]->MyScale(Eigen::Vector3d(amt, amt, amt));
		}
	}
}

Assignment3::~Assignment3(void)
{
}

