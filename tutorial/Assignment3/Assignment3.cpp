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
	isActive = true;
	offset = wallSize - CUBE_SIZE;
	GenerateCubes();
}

void Assignment3::GenerateCubes() {
	bool flag = false;
	for (int x = -offset; x <= offset; x += CUBE_SIZE * 2) {
		for (int y = -offset; y <= offset; y += CUBE_SIZE * 2) {
			for (int z = -offset; z <= offset; z += CUBE_SIZE * 2) {
				if (abs(x) != offset && abs(y) != offset && abs(z) != offset) {
					// internal cube. can skip rendering it.
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
				Eigen::Vector3d toCenter = Eigen::Vector3d(x, y, z).array() / 2.0;
				data()->MyTranslate(toCenter, 1);
				data()->SetCenterOfRotation(-toCenter);
				cubesData.push_back(new CubeData(newShapeIndex, Eigen::Vector3i(x, y, z)));
				if (x == -offset && false) {
					Eigen::Matrix3d rotMat = Eigen::AngleAxisd(EIGEN_PI / 4.0, Eigen::Vector3d(1, 0, 0)).toRotationMatrix();
					data()->MyRotate(rotMat);
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
		if (!actionsQueue.empty()) {
			Action* currentAction = actionsQueue.front();
			double newProgress = currentAction->progress + SPEED * EIGEN_PI/2.0;
			double angleDelta = std::fmin(SPEED* EIGEN_PI / 2.0, (EIGEN_PI / 2.0) - currentAction->progress);
			for (CubeData* cube : cubesData) {
				if (cube->GetIndexes()(currentAction->axis) == currentAction->targetIndex) {
					data_list[cube->GetMeshId()]->RotateInSystem(currentAction->GetRotationAxis(), angleDelta);
				}
			}
			currentAction->progress = currentAction->progress + angleDelta;
			if (newProgress >= EIGEN_PI / 2.0) {
				Eigen::Matrix3d fullRotationMatrix = Eigen::AngleAxisd(EIGEN_PI / 2.0, currentAction->GetRotationAxis()).toRotationMatrix();
				for (CubeData* cube : cubesData) {
					if (cube->GetIndexes()(currentAction->axis) == currentAction->targetIndex) {
						cube->RotateIndexes(fullRotationMatrix.cast<int>());
					}
				}
				delete currentAction;
				actionsQueue.pop();
			}
		}
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
	for (CubeData* cube : cubesData) {
		delete cube;
	}
	while (!actionsQueue.empty())
	{
		Action* action = actionsQueue.front();
		actionsQueue.pop();
		delete action;
	}
}

void Assignment3::AddAction(int axis, int targetIndex)
{
	Action* action = new Action(axis, targetIndex);
	actionsQueue.push(action);
}