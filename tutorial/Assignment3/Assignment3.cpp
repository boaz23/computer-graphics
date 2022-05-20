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
	speed = SPEED;
	rotationDirection = 1;
	std::random_device rd;
	randomizer = std::mt19937(rd());
	isActive = true;
	offset = wallSize - CUBE_SIZE;
	GenerateCubes();
	ShuffleCubesInitial();
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
			}
		}
	}
}

void Assignment3::ShuffleCubesInitial() {
	std::uniform_int_distribution<int> axisDist(0, 2);
	std::uniform_int_distribution<int> targetDist(0, wallSize - 1);
	for (int i = 0; i < 10; i++) {
		int targetAxis = axisDist(randomizer);
		int targetOffset = -offset + CUBE_SIZE * 2 * targetDist(randomizer);
		Action action(targetAxis, targetOffset, rotationDirection);
		Eigen::Matrix3d fullRotationMatrix = Eigen::AngleAxisd(EIGEN_PI / 2.0, action.GetRotationAxis()).toRotationMatrix();
		for (CubeData* cube : cubesData) {
			if (cube->GetIndexes()(action.axis) == action.targetIndex) {
				data_list[cube->GetMeshId()]->RotateInSystem(action.GetRotationAxis(), EIGEN_PI / 2.0);
				cube->RotateIndexes(fullRotationMatrix.cast<int>());
			}
		}
	}
}
void Assignment3::Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx)
{
	cachedProj = Proj;
	cachedView = View;
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
			if (currentAction->isPickingAction) {
				Picking(currentAction->x, currentAction->y);				
				delete currentAction;
				actionsQueue.pop();
				return;
			}
			double newProgress = currentAction->progress + speed * EIGEN_PI/2.0;
			double angleDelta = std::fmin(speed* EIGEN_PI / 2.0, (EIGEN_PI / 2.0) - currentAction->progress);
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

void Assignment3::AddRotationAction(int axis, int targetIndex)
{
	Action* action = new Action(axis, targetIndex, rotationDirection);
	actionsQueue.push(action);
}

void Assignment3::AddPickingAction(double x, double y)
{
	Action* action = new Action(x, y);
	actionsQueue.push(action);
}

void Assignment3::ProjectScreenCoordToScene(double x, double y, float m_viewport[], const Eigen::Matrix4f& sceneViewInverse, Eigen::Vector3f& sceneSourceCoord, Eigen::Vector3f& sceneDir) {
	Eigen::Vector4f sceneCoord, destCoord;
	sceneCoord << 
		2 * ((x - m_viewport[0]) / m_viewport[2]) - 1.0,
		2 * ((y - m_viewport[1]) / m_viewport[3]) - 1.0,
		0, 1
		;
	destCoord <<
		2 * ((x - m_viewport[0]) / m_viewport[2]) - 1.0,
		2 * ((y - m_viewport[1]) / m_viewport[3]) - 1.0,
		1, 1
		;
	sceneCoord = sceneViewInverse * sceneCoord;
	destCoord = sceneViewInverse * destCoord;
	sceneCoord = sceneCoord / sceneCoord(3);
	destCoord = destCoord / destCoord(3);
	sceneSourceCoord = sceneCoord.head(3);
	sceneDir = (destCoord - sceneCoord).normalized().head(3);
}

bool Assignment3::TriangleIntersection(const Eigen::Vector3f& source, const Eigen::Vector3f& dir, const Eigen::Matrix3f& vertices, const Eigen::Vector3f& normal, Eigen::Vector3f& intersectionPoint) {
	Eigen::Vector3f vv1 = vertices.row(0).transpose() - source;
	float dist = (normal.dot(vv1)) / normal.dot(dir);
	if (dist > FLT_EPSILON && dist < INFINITY) {
		Eigen::Vector3f P = source + dir * dist;
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3f v1 = vertices.row(i);
			Eigen::Vector3f v2 = vertices.row((i + 1) % 3);
			Eigen::Vector3f e1 = v2 - v1;
			Eigen::Vector3f e2 = P - v1;
			Eigen::Vector3f N = e1.cross(e2).normalized();
			if (normal.dot(N) < 0) {
				return false;
			}
		}
		intersectionPoint = P;
		return true;
	}
	return false;
}

bool Assignment3::GetClosestIntersectingFace(CubeData* cubeData, const Eigen::Matrix4f& cubeView, const Eigen::Vector3f& sourcePointScene, const Eigen::Vector3f& dirToScene, Eigen::Vector3f& closestIntersectionToRet, int& closestFaceIndexToRet) {
	igl::opengl::ViewerData* cubeMesh = data_list[cubeData->GetMeshId()];
	int closestShapeFaceIndex = -1;
	float closestShapeFaceDist = -1;
	Eigen::Vector3f closestIntersectionPoint = Eigen::Vector3f::Ones();
	for (int i = 0; i < cubeMesh->F.rows(); i++) {
		Eigen::Matrix3f vertices;
		Eigen::Vector3i vIndexes = cubeMesh->F.row(i);
		vertices.row(0) = cubeMesh->V.row(vIndexes(0)).cast<float>();
		vertices.row(1) = cubeMesh->V.row(vIndexes(1)).cast<float>();
		vertices.row(2) = cubeMesh->V.row(vIndexes(2)).cast<float>();
		Eigen::Vector3f intersectionPoint;
		bool intersect = TriangleIntersection(sourcePointScene, dirToScene, vertices, cubeMesh->F_normals.row(i).cast<float>(), intersectionPoint);
		if (intersect) {
			float dist = (intersectionPoint - sourcePointScene).norm();
			if (closestShapeFaceIndex == -1 || dist < closestShapeFaceDist) {
				closestShapeFaceIndex = i;
				closestShapeFaceDist = dist;
				closestIntersectionPoint = intersectionPoint;
			}
		}
	}
	closestIntersectionToRet = closestIntersectionPoint;
	closestFaceIndexToRet = closestShapeFaceIndex;
	return closestShapeFaceIndex != -1;
}

void Assignment3::Picking(double x, double y) {
	int viewporti[4];
	glGetIntegerv(GL_VIEWPORT, viewporti);
	y = viewporti[3] - y;
	float viewport[] = { viewporti[0], viewporti[1], viewporti[2], viewporti[3] };
	Eigen::Matrix4f sceneView = cachedView * MakeTransScale();
	int closestShapeIndex = -1;
	int closestFaceIndex = -1;
	float closestFaceDist = -1;
	for (CubeData* cube : cubesData) {
		igl::opengl::ViewerData* cubeMesh = data_list[cube->GetMeshId()];
		Eigen::Matrix4f cubeView = sceneView * cubeMesh->MakeTransScale();
		Eigen::Vector3f sourcePointScene, dirToScene;
		ProjectScreenCoordToScene(x, y, viewport, (cachedProj * cubeView).inverse(), sourcePointScene, dirToScene);
		Eigen::Vector3f closestIntersectionPoint;
		int closestShapeFaceIndex;
		if (GetClosestIntersectingFace(cube, cubeView, sourcePointScene, dirToScene, closestIntersectionPoint, closestShapeFaceIndex)) {
			Eigen::Vector4f transformed = cubeView * closestIntersectionPoint.homogeneous();
			if (closestShapeIndex == -1 || transformed(2) > closestFaceDist) {
				closestShapeIndex = cube->GetMeshId();
				closestFaceIndex = closestShapeFaceIndex;
				closestFaceDist = transformed(2);
			}
		}
	}
	if (closestShapeIndex != -1) {
		RotateCubeByFace(cubesData[closestShapeIndex], closestFaceIndex);
	}
}

void Assignment3::RotateCubeByFace(CubeData* cube, int faceIndex) {
	Eigen::Vector3d faceNormal = data_list[cube->GetMeshId()]->F_normals.row(faceIndex);
	faceNormal = data_list[cube->GetMeshId()]->MakeTransd().block(0, 0, 3, 3).matrix() * faceNormal;
	if ((faceNormal - Eigen::Vector3d(1, 0, 0)).norm() <= 0.01) {
		AddRotationAction(0, offset);
	}
	else if ((faceNormal - Eigen::Vector3d(-1, 0, 0)).norm() <= 0.01) {
		AddRotationAction(0, -offset);
	}
	else if ((faceNormal - Eigen::Vector3d(0, 1, 0)).norm() <= 0.01) {
		AddRotationAction(1, offset);
	}
	else if ((faceNormal - Eigen::Vector3d(0, -1, 0)).norm() <= 0.01) {
		AddRotationAction(1, -offset);
	}
	else if ((faceNormal - Eigen::Vector3d(0, 0, 1)).norm() <= 0.01) {
		AddRotationAction(2, offset);
	}
	else if ((faceNormal - Eigen::Vector3d(0, 0, -1)).norm() <= 0.01) {
		AddRotationAction(2, -offset);
	}
	else {
		std::cout << "ilegal normal!!!!!!!!!!!!: " << faceNormal << std::endl;
	}
}

void Assignment3::FlipDirection() {
	rotationDirection = rotationDirection * -1;
}

void Assignment3::IncreaseSpeed() {
	speed += SPEED_STEP;
	if (speed > 1.0) {
		speed = 1.0;
	}
}

void Assignment3::DecreaseSpeed() {
	speed -= SPEED_STEP;
	if (speed <= 0.001) {
		speed = 0.001;
	}
}

void Assignment3::RandomMix() {
	std::uniform_int_distribution<int> axisDist(0, 2);
	std::uniform_int_distribution<int> targetDist(0, wallSize - 1);
	for (int i = 0; i < 10; i++) {
		int targetAxis = axisDist(randomizer);
		int targetOffset = -offset + CUBE_SIZE * 2 * targetDist(randomizer);
		AddRotationAction(targetAxis, targetOffset);
	}
}