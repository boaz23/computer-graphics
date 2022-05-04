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
	AddShape(Plane, -1, TRIANGLES, 0);
	xOld = 0;
	yOld = 0;
	xRel = 0;
	yRel = 0;

	SetShapeShader(0, 0);
	pickedShape = 0;
	SetShapeStatic(0);
	radius = 1.0f;
	upDownAngle = 0.0;
	leftRightAngle = 0.0;
	cameraCenter = Eigen::Vector4f(0, 0, 0, 0);
	isUp = 1.0f;
	ComputeAngleFromEye();
	ComputeLookAtMatrix();
}

void Assignment2::Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((shapeIndx+1) & 0x000000FF) >>  0;
	int g = ((shapeIndx+1) & 0x0000FF00) >>  8;
	int b = ((shapeIndx+1) & 0x00FF0000) >> 16;
	s->Bind();
	Eigen::Vector4f eyeToUpload = sceneData.eye;
	s->SetUniform4f("eye", eyeToUpload(0), eyeToUpload(1), eyeToUpload(2), eyeToUpload(3));	
	s->SetUniform4f("ambient", sceneData.ambient(0), sceneData.ambient(1), sceneData.ambient(2), sceneData.ambient(3));
	s->SetUniform4fv("objects", sceneData.objects.data(), sceneData.objects.size());
	s->SetUniform4fv("objColors", sceneData.colors.data(), sceneData.colors.size());
	s->SetUniform4fv("lightsDirection", sceneData.directions.data(), sceneData.directions.size());
	s->SetUniform4fv("lightsIntensity", sceneData.intensities.data(), sceneData.intensities.size());
	s->SetUniform4fv("lightsPosition", sceneData.lights.data(), sceneData.lights.size());
	s->SetUniform4i("sizes", sceneData.sizes(0), sceneData.sizes(1), sceneData.sizes(2), sceneData.sizes(3));
	s->SetUniform4f("cameraCenter", cameraCenter(0), cameraCenter(1), cameraCenter(2), cameraCenter(3));
	s->SetUniformMat4f("Proj", Proj);
	s->SetUniformMat4f("LookAt", lookAtMatrix);
	s->SetUniformMat4f("View", View);
	s->SetUniformMat4f("Model", Model);
	s->Unbind();
}


void Assignment2::WhenRotate()
{
}

float Assignment2::UpdatePosition(float xpos, float ypos)
{
	xRel = -xOld + xpos / 1200;
	yRel = yOld - ypos / 800;
	xOld += xRel;
	yOld -= yRel;
	return yRel;
}

void Assignment2::WhenTranslate()
{
	Eigen::Vector3f eye = sceneData.eye.head(3);
	Eigen::RowVector3f forwardVector = (cameraCenter.head(3) - eye).normalized();
	Eigen::RowVector3f rightVector = forwardVector.cross(Eigen::Vector3f(0, isUp, 0)).normalized();
	Eigen::RowVector3f upVector = rightVector.cross(forwardVector).normalized();
	float dx = xRel * 2;
	float dy = yRel * 2;
	Eigen::RowVector4f translation;
	translation << (rightVector * dx + upVector * dy), 0;
	cameraCenter = cameraCenter + translation.transpose();
	ComputeEyeFromAngle();
}

void Assignment2::ChangeZoomBy(float dz)
{
	Eigen::Vector3f eye = sceneData.eye.head(3);
	Eigen::RowVector3f forwardVector = (cameraCenter.head(3) - eye).normalized();
	Eigen::RowVector4f translation;
	translation << (forwardVector * SCALE_SENSITIVITY * dz), 0;
	cameraCenter = cameraCenter + translation.transpose();
	ComputeEyeFromAngle();
}


void Assignment2::RotateScene(int unitsUp, int unitsRight) {
	if (isUp > 0.0f) {
		leftRightAngle += unitsRight * ANGLE_STEP;
	}
	else {
		leftRightAngle -= unitsRight * ANGLE_STEP;
	}
	if (leftRightAngle > 2 * EIGEN_PI) {
		leftRightAngle -= 2 * EIGEN_PI;
	}
	else if (leftRightAngle < -2 * EIGEN_PI) {
		leftRightAngle += 2 * EIGEN_PI;
	}
	upDownAngle -= unitsUp * ANGLE_STEP;
	if (upDownAngle > 2 * EIGEN_PI) {
		upDownAngle -= 2 * EIGEN_PI;
	}
	else if (upDownAngle < -2 * EIGEN_PI) {
		upDownAngle += 2 * EIGEN_PI;
	}
	if ((upDownAngle > 0 && upDownAngle < EIGEN_PI) || (upDownAngle < -EIGEN_PI && upDownAngle > -2 * EIGEN_PI)) {
		isUp = 1.0f;
	}
	else {
		isUp = -1.0f;
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

float Assignment2::intersectionWithObject(Eigen::Vector4f objectPos, Eigen::Vector3f sourcePoint, Eigen::Vector3f dir) {
	float dist = -1.0;
	if (objectPos(3) <= 0) {
		//plane
		float d = objectPos(3);
		Eigen::Vector3f perpendicular = objectPos.head(3);
		dist = -(perpendicular.dot(sourcePoint) + d) / perpendicular.dot(dir);
	}
	else {
		//sphere
		Eigen::Vector3f o = objectPos.head(3);
		float r = objectPos(3);
		Eigen::Vector3f L = o - sourcePoint;
		float tm = L.dot(dir);
		float dSquared = (L.norm() * L.norm()) - (tm * tm);
		if (dSquared <= r * r) {
			float th = sqrt(r * r - dSquared);
			float t1 = tm - th, t2 = tm + th;
			if (t1 > 0.001) {
				dist = t1;
			}
			else if (t2 > 0.001) {
				dist = t2;
			}
			else {
				dist = -1.0;
			}
		}
	}
	return dist;
}
void Assignment2::intersection(Eigen::Vector3f cursorPoint) {
	Eigen::Vector4f eye = sceneData.eye;
	Eigen::Vector3f v = (cursorPoint + cameraCenter.head(3) - eye.head(3)).normalized();
	cursorPoint = eye.head(3);
	int minIndex = -1;
	float minDist = -1;
	for (int i = 0; i < sceneData.sizes[0]; i++) {
		Eigen::Vector4f curObject = sceneData.objects[i];
		float dist = intersectionWithObject(curObject, cursorPoint, v);
		if (true
			&& dist >= 0.001
			&& (minIndex == -1 || dist < minDist)
			) {
			minIndex = i;
			minDist = dist;
		}
	}
	pickedObjectIndex = minIndex;
}
Assignment2::~Assignment2(void)
{

}

void Assignment2::ComputeEyeFromAngle() {
	float x = radius * sinf(upDownAngle) * sinf(leftRightAngle) + cameraCenter.x();
	float y = radius * cosf(upDownAngle) + cameraCenter.y();
	float z = radius * sinf(upDownAngle) * cosf(leftRightAngle) + cameraCenter.z();
	sceneData.eye.x() = x;
	sceneData.eye.y() = y;
	sceneData.eye.z() = z;
	ComputeLookAtMatrix();
}

void Assignment2::ComputeAngleFromEye() {
	Eigen::Vector4f eyeToCamera = sceneData.eye - cameraCenter;
	leftRightAngle = eyeToCamera.z() != 0 ? atanf(eyeToCamera.x() / eyeToCamera.z()) : 0.0;
	upDownAngle = acosf(eyeToCamera.y() / radius);
}

void Assignment2::TransformObject() {
	if (pickedObjectIndex != -1) {
		Eigen::Vector3f eye = sceneData.eye.head(3);
		Eigen::RowVector3f forwardVector = (cameraCenter.head(3) - eye).normalized();
		Eigen::RowVector3f rightVector = forwardVector.cross(Eigen::Vector3f(0, isUp, 0)).normalized();
		Eigen::RowVector3f upVector = rightVector.cross(forwardVector).normalized();
		float dx = xRel;
		float dy = yRel;
		Eigen::RowVector4f translation;
		translation << (rightVector * dx + upVector * dy), 0;
		if (sceneData.objects[pickedObjectIndex](3) > 0) {
			//sphere
			sceneData.objects[pickedObjectIndex] += translation.transpose() * 8;
		}
		else {
			//plane
			sceneData.objects[pickedObjectIndex] -= translation.transpose();
		}
	}
}

void Assignment2::ComputeLookAtMatrix() {
	Eigen::Vector3f eye = sceneData.eye.head(3);
	Eigen::RowVector3f forwardVector = (cameraCenter.head(3) - eye.head(3)).normalized();
	Eigen::RowVector3f rightVector = forwardVector.cross(Eigen::Vector3f(0, isUp, 0)).normalized();
	Eigen::RowVector3f upVector = rightVector.cross(forwardVector).normalized();
	Eigen::Matrix4f toRet = Eigen::Matrix4f::Identity();
	toRet <<
		rightVector, 0,
		upVector, 0,
		forwardVector, 0,
		-rightVector * eye, -upVector * eye, forwardVector* eye, 1;
	lookAtMatrix = toRet;
}
