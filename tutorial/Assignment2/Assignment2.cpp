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
	radius = 1.0f;
	upDownAngle = 0.0;
	leftRightAngle = 0.0;
	cameraCenter = Eigen::Vector4f(0, 0, 0, 0);
	isUp = 1.0f;
	ComputeAngleFromEye();
	//ComputeLookAtMatrix();
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
	Eigen::Vector4f eyeToUpload = ComputeEyePosition();
	s->SetUniform4f("eye", eyeToUpload(0), eyeToUpload(1), eyeToUpload(2), eyeToUpload(3));	
	s->SetUniform4f("ambient", sceneData.ambient(0), sceneData.ambient(1), sceneData.ambient(2), sceneData.ambient(3));
	s->SetUniform4fv("objects", sceneData.objects.data(), sceneData.objects.size());
	s->SetUniform4fv("objColors", sceneData.colors.data(), sceneData.colors.size());
	s->SetUniform4fv("lightsDirection", sceneData.directions.data(), sceneData.directions.size());
	s->SetUniform4fv("lightsIntensity", sceneData.intensities.data(), sceneData.intensities.size());
	s->SetUniform4fv("lightsPosition", sceneData.lights.data(), sceneData.lights.size());
	s->SetUniform4i("sizes", sceneData.sizes(0), sceneData.sizes(1), sceneData.sizes(2), sceneData.sizes(3));
	s->SetUniform1f("upDownAngle", upDownAngle);
	s->SetUniform1f("leftRightAngle", leftRightAngle);
	s->SetUniform4f("cameraCenter", cameraCenter(0), cameraCenter(1), cameraCenter(2), cameraCenter(3));
	s->SetUniformMat4f("Proj", Proj);
	s->SetUniformMat4f("MyProj", GetProjectionMatrix());
	s->SetUniformMat4f("LookAt", GetLookAtMatrix());
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
	Eigen::Vector4f forwardVector = (cameraCenter - ComputeEyePosition()).normalized();
	Eigen::Vector4f rightVector = forwardVector.cross3(Eigen::Vector4f(0, isUp, 0, 0)).normalized();
	Eigen::Vector4f upVector = rightVector.cross3(forwardVector).normalized();
	float dx = xRel * 2;
	float dy = yRel * 2;
	Eigen::Vector4f translation = rightVector * dx + upVector * dy;
	cameraCenter = cameraCenter + translation;
	//sceneData.eye = sceneData.eye + translation;
	//ComputeAngleFromEye();
}

void Assignment2::TranslateX(float dx)
{
	sceneData.eye.x() += xRel * 2 * radius;
}

void Assignment2::TranslateY(float dy)
{
	sceneData.eye.y() += yRel * 2 * radius;
}

void Assignment2::ChangeZoomBy(float dz)
{
	//float nextZoom = radius / (1 + SCALE_SENSITIVITY * dz);
	//if (nextZoom >= 0) {
	//	radius = nextZoom;
	//	std::cout << ComputeEyePosition() << std::endl;
	//	std::cout << "**********" << std::endl;
	//	//ComputeEyeFromAngle();
	//}
	Eigen::Vector4f forwardVector = (cameraCenter - ComputeEyePosition()).normalized();
	Eigen::Vector4f translation = forwardVector * SCALE_SENSITIVITY * dz;
	cameraCenter = cameraCenter + translation;
}


void Assignment2::RotateScene(int unitsUp, int unitsRight) {
	if (isUp > 0.0f) {
		leftRightAngle -= unitsRight * ANGLE_STEP;
	}
	else {
		leftRightAngle += unitsRight * ANGLE_STEP;
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
	//upDownAngle += unitsUp * ANGLE_STEP;
	//leftRightAngle += unitsRight * ANGLE_STEP;
	//if (upDownAngle <= 0) {
	//	upDownAngle += 2*EIGEN_PI;
	//}
	//else if (upDownAngle >= 2*EIGEN_PI) {
	//	upDownAngle -= 2*EIGEN_PI;
	//}
	//if (leftRightAngle <= 0) {
	//	leftRightAngle += 2*EIGEN_PI;
	//}
	//else if (leftRightAngle >= 2*EIGEN_PI) {
	//	leftRightAngle -= 2*EIGEN_PI;
	//}
	//ComputeEyeFromAngle();
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
	//vec3 vRay = normalize(position0 + cameraCenter.xyz - eye.xyz);
	//StraightLine ray = StraightLine(eye.xyz, vRay);
	Eigen::Vector4f eye = ComputeEyePosition();
	Eigen::Vector3f v = (cursorPoint + cameraCenter.head(3) - eye.head(3)).normalized();
	// burn with fire
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

Eigen::Vector4f Assignment2::ComputeEyePosition() {
	float m_radius = radius;
	return Eigen::Vector4f(m_radius * sinf(upDownAngle) * sinf(leftRightAngle),
		m_radius * cosf(upDownAngle),
		m_radius * sinf(upDownAngle) * cosf(leftRightAngle),
		1) + cameraCenter;
}
void Assignment2::ComputeEyeFromAngle() {
	float x = radius * sinf(upDownAngle) * sinf(leftRightAngle) + cameraCenter.x();
	float y = radius * cosf(upDownAngle) + cameraCenter.y();
	float z = radius * sinf(upDownAngle) * cosf(leftRightAngle) + cameraCenter.z();
	sceneData.eye.x() = x;
	sceneData.eye.y() = y;
	sceneData.eye.z() = z;
	ComputeLookAtMatrix();
	//float dis = sqrt(x * x + y * y + z * z);
	//std::cout << "x: " << x << ", ";
	//std::cout << "y: " << y << ", ";
	//std::cout << "z: " << z << ", ";
	//std::cout << "dis: " << dis << ", ";
	//std::cout << "radius: " << radius << ", ";
	//std::cout << "updown: " << upDownAngle << ", ";
	//std::cout << "leftright: " << leftRightAngle << std::endl;
}

void Assignment2::ComputeAngleFromEye() {
	Eigen::Vector4f eyeToCamera = sceneData.eye - cameraCenter;
	//radius = eyeToCamera.norm();
	leftRightAngle = eyeToCamera.z() != 0 ? atanf(eyeToCamera.x() / eyeToCamera.z()) : 0.0;
	upDownAngle = acosf(eyeToCamera.y() / radius);
	//ComputeLookAtMatrix();
}

void Assignment2::TransformObject() {
	if (pickedObjectIndex != -1) {
		Eigen::Vector4f forwardVector = (cameraCenter - ComputeEyePosition()).normalized();
		Eigen::Vector4f rightVector = forwardVector.cross3(Eigen::Vector4f(0, isUp, 0, 0)).normalized();
		Eigen::Vector4f upVector = rightVector.cross3(forwardVector).normalized();
		float dx = xRel;
		float dy = yRel;
		Eigen::Vector4f translation = rightVector * dx + upVector * dy;
		if (sceneData.objects[pickedObjectIndex](3) > 0) {
			//sphere
			//sceneData.objects[pickedObjectIndex](0) += xRel * 2;
			//sceneData.objects[pickedObjectIndex](1) += yRel * 2;
			sceneData.objects[pickedObjectIndex] += translation * 4;
		}
		else {
			//plane
			//sceneData.objects[pickedObjectIndex](0) -= xRel * 2 * radius;
			//sceneData.objects[pickedObjectIndex](1) -= yRel * 2 * radius;
			sceneData.objects[pickedObjectIndex] -= translation;
		}
	}
}

void Assignment2::ComputeLookAtMatrix() {
	Eigen::Vector4f forwardVector = (cameraCenter - sceneData.eye).normalized();
	Eigen::Vector4f rightVector = forwardVector.cross3(Eigen::Vector4f(0, isUp, 0, 0)).normalized();
	Eigen::Vector4f upVector = rightVector.cross3(forwardVector).normalized();
	lookAtMatrix.row(0) = rightVector;
	lookAtMatrix.row(1) = upVector;
	lookAtMatrix.row(2) = forwardVector;
	//lookAtMatrix.row(3) = Eigen::Vector4f(
	//	-rightVector.transpose() * sceneData.eye, 
	//	-upVector.transpose() * sceneData.eye,
	//	-forwardVector.transpose() * sceneData.eye,
	//	1);
	lookAtMatrix.row(3) = sceneData.eye;
}

Eigen::Matrix4f Assignment2::GetLookAtMatrix() {
	Eigen::Vector3f eye = ComputeEyePosition().head(3);
	Eigen::RowVector3f forwardVector = (cameraCenter.head(3) - eye.head(3)).normalized();
	Eigen::RowVector3f rightVector = forwardVector.cross(Eigen::Vector3f(0, isUp, 0)).normalized();
	Eigen::RowVector3f upVector = rightVector.cross(forwardVector).normalized();
	Eigen::Matrix4f toRet = Eigen::Matrix4f::Identity();
	toRet << 
		rightVector, 0,
		upVector, 0,
		forwardVector, 0,
		-rightVector * eye, -upVector * eye, forwardVector * eye, 1;
	//toRet.row(0) = rightVector;
	//toRet.row(1) = upVector;
	//toRet.row(2) = forwardVector;
	//toRet.row(3) = Eigen::Vector4f(
	//	-rightVector.transpose() * eye, 
	//	-upVector.transpose() * eye,
	//	-forwardVector.transpose() * eye,
	//	1);
	return toRet;
	//lookAtMatrix.row(3) = sceneData.eye;
}

Eigen::Matrix4f Assignment2::GetProjectionMatrix() {
	Eigen::Matrix4f toRet;
	float l = 0.0f, r = 1200.0f,
		b = 0.0f, t = 800.0f,
		n = 1.0f, f = 120.0f;
	toRet << 2.0f * n / (r - l), 0, (r + l) / (r - l), 0,
		0, 2.0f * n / (t - b), (t + b) / (t - b), 0,
		0, 0, -(f + n) / (f - n), -2.0f * f * n / (f - n),
		0, 0, -1, 0;
	return toRet;
}