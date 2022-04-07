#include "Assignment1.h"
#include <iostream>


static void printMat(const Eigen::Matrix4d& mat)
{
	std::cout << " matrix:" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << mat(j, i) << " ";
		std::cout << std::endl;
	}
}

Assignment1::Assignment1()
{
	time = 0;
	coeffs = Eigen::Vector4f::Zero();
}

//Assignment1::Assignment1(float angle ,float relationWH, float near, float far) : Scene(angle,relationWH,near,far)
//{ 	
//}

void Assignment1::Init()
{
	//unsigned int texIDs[1] = { 0 };
	//unsigned int slots[1] = { 0 };

	AddShader("shaders/pickingShader");
	AddShader("shaders/newtonShader");
	//unsigned char* data = CreateTexture();
	//AddTexture(1200, 800, data, COLOR);

	////AddTexture("textures/box0.bmp", 2);

	//AddMaterial(texIDs, slots, 1);
	coeffs = Eigen::Vector4f(1, 1, 0, -1);
	roots = FindCubicRoots();

	AddShape(Plane, -1, TRIANGLES, 0);
	SetShapeShader(0, 1);
	//SetShapeMaterial(0, 0);
	// pickedShape = 0;
	// ShapeTransformation(zTranslate,-5,0);
	// pickedShape = -1;
	SetShapeStatic(0);

}

void Assignment1::Update(const Eigen::Matrix4f& Proj, const Eigen::Matrix4f& View, const Eigen::Matrix4f& Model, unsigned int  shaderIndx, unsigned int shapeIndx)
{
	Shader* s = shaders[shaderIndx];
	int r = ((shapeIndx + 1) & 0x000000FF) >> 0;
	int g = ((shapeIndx + 1) & 0x0000FF00) >> 8;
	int b = ((shapeIndx + 1) & 0x00FF0000) >> 16;
	s->SetUniform1i("iterationNum", 20);
	s->SetUniform4f("coeffs", coeffs(0), coeffs(1), coeffs(2), coeffs(3));
	s->SetUniform4f("rootA", roots(0).real(), roots(0).imag(), 0, 0);
	s->SetUniform4f("rootB", roots(1).real(), roots(1).imag(), 0, 0);
	s->SetUniform4f("rootC", roots(2).real(), roots(2).imag(), 0, 0);
	s->SetUniform4f("colorA", 1, 0, 0, 1);
	s->SetUniform4f("colorB", 0, 1, 0, 1);
	s->SetUniform4f("colorC", 0, 0, 1, 1);

	s->SetUniform1f("time", time);
	s->SetUniform1f("x", x);
	s->SetUniform1f("y", y);
	s->Bind();
	s->SetUniformMat4f("Proj", Proj);
	s->SetUniformMat4f("View", View);
	s->SetUniformMat4f("Model", Model);
	if (data_list[shapeIndx]->GetMaterial() >= 0 && !materials.empty())
	{
		//		materials[shapes[pickedShape]->GetMaterial()]->Bind(textures);
		BindMaterial(s, data_list[shapeIndx]->GetMaterial());
	}
	if (shaderIndx == 0)
		s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 0.0f);
	else
		s->SetUniform4f("lightColor", time / 10.0f, 60 / 100.0f, 99 / 100.0f, 0.5f);
	//textures[0]->Bind(0);
	//s->SetUniform1i("sampler2", materials[shapes[pickedShape]->GetMaterial()]->GetSlot(1));
	//s->SetUniform4f("lightDirection", 0.0f , 0.0f, -1.0f, 0.0f);
//	if(shaderIndx == 0)
//		s->SetUniform4f("lightColor",r/255.0f, g/255.0f, b/255.0f,1.0f);
//	else 
//		s->SetUniform4f("lightColor",0.7f,0.8f,0.1f,1.0f);
	s->Unbind();
}


void Assignment1::WhenRotate()
{
}

void Assignment1::WhenTranslate()
{
}

void Assignment1::Animate() {
	if (isActive)
	{

		time += 0.01f;
	}
}

void Assignment1::ScaleAllShapes(float amt, int viewportIndx)
{
	for (int i = 1; i < data_list.size(); i++)
	{
		if (data_list[i]->Is2Render(viewportIndx))
		{
			data_list[i]->MyScale(Eigen::Vector3d(amt, amt, amt));
		}
	}
}
Eigen::Vector3cf Assignment1::FindCubicRoots()
{
	Eigen::Vector2cf reduceCoeffs = Eigen::Vector2cf::Zero();
	Eigen::Vector3cf roots;
	std::complex<float> bOver3a = (coeffs[1] / coeffs[0]) / 3.0f;
	reduceCoeffs[0] = coeffs[2] / coeffs[0] - 3.0f * bOver3a * bOver3a;
	reduceCoeffs[1] = coeffs[2] / coeffs[0] * bOver3a - coeffs[3] / coeffs[0] - 2.0f * bOver3a * bOver3a * bOver3a;
	std::cout << "reduced\n" << reduceCoeffs << std::endl;
	if (reduceCoeffs.norm() > 0.000001)
	{
		roots = FindRootsOfReduceEquation(reduceCoeffs);
		roots[0] -= bOver3a;
		roots[1] -= bOver3a;
		roots[2] -= bOver3a;
	}
	else
	{
		roots[0] = -1.0f * bOver3a;
		roots[1] = std::complex<float>(std::cosf(3.14159f / 3.0f), std::sinf(3.14159f / 3.0f)) * bOver3a;
		roots[2] = std::complex<float>(std::cosf(2.0f * 3.14159f / 3.0f), std::sinf(2 * 3.14159f / 3.0f)) * bOver3a;
	}

	return roots;
}

std::complex<float> Assignment1::NewtonCubicRoot(std::complex<float> num)
{
	std::complex<float> root = num;
	const int iter = 9;
	bool isSmall = false;
	if (std::abs(num) < 1e-3)
	{
		if (std::abs(num) == 0)
			return num;
		isSmall = true;
		num = num * 1e6f;
		root = num;
	}
	else
		if (std::abs(num) < 0.9f)
			root = 1;
	for (int k = 0; k < iter; k++)
	{
		root = (2.0f * root * root * root + num) / root / root / 3.0f;
	}
	if (isSmall)
		root = root / 100.0f;
	return root;
}


Eigen::Vector3cf Assignment1::FindRootsOfReduceEquation(Eigen::Vector2cf reduceCoeffs)
{
	Eigen::Vector3cf roots = Eigen::Vector3cf::Zero();
	std::complex<float> sqroot = std::sqrt(reduceCoeffs[0] * reduceCoeffs[0] * reduceCoeffs[0] / 27.0f + reduceCoeffs[1] * reduceCoeffs[1] / 4.0f);
	std::complex<float> p = NewtonCubicRoot(reduceCoeffs[1] / 2.0f + sqroot);
	std::complex<float> n = NewtonCubicRoot(reduceCoeffs[1] / 2.0f - sqroot);
	roots[0] = p + n;
	roots[1] = p * std::complex<float>(std::cosf(2.0f * 3.14159f / 3.0f), std::sinf(2 * 3.14159f / 3.0f)) - n * std::complex<float>(std::cosf(1.0f * 3.14159f / 3.0f), std::sinf(1 * 3.14159f / 3.0f));
	roots[2] = -p * std::complex<float>(std::cosf(1.0f * 3.14159f / 3.0f), std::sinf(1 * 3.14159f / 3.0f)) + n * std::complex<float>(std::cosf(2.0f * 3.14159f / 3.0f), std::sinf(2 * 3.14159f / 3.0f));
	return roots;
}

unsigned char* Assignment1::CreateTexture() {
	std::vector<Eigen::Vector4f> colors = {
		Eigen::Vector4f(1, 0, 0, 1),
		Eigen::Vector4f(0, 1, 0, 1),
		Eigen::Vector4f(0, 0, 1, 1)
	};
	int width = 1200;
	int height = 800;
	unsigned char* data = new unsigned char[width * height * 4];
	coeffs = Eigen::Vector4f(1, 1, 1, 1);
	//Eigen::Vector3cf roots = FindCubicRoots();
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			//Eigen::Vector4f pixelColor = ComputePixelColor(colors, Eigen::Vector2f(i, j), coeffs, roots, 20);
			//data[i * width + j] = (unsigned char) pixelColor(0) * 255;
			//data[i * width + j + 1] = (unsigned char)pixelColor(1) * 255;
			//data[i * width + j + 2] = (unsigned char)pixelColor(2) * 255;
			//data[i * width + j + 3] = (unsigned char)pixelColor(3) * 255;
			unsigned char a = 255, b = 0, c = 0, d = 255;
			data[i * width + j*4] = (unsigned char)0xff;
			data[i * width + j*4 + 1] = (unsigned char)0x0;
			data[i * width + j*4 + 2] = (unsigned char)0x0;
			data[i * width + j*4 + 3] = (unsigned char) 0xff;
		}
	}
	unsigned char a = data[0];
	a = data[1];
	a = data[2];
	a = data[3];
	return data;
}

Eigen::Vector4f Assignment1::ComputePixelColor(std::vector<Eigen::Vector4f> colors,  Eigen::Vector2f coordinates, Eigen::Vector4f coeef, Eigen::Vector3cf roots, int iterationNum) {
	float x = coordinates(0);
	float y = coordinates(1);
	std::complex<float> z(x, y);
	for (int i = 0; i < iterationNum; i++) {
		std::complex<float> dx = 3 * coeef(3) * z * z + 2 * coeef(2) * z + coeef(3);
		std::complex<float> currentValue = coeef(3) * z * z * z + coeef(2) * z * z + coeef(1) * z + coeef(0);
		z = z - currentValue / dx;
	}
	int currentMinIndex = 0;
	for (int i = 1; i < 3; i++) {
		if(std::abs(z - roots[i]) < std::abs(z - roots[currentMinIndex])) {
			currentMinIndex = i;
		}
	}
	return colors[currentMinIndex];
}

Assignment1::~Assignment1(void)
{
}
