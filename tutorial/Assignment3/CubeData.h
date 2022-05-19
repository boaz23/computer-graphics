#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>

class CubeData {
public:
	CubeData(int shapeIndex, Eigen::Vector3d indexes) :
		shapeIndex(shapeIndex),
		indexes(indexes) {};

	void RotateIndexes(Eigen::Matrix3i rotationMat) {
		indexes = rotationMat.cast<double>() * indexes;
	}

private:
	int shapeIndex, offset;
	Eigen::Vector3d indexes;
};