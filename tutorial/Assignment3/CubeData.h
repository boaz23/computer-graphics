#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>

class CubeData {
public:
	CubeData(int shapeIndex, Eigen::Vector3i indexes) :
		shapeIndex(shapeIndex),
		originalIndexes(indexes),
		indexes(indexes) {};

	void RotateIndexes(Eigen::Matrix3i rotationMat) {
		indexes = rotationMat * indexes;
	}

	Eigen::Vector3i GetIndexes() {
		return indexes;
	}

	int GetMeshId() {
		return shapeIndex;
	}

private:
	int shapeIndex;
	Eigen::Vector3i originalIndexes, indexes;
};