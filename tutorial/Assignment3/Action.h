

class Action {
public:
	Action(int axis, int targetIndex) :
		axis(axis),
		targetIndex(targetIndex),
		progress(0.0)
	{};
	Eigen::Vector3d GetRotationAxis() {
		switch (axis) {
		case 0:
			return Eigen::Vector3d(1, 0, 0);
		case 1:
			return Eigen::Vector3d(0, 1, 0);
		case 2:
			return Eigen::Vector3d(0, 0, 1);
		}
	}
	int axis, targetIndex;
	double progress;
};