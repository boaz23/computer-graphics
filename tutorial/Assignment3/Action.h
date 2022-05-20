

class Action {
public:
	Action(int axis, int targetIndex, int rotationDirection) :
		axis(axis),
		targetIndex(targetIndex),
		rotationDirection(rotationDirection),
		progress(0.0),
		x(-1),
		y(-1),
		isPickingAction(false)
	{};
	Action(double x, double y) :
		axis(-1),
		targetIndex(-1),
		rotationDirection(-1),
		progress(-1),
		x(x),
		y(y),
		isPickingAction(true)
	{};
	Eigen::Vector3d GetRotationAxis() {
		switch (axis) {
		case 0:
			return Eigen::Vector3d(rotationDirection, 0, 0);
		case 1:
			return Eigen::Vector3d(0, rotationDirection, 0);
		case 2:
			return Eigen::Vector3d(0, 0, rotationDirection);
		}
	}
	int axis, targetIndex, rotationDirection;
	double x, y;
	double progress;
	bool isPickingAction;
};