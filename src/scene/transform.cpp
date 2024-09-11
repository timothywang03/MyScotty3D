
#include "transform.h"

Mat4 Transform::local_to_parent() const {
	return Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
}

Mat4 Transform::parent_to_local() const {
	return Mat4::scale(1.0f / scale) * rotation.inverse().to_mat() * Mat4::translate(-translation);
}

Mat4 Transform::local_to_world() const {
	// A1T1: local_to_world
	//don't use Mat4::inverse() in your code.

	Mat4 r = rotation.to_mat();

	// define the default as the identity matrix
	Mat4 m = Mat4::I;

	if (std::shared_ptr<Transform> parent_ = parent.lock()) {
		m = parent_->local_to_world();
	}

	m = m * Mat4::translate(translation) * r * Mat4::scale(scale);
	return m; //<-- wrong, but here so code will compile
}

Mat4 Transform::world_to_local() const {
	// A1T1: world_to_local
	//don't use Mat4::inverse() in your code.

	// apply the inverse transformations in reverse order
	Vec3 t = -translation;
	Mat4 r = rotation.inverse().to_mat();
	Vec3 s = 1.0f / scale;

	// define the default as the identity matrix
	Mat4 m = Mat4::I;

	if (std::shared_ptr<Transform> parent_ = parent.lock()) {
		m = parent_->world_to_local();
	}

	m = Mat4::scale(s) * r * Mat4::translate(t) * m;
	return m; //<-- wrong, but here so code will compile
}

bool operator!=(const Transform& a, const Transform& b) {
	return a.parent.lock() != b.parent.lock() || a.translation != b.translation ||
	       a.rotation != b.rotation || a.scale != b.scale;
}
