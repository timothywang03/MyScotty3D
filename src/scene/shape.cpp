
#include "shape.h"
#include "../geometry/util.h"
#include <iostream>

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const {
	//A3T2 - sphere hit

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

	Vec3 origin = ray.point;
	Vec3 dir = ray.dir;

	float discriminant = 4 * pow(dot(origin, dir), 2) - 4 * (pow(dir.norm(), 2) * (pow(origin.norm(), 2) - pow(radius, 2)));

	PT::Trace ret;
	ret.origin = ray.point;
	ret.hit = false;       // was there an intersection?
	ret.distance = 0.0f;   // at what distance did the intersection occur?
	ret.position = Vec3{}; // where was the intersection?
	ret.normal = Vec3{};   // what was the surface normal at the intersection?
	ret.uv = Vec2{}; 	   // what was the uv coordinates at the intersection? (you may find Sphere::uv to be useful)

	if (discriminant > 0) {
		float t1 = (-2 * dot(origin, dir) - sqrt(discriminant)) / (2 * pow(dir.norm(), 2));
		float t2 = (-2 * dot(origin, dir) + sqrt(discriminant)) / (2 * pow(dir.norm(), 2));

		if (t1 >= ray.dist_bounds[0] && t1 <= ray.dist_bounds[1]) {
			ret.hit = true;
			ret.position = origin + t1 * dir;
			ret.distance = (ret.position - origin).norm();
			ret.normal = (ret.position - Vec3{0, 0, 0}).unit();
			ret.uv = uv(ret.normal);
			return ret;
		} else if (t2 >= ray.dist_bounds[0] && t2 <= ray.dist_bounds[1]) {
			ret.hit = true;
			ret.position = origin + t2 * dir;
			ret.distance = (ret.position - origin).norm();
			ret.normal = (ret.position - Vec3{0, 0, 0}).unit();
			ret.uv = uv(ret.normal);
			return ret;
		} else {
			return ret;
		}

		return ret;
	} else if (discriminant == 0) {
		float t = -2 * dot(origin, dir) / (2 * pow(dir.norm(), 2));
		if (t > ray.dist_bounds.y || t < ray.dist_bounds.x) {
			return ret;
		} else {
			ret.hit = true;
			ret.position = origin + t * dir;
			ret.distance = (ret.position - origin).norm();
			ret.normal = (ret.position - Vec3{0, 0, 0}).unit();
			ret.uv = uv(ret.normal);
			return ret;
		}
	} else {
		return ret;
	}
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}
