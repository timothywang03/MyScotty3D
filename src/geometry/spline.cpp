
#include "../geometry/spline.h"
#include <iostream>

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

	// edge case handling
	if (knots.empty()) return T();
	if (knots.size() == 1) return knots.begin()->second;
	if (knots.begin()->first > time) return knots.begin()->second;
	auto upperbound_knot = knots.upper_bound(time);
	if (upperbound_knot == knots.end()) return std::prev(knots.end())->second;

	auto lower_iter = std::prev(knots.upper_bound(time));
	auto upper_iter = knots.upper_bound(time);

	float t1 = lower_iter->first;
	float t2 = upper_iter->first;
	float t0, t3;

	T p1 = lower_iter->second;
	T p2 = upper_iter->second;
	T p0, p3;

	if (lower_iter == knots.begin()) {
		t0 = t1 - (t2 - t1);
		p0 = p1 - (p2 - p1);
	} else {
		t0 = std::prev(lower_iter)->first;
		p0 = std::prev(lower_iter)->second;
	}

	if (std::next(upper_iter) == knots.end()) {
		t3 = t2 + (t2 - t1);
		p3 = p2 + (p2 - p1);
	} else {
		t3 = std::next(upper_iter)->first;
		p3 = std::next(upper_iter)->second;
	}

	T m0 = (p2 - p0) / (t2 - t0);
	T m1 = (p3 - p1) / (t3 - t1);

	// std::cout << "a: " << (t2 - t0) << std::endl;
	// std::cout << "b: " << (t3 - t1) << std::endl;
	// std::cout << "c: " << (t2 - t1) << std::endl;

	// std::cout << float((time - t1)) / float((t2 - t1)) << std::endl
	return cubic_unit_spline((time - t1) / (t2 - t1), p1, p2, m0, m1);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.

	float h00 = 2 * time * time * time - 3 * time * time + 1;
	float h10 = time * time * time - 2 * time * time + time;
	float h01 = -2 * time * time * time + 3 * time * time;
	float h11 = time * time * time - time * time;

	return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
