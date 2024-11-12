
#include "material.h"
#include "../util/rand.h"
#include <iostream>

namespace Materials {

Vec3 reflect(Vec3 dir) {
	//A3T5 Materials - reflect helper

    // Return direction to incoming light that would be
	// reflected out in direction dir from surface
	// with normal (0,1,0)

	Vec3 normal = Vec3(0, 1, 0);
	// std::cout << "dir: " << dir << std::endl;
	float costheta = dot(dir, normal) / dir.norm();
	Vec3 reflection = -dir + 2 * costheta * normal;

    return reflection;
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {
    // The surface normal is (0,1,0)
    Vec3 normal(0, 1, 0);
    
    float etaA = 1, etaB = index_of_refraction;
    float etaI, etaT;
    
    // Determine if we're entering or exiting the material
    if (out_dir.y > 0) {
        etaI = etaA;
        etaT = etaB;
    } else {
        etaI = etaB;
        etaT = etaA;
        normal = -normal; // Flip normal if we're exiting
    }
    
    float eta = etaI / etaT;
    float costheta = std::abs(dot(out_dir, normal));
    float internal = (eta * eta) * (1.0f - costheta * costheta);
    
    // Check for total internal reflection
    if (internal > 1.0f) {
        was_internal = true;
        return Vec3(); 
    }
    
    was_internal = false;
    float theta_out = std::sqrt(1.0f - internal);
    return eta * -out_dir + (eta * costheta - theta_out) * normal;
}

float schlick(Vec3 in_dir, float index_of_refraction) {
	//A3T5 Materials - Schlick's approximation helper

	// Implement Schlick's approximation of the Fresnel reflection factor.

	Vec3 norm = Vec3(0, 1, 0);
	float costheta = abs(dot(in_dir, norm) / in_dir.norm());
	float n1 = 1.0f;
	float n2 = index_of_refraction;
	float R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
	float R = R0 + (1 - R0) * pow(1 - costheta, 5);

	return R;
}

Spectrum Lambertian::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF evaluation

    // Compute the ratio of outgoing/incoming radiance when light from in_dir
    // is reflected through out_dir: (albedo / PI_F) * cos(theta).
    // Note that for Scotty3D, y is the 'up' direction.
	Spectrum albedo = Lambertian::albedo.lock()->evaluate(uv);
	float theta = dot(in, Vec3(0, 1, 0)) / in.norm();
	if (theta <= 0) return Spectrum();
    return (albedo / PI_F) * in.y;
}

Scatter Lambertian::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF scattering
	//Select a scattered light direction at random from the Lambertian BSDF

	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this will be useful

	Scatter ret;
	Samplers::Hemisphere::Cosine cosine_sampler;
	ret.direction = cosine_sampler.sample(rng);

	ret.attenuation = evaluate(out, ret.direction, uv);

	return ret;
}

float Lambertian::pdf(Vec3 out, Vec3 in) const {
	//A3T4: Materials - Lambertian BSDF probability density function
    // Compute the PDF for sampling in_dir from the cosine-weighted hemisphere distribution.
	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this might be handy!

    return sampler.pdf(in);
}

Spectrum Lambertian::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Lambertian::display() const {
	return albedo;
}

void Lambertian::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(albedo);
}

Spectrum Mirror::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return reflectance.lock()->evaluate(uv);
}

Scatter Mirror::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5: mirror

	// Use reflect to compute the new direction
	// Don't forget that this is a discrete material!
	// Similar to albedo, reflectance represents the ratio of incoming light to reflected light

    Scatter ret;
    ret.direction = reflect(out);
    ret.attenuation = evaluate(out, ret.direction, uv);
    return ret;
}

float Mirror::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Mirror::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Mirror::display() const {
	return reflectance;
}

void Mirror::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
}

Spectrum Refract::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	Spectrum spec = transmittance.lock()->evaluate(uv);
    float new_ior = ior;
    if (out.y < 0) {
        new_ior = 1/ior;
    }
    return (new_ior * new_ior) * spec;
}

Scatter Refract::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - refract

	// Use refract to determine the new direction - what happens in the total internal reflection case?
    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
	// Don't forget that this is a discrete material!
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

	bool internal = false;	
	Vec3 refraction = refract(out, ior, internal);

    Scatter ret;
    if (internal) {
		ret.direction = reflect(out);
	} else {
		ret.direction = refraction;
	}
    ret.attenuation = evaluate(out, ret.direction, uv);
    return ret;
}

float Refract::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Refract::emission(Vec2 uv) const {
	return {};
}

bool Refract::is_emissive() const {
	return false;
}

bool Refract::is_specular() const {
	return true;
}

bool Refract::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Refract::display() const {
	return transmittance;
}

void Refract::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(transmittance);
}

Spectrum Glass::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return transmittance.lock()->evaluate(uv);
}

Scatter Glass::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - glass

    // (1) Compute Fresnel coefficient. Tip: Schlick's approximation.
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
    // What happens upon total internal reflection?
    // When debugging Glass, it may be useful to compare to a pure-refraction BSDF
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

	float ratio = schlick(out, ior);
	bool internal = false;
	Vec3 refraction = refract(out, ior, internal);

    Scatter ret;
	if (internal) {
		ret.direction = reflect(out);
	} else {
		if (rng.coin_flip(ratio)) {
			ret.direction = reflect(out);
		} else {
			ret.direction = refraction;
		}
	}
    ret.attenuation = evaluate(out, ret.direction, uv);
    return ret;
}

float Glass::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Glass::emission(Vec2 uv) const {
	return {};
}

bool Glass::is_emissive() const {
	return false;
}

bool Glass::is_specular() const {
	return true;
}

bool Glass::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Glass::display() const {
	return transmittance;
}

void Glass::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
	f(transmittance);
}

Spectrum Emissive::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Emissive::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	Scatter ret;
	ret.direction = {};
	ret.attenuation = {};
	return ret;
}

float Emissive::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Emissive::emission(Vec2 uv) const {
	return emissive.lock()->evaluate(uv);
}

bool Emissive::is_emissive() const {
	return true;
}

bool Emissive::is_specular() const {
	return true;
}

bool Emissive::is_sided() const {
	return false;
}

std::weak_ptr<Texture> Emissive::display() const {
	return emissive;
}

void Emissive::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(emissive);
}

} // namespace Materials

bool operator!=(const Materials::Lambertian& a, const Materials::Lambertian& b) {
	return a.albedo.lock() != b.albedo.lock();
}

bool operator!=(const Materials::Mirror& a, const Materials::Mirror& b) {
	return a.reflectance.lock() != b.reflectance.lock();
}

bool operator!=(const Materials::Refract& a, const Materials::Refract& b) {
	return a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Glass& a, const Materials::Glass& b) {
	return a.reflectance.lock() != b.reflectance.lock() ||
	       a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Emissive& a, const Materials::Emissive& b) {
	return a.emissive.lock() != b.emissive.lock();
}

bool operator!=(const Material& a, const Material& b) {
	if (a.material.index() != b.material.index()) return false;
	return std::visit(
		[&](const auto& material) {
			return material != std::get<std::decay_t<decltype(material)>>(b.material);
		},
		a.material);
}
