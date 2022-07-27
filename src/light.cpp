#include <utility>
#include "light.hpp"
#include "geometry.hpp"
#include "utils.hpp"


/**
 * Light class
 */

Light::Light(Eigen::Vector3f pos, Eigen::Vector3f power)
    : position(pos)
    , radiance(power)
{
}

Eigen::Vector3f Light::getPosition() const
{
    return position;
}

Eigen::Vector3f Light::getRadiance() const
{
    return radiance;
}


/**
 * AreaLight class
 */
AreaLight::AreaLight(Eigen::Vector3f pos, Eigen::Vector3f power, Eigen::Vector2f size)
    : Light(pos, power)
    , area_size(size)
    , geometry_delegation(
        pos - Eigen::Vector3f(size[0], 0, size[1]) / 2,
        Eigen::Vector3f(size[0], 0, 0),
        Eigen::Vector3f(0, 0, size[1]),
        Eigen::Vector3f(0, -1, 0), nullptr)
{
}

Eigen::Vector3f AreaLight::emission(Eigen::Vector3f pos, Eigen::Vector3f dir)
{
    return std::max(dir.dot(geometry_delegation.getNormal()), 0.0f) * radiance;
}

float AreaLight::samplePdf(const Interaction& ref_it, Eigen::Vector3f pos)
{
    return (1.0f / (area_size(0) * area_size(1)));
}

Eigen::Vector3f AreaLight::sample(Interaction& ref_it, float* pdf)
{
    // By the design of framework, arealight is always aligned with xyz coordinate, 
    // and points to y direction (), so the following operation is safe
    Eigen::Vector2f center(position[0], position[2]);

    Eigen::Vector2f ori_pos = center - (area_size / 2.0f);
    Eigen::Vector2f des_pos = center + (area_size / 2.0f);

    Eigen::Vector3f sampled_pos;
    sampled_pos(0) = mathutils::unif(ori_pos(0), des_pos(0), 1)[0];
    sampled_pos(1) = position[1];
    sampled_pos(2) = mathutils::unif(ori_pos(1), des_pos(1), 1)[0];

    *pdf = samplePdf(ref_it, sampled_pos);
    return sampled_pos;
}

bool AreaLight::rayIntersection(Interaction& interaction, const Ray& ray)
{
    bool intersection = geometry_delegation.rayIntersection(interaction, ray);
    interaction.type = Interaction::Type::LIGHT;
    return intersection;
}
