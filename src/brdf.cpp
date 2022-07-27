#include <vector>
#include "Eigen/Dense"
#include "interaction.hpp"
#include "utils.hpp"
#include "constant.hpp"
#include "brdf.hpp"


/**
 * BRDF class
 */

BRDF::BRDF()
{
}


/**
 * IdealDiffusion class
 */

IdealDiffusion::IdealDiffusion(Eigen::Vector3f diff)
    : reflectivity(diff)
{
}
    
Eigen::Vector3f IdealDiffusion::eval(const Interaction& interact)
{
    return samplePdf(interact) * reflectivity;
}

float IdealDiffusion::samplePdf(const Interaction& interact)
{
    return 1.0f / (2 * PI);
}

float IdealDiffusion::sample(Interaction& interact, std::vector<float>& cdf)
{
    // Sample on a flat dist
    std::vector<float> iids = mathutils::unif(0, 1, 2);
    float r = sqrt(iids[0]);
    float theta = 2 * PI * iids[1];

    // Project it to local hemisphere
    Eigen::Vector3f wi_local = Eigen::Vector3f(r * cos(theta), r * sin(theta), 1 - pow(r, 2));

    // Rotate the local hemishpere to get global wi
    Eigen::Matrix3f rotate = Eigen::Quaternionf::FromTwoVectors(Eigen::Vector3f(0, 0, 1.0f), interact.normal).toRotationMatrix();
    interact.wi = rotate * (-wi_local);

    interact.exit_point = interact.entry_point;

    // Calculate pdf
    return samplePdf(interact);
}

bool IdealDiffusion::isDelta() const
{
    return false;
}


/**
 * IdealSpecular class
 */

IdealSpecular::IdealSpecular(Eigen::Vector3f spec)
    : reflectivity(spec)
{
}

Eigen::Vector3f IdealSpecular::eval(const Interaction& interact)
{
    return samplePdf(interact) * reflectivity;
}

float IdealSpecular::samplePdf(const Interaction& interact)
{
    Eigen::Vector3f ref = -(2 * interact.normal.dot(interact.wo) * interact.normal - interact.wo);
    if ((interact.wi - ref).norm() <= DELTA) return 1.0f;
    else return 0.0f;
}

float IdealSpecular::sample(Interaction& interact, std::vector<float>& cdf)
{
    interact.wi = -(2 * interact.normal.dot(interact.wo) * interact.normal - interact.wo);
    interact.exit_point = interact.entry_point;
    return 1.0f;
}

bool IdealSpecular::isDelta() const
{
    return true;
}


/**
 * IdealTransmission class
 */

IdealTransmission::IdealTransmission(Eigen::Vector3f reflect, float idx_refract)
    : reflectivity(reflect)
    , ior(idx_refract)
{

}

Eigen::Vector3f refract(Eigen::Vector3f v1, Eigen::Vector3f n, float ior)
{
    float cos_theta_i = v1.dot(n);
    float eta = 1 / ior;
    if (cos_theta_i < 0)
    {
        eta = 1 / eta;
        n *= -1;
        cos_theta_i *= -1;
    }
    float c = sqrt(1 - eta * eta * (1 - cos_theta_i * cos_theta_i));
    return - (c < 0 ? Eigen::Vector3f(0, 0, 0) : (eta * (-v1) + (eta * cos_theta_i - c) * n));
}

Eigen::Vector3f IdealTransmission::eval(const Interaction& interact)
{
    return samplePdf(interact) * reflectivity;
}

float IdealTransmission::samplePdf(const Interaction& interact)
{
    Eigen::Vector3f normal = interact.normal;
    float cos_alpha = normal.dot(interact.wo);
    float new_ior = ior;

    if (cos_alpha < 0) {
        normal = -interact.normal;
        cos_alpha = normal.dot(interact.wo);
        new_ior = 1 / ior;
    }

    float sin_alpha = sqrt(1 - pow(cos_alpha, 2));
    float sin_beta = sin_alpha / new_ior;
    float cos_beta = sqrt(1 - pow(sin_beta, 2));

    Eigen::Vector3f ref;
    if (cos_beta < 0) {
        ref = Eigen::Vector3f::Zero();
    }
    else {
        Eigen::Vector3f wi1 = (interact.wo - cos_alpha * normal) * (1.0f / new_ior);
        Eigen::Vector3f wi2 = cos_beta * normal;
        ref = wi1 + wi2;
    }

    if ((interact.wi - ref).norm() <= DELTA) return 1.0f;
    else return 0.0f;
}

float IdealTransmission::sample(Interaction& interact, std::vector<float>& cdf)
{
    interact.exit_point = interact.entry_point;

    Eigen::Vector3f normal = interact.normal;
    float cos_alpha = normal.dot(interact.wo);
    float new_ior = ior;

    if (cos_alpha < 0) {
        normal = -interact.normal;
        cos_alpha = normal.dot(interact.wo);
        new_ior = 1 / ior;
    }

    float sin_alpha = sqrt(1 - pow(cos_alpha, 2));
    float sin_beta = sin_alpha / new_ior;
    float cos_beta_squ = 1 - pow(sin_beta, 2);

    if (cos_beta_squ < 0) {
        interact.wi = Eigen::Vector3f::Zero();
        return 1.0f;
    }

    float cos_beta = sqrt(1 - pow(sin_beta, 2));
    Eigen::Vector3f wi1 = (interact.wo - cos_alpha * normal) * (1.0f / new_ior);
    Eigen::Vector3f wi2 = cos_beta * normal;
    interact.wi = wi1 + wi2;

    return 1.0f;
}

bool IdealTransmission::isDelta() const
{
    return true;
}

/**
 * IdealSubsurfaceScatter class
 */

IdealSubsurfaceScatter::IdealSubsurfaceScatter(Eigen::Vector3f reflec, float idx_refract, float scatter_d)
    : reflectivity(reflec)
    , ior(idx_refract)
    , scatter_dist(scatter_d)
{
}

// Performs sampling of a Normalized Burley diffusion profile in polar coordinates.
// 'u' is the random number (the value of the CDF): [0, 1).
// rcp(s) = 1 / ShapeParam = ScatteringDistance.
// 'r' is the sampled radial distance, s.t. (u = 0 -> r = 0) and (u = 1 -> r = Inf).
// rcp(Pdf) is the reciprocal of the corresponding PDF value.
void IdealSubsurfaceScatter::sampleBurleyDiffusionProfile(float u, float rcpS, float& r, float& rcpPdf)
{
    u = 1 - u; // Convert CDF to CCDF; the resulting value of (u != 0)

    float g = 1 + (4 * u) * (2 * u + sqrt(1 + (4 * u) * u));
    float n = exp2(log2(g) * (-1.0 / 3.0));                    // g^(-1/3)
    float p = (g * n) * n;                                     // g^(+1/3)
    float c = 1 + p + n;                                       // 1 + g^(+1/3) + g^(-1/3)
    float x = (3 / log2(2.7182818)) * log2(c / (4 * u));       // 3 * Log[c / (4 * u)]

    // x      = s * r
    // exp_13 = Exp[-x/3] = Exp[-1/3 * 3 * Log[c / (4 * u)]]
    // exp_13 = Exp[-Log[c / (4 * u)]] = (4 * u) / c
    // exp_1  = Exp[-x] = exp_13 * exp_13 * exp_13
    // expSum = exp_1 + exp_13 = exp_13 * (1 + exp_13 * exp_13)
    // rcpExp = rcp(expSum) = c^3 / ((4 * u) * (c^2 + 16 * u^2))
    float rcpExp = ((c * c) * c) / ((4 * u) * ((c * c) + (4 * u) * (4 * u)));

    r = x * rcpS;
    rcpPdf = (8 * PI * rcpS) * rcpExp; // (8 * Pi) / s / (Exp[-s * r / 3] + Exp[-s * r])
}

float IdealSubsurfaceScatter::calFresnel(Eigen::Vector3f normal, Eigen::Vector3f w)
{
    float F_w = 0;

    float cos_w = w.dot(normal);
    float sin_w = sqrt(1 - pow(cos_w, 2));
    float sin_wt = sin_w / ior;
    if (1 - pow(sin_wt, 2) < 0) {
        F_w = 1;
    }
    else {
        float cos_wt = sqrt(1 - pow(sin_wt, 2));
        float Rs = pow((cos_w - ior * cos_wt) / (cos_w + ior * cos_wt), 2);
        float Rp = pow((cos_wt - ior * cos_w) / (cos_wt + ior * cos_w), 2);
        float F_wi = 0.5 * (Rs + Rp);
    }

    return F_w;
}

Eigen::Vector3f IdealSubsurfaceScatter::eval(const Interaction& interact)
{
    float F_wi = 1 - calFresnel(interact.normal_wi, -interact.wi);
    float F_wo = 1 - calFresnel(interact.normal, interact.wo);

    Eigen::Vector3f result = F_wi * F_wo * samplePdf(interact) * reflectivity / PI;
    return result;
}

float IdealSubsurfaceScatter::samplePdf(const Interaction& interact)
{
    /*float pdf = (exp(-(interact.r / this->scatter_dist)) + exp(-(interact.r / (3.0f * this->scatter_dist)))) / (8 * PI * interact.r * this->scatter_dist) / (2 * PI);*/
    float pdf = (exp(-(interact.r / this->scatter_dist)) + exp(-(interact.r / (3.0f * this->scatter_dist)))) / (8 * PI * interact.r * this->scatter_dist);
    return pdf;
}

float IdealSubsurfaceScatter::sample(Interaction& interact, std::vector<float>& cdf)
{
    // Sample radius
    float u = mathutils::unif(0, 0.36, 1)[0];
    for (int i = 0; i < 100; i++) {
        if (u < cdf[i]) {
            interact.r = i / 100.0f;
            break;
        }
    }
    float pdf = (exp(-(interact.r / 1.0f)) + exp(-(interact.r / 3.0f))) / (8 * PI * interact.r);
    //std::cout << interact.r << " " << pdf << "\n";

    // Sample angle
    float angle = mathutils::unif(0, 2 * PI, 1)[0];

    // Find exit point on plate
    Eigen::Vector3f exit_point_local = Eigen::Vector3f(interact.r * cos(angle), interact.r * sin(angle), 0);
    Eigen::Matrix3f rotate_matrix = Eigen::Quaternionf::FromTwoVectors(Eigen::Vector3f(0, 0, 1.0f), interact.normal).toRotationMatrix();
    Eigen::Vector3f exit_point_plate = rotate_matrix * exit_point_local + interact.entry_point;

    // Project exit point to a point on the surface
    Eigen::Vector3f sphere_center = Eigen::Vector3f(-2, -4, -6);
    float sphere_radius = 1.5f;
    interact.normal_wi = (exit_point_plate - sphere_center).normalized();
    Eigen::Vector3f proj_vec = (sphere_center - exit_point_plate) + interact.normal_wi * sphere_radius;
    interact.exit_point = exit_point_plate + proj_vec;
    if ((interact.exit_point - sphere_center).norm() - sphere_radius > DELTA) {
        std::cout << "Norm: " << (interact.exit_point - sphere_center).norm() << "   ";
    }

    // Sample wi
    // Sample on a flat dist
    std::vector<float> iids = mathutils::unif(0, 1, 2);
    float r = sqrt(iids[0]);
    float theta = 2 * PI * iids[1];
    // Project it to local hemisphere
    Eigen::Vector3f wi_local = Eigen::Vector3f(r * cos(theta), r * sin(theta), 1 - pow(r, 2));
    // Rotate the local hemishpere to get global wi
    interact.wi = rotate_matrix * (-wi_local);

    //return pdf / (2.0f * PI);
    return pdf;


    //// Sample radius
    //float u = mathutils::unif(0, 1, 1)[0];
    //float rcpS = 1 / scatter_dist;
    //float rcpPdf;
    //sampleBurleyDiffusionProfile(u, rcpS, interact.r, rcpPdf);

    //// Sample angle
    //float angle = mathutils::unif(0, 2 * PI, 1)[0];
    //
    //// Find exit point on plate
    //Eigen::Vector3f exit_point_local = Eigen::Vector3f(interact.r * cos(angle), interact.r * sin(angle), 0);
    //Eigen::Matrix3f rotate_matrix = Eigen::Quaternionf::FromTwoVectors(Eigen::Vector3f(0, 0, 1.0f), interact.normal).toRotationMatrix();
    //Eigen::Vector3f exit_point_plate = rotate_matrix * exit_point_local + interact.entry_point;

    //// Project exit point to a point on the surface
    //Eigen::Vector3f sphere_center = Eigen::Vector3f(-2, -4, -6);
    //float sphere_radius = 1.5f;
    //Eigen::Vector3f proj_dir = (exit_point_plate - sphere_center).normalized();
    //Eigen::Vector3f proj_dist = (exit_point_plate - sphere_center) - proj_dir * sphere_radius;
    //interact.exit_point = exit_point_plate + proj_dist;

    //// Sample wi
    //// Sample on a flat dist
    //std::vector<float> iids = mathutils::unif(0, 1, 2);
    //float r = sqrt(iids[0]);
    //float theta = 2 * PI * iids[1];
    //// Project it to local hemisphere
    //Eigen::Vector3f wi_local = Eigen::Vector3f(r * cos(theta), r * sin(theta), 1 - pow(r, 2));
    //// Rotate the local hemishpere to get global wi
    //Eigen::Matrix3f rotate = Eigen::Quaternionf::FromTwoVectors(Eigen::Vector3f(0, 0, 1.0f), interact.normal).toRotationMatrix();
    //interact.wi = rotate * (-wi_local);

    //return rcpPdf / (2.0f * PI);
}

bool IdealSubsurfaceScatter::isDelta() const
{
    return true;
}