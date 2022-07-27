#ifndef NO_OMP
#include <omp.h>
#endif
#include "progressbar.hpp"
#include "integrator.hpp"
#include "constant.hpp"
#include "light.hpp"
#include "utils.hpp"
#include "config.hpp"
#include "brdf.hpp"



extern Config conf;


//std::vector<float> cdf;

/**
 * Integrator class
 */

Integrator::Integrator(Scene* scn, Camera* cam)
    : scene(scn)
    , camera(cam)
{
}


/**
 * PathTracingIntegrator class
 */

PathTracingIntegrator::PathTracingIntegrator(Scene* scene, Camera* camera)
    : Integrator(scene, camera)
{  
}

void precompute_cdf_inv(std::vector<float>& cdf)
{
    for (int r = 0; r < 100; r++) {
        cdf[r] = 1 - 0.25 * exp(-(r / 100.0f)) - 0.75 * exp(-(r / 300.0f));
        //std::cout << r << " " << cdf[r] << "\n";
    }
}

void PathTracingIntegrator::render()
{
    int dx, dy;
    int res_x = camera->getFilm().resolution.x(), res_y = camera->getFilm().resolution.y();

    /* Initialize a progress bar */
    progressbar progress_bar(res_x * res_y);

    /* Precompute CDF^-1 */
    std::vector<float> cdf(100, 0);
    precompute_cdf_inv(cdf);

#ifndef NO_OMP
    #pragma omp parallel for private(dy)
#endif
    for (dx = 0; dx < res_x; ++dx)
    {
        for (dy = 0; dy < res_y; ++dy)
        {
            //if (dx == 200 && dy == 140) {
            //    Eigen::Vector3f color(0, 0, 0);
            //    for (int j = 0; j < conf.spp; j++) {
            //        Ray ray = camera->generateRay(dx, dy);
            //        color += radiance(ray, conf.max_depth, cdf);
            //    }
            //    camera->setPixel(dx, dy, color / conf.spp);

            //    /*Eigen::Vector3f color(1, 1, 1);
            //    camera->setPixel(dx, dy, color / conf.spp);*/

            //    //camera->setPixel(dx, dy, Eigen::Vector3f(1.0f, 1.0f, 1.0f));
            //}
            //else {
            //    Eigen::Vector3f color(0, 0, 0);
            //    camera->setPixel(dx, dy, color / conf.spp);
            //}
            
            Eigen::Vector3f color(0, 0, 0);
            for (int j = 0; j < conf.spp; j++) {
                std::vector<float> random_xy = mathutils::unif(0, 1, 2);
                Ray ray = camera->generateRay(dx + random_xy[0], dy + random_xy[1]);
                color += radiance(ray, conf.max_depth, cdf);
            }
            camera->setPixel(dx, dy, color / conf.spp);


#ifndef NO_OMP
            //#pragma omp critical
#endif
            //progress_bar.update();
        }
        std::cout << dx << " ";
    }
    std::cout << "done.\n";
}

Eigen::Vector3f PathTracingIntegrator::radiance(Ray ray, int depth, std::vector<float>& cdf)
{
    if (depth == 0) return Eigen::Vector3f::Zero();

    // init color for current pixel
    Eigen::Vector3f dir_light = Eigen::Vector3f::Zero();
    Eigen::Vector3f indir_light = Eigen::Vector3f::Zero();

    Interaction inter;
    if (!scene->intersection(ray, inter)) {
        return Eigen::Vector3f::Zero();
    }
    inter.wo = -ray.direction;

    if (inter.type == Interaction::Type::NONE) {
        return Eigen::Vector3f::Zero();
    }
    else if (inter.type == Interaction::Type::LIGHT) {
        return scene->getLight()->emission(inter.entry_point, inter.wo);
    }
    else {
        // Sample a point from area light and shoot a ray
        float pdf_light = 0.0f;
        Eigen::Vector3f sampled_pos_light = scene->getLight()->sample(inter, &pdf_light);
        Ray light_ray(inter.entry_point, sampled_pos_light - inter.entry_point, DELTA);
        inter.wi = -light_ray.direction;

        // Add contribution from light if not shadowed
        if (!scene->isShadowed(light_ray)) {
            float dist = (sampled_pos_light - inter.entry_point).norm();
            Eigen::Vector3f Li_cos = scene->getLight()->emission(sampled_pos_light, -light_ray.direction);

            float cos_theta = inter.normal.dot(-inter.wi);

            BRDF* inter_brdf = (BRDF*) (inter.material);
            Eigen::Vector3f fr = inter_brdf->eval(inter);

            dir_light = Li_cos.cwiseProduct(fr) * cos_theta / pow(dist, 2) / pdf_light;
        }

        // Sample a new ray randomly
        BRDF* inter_brdf = (BRDF*)(inter.material);
        float pdf_wi = inter_brdf->sample(inter, cdf);
        Eigen::Vector3f fr = inter_brdf->eval(inter);
        Ray new_ray(inter.exit_point, -inter.wi, DELTA);

        Interaction new_inter;
        if (scene->intersection(new_ray, new_inter)) {
            if (new_inter.type != Interaction::Type::NONE && new_inter.type != Interaction::Type::LIGHT) {
                Eigen::Vector3f Li_cos = radiance(new_ray, depth - 1, cdf);
                // Abs is for transmission
                float cos_theta = abs(inter.normal.dot(-inter.wi));
                indir_light = Li_cos.cwiseProduct(fr) * cos_theta / pdf_wi;
            }
        }
    }
    //std::cout << "L: " << dir_light(0) << " " << dir_light(1) << " " << dir_light(2) << "   " << indir_light(0) << " " << indir_light(1) << " " << indir_light(2) << "\n";
    return dir_light + indir_light;
}

//Eigen::Vector3f PathTracingIntegrator::radiance(Ray ray, int depth)
//{
//    if (depth == 0) return Eigen::Vector3f::Zero();
//
//    // init color for current pixel
//    Eigen::Vector3f dir_light = Eigen::Vector3f::Zero();
//    Eigen::Vector3f indir_light = Eigen::Vector3f::Zero();
//
//    Interaction inter;
//    if (!scene->intersection(ray, inter)) {
//        return Eigen::Vector3f::Zero();
//    }
//    inter.wo = -ray.direction;
//
//    if (inter.type == Interaction::Type::NONE) {
//        return Eigen::Vector3f::Zero();
//    }
//    else if (inter.type == Interaction::Type::LIGHT) {
//        return scene->getLight()->emission(inter.entry_point, inter.wo);
//    }
//    else {
//        // Sample a point from area light
//        float pdf_light = 0.0f;
//        Eigen::Vector3f sampled_pos_light = scene->getLight()->sample(inter, &pdf_light);
//        Ray light_ray(inter.entry_point, sampled_pos_light - inter.entry_point, DELTA);
//        inter.wi = -light_ray.direction;
//
//        // Add contribution from light if not shadowed
//        if (!scene->isShadowed(light_ray)) {
//            float dist = (sampled_pos_light - inter.entry_point).norm();
//            Eigen::Vector3f Li_cos = scene->getLight()->emission(sampled_pos_light, -light_ray.direction);
//
//            float cos_theta = inter.normal.dot(-inter.wi);
//
//            BRDF* inter_brdf = (BRDF*)(inter.material);
//            Eigen::Vector3f fr = inter_brdf->eval(inter);
//
//            dir_light = Li_cos.cwiseProduct(fr) * cos_theta / pow(dist, 2) / pdf_light;
//        }
//
//        // Sample a new ray randomly
//        BRDF* inter_brdf = (BRDF*)(inter.material);
//        float pdf_wi = inter_brdf->sample(inter);
//        Eigen::Vector3f fr = inter_brdf->eval(inter);
//        Ray new_ray(inter.entry_point, -inter.wi, DELTA);
//
//        Interaction new_inter;
//        if (scene->intersection(new_ray, new_inter)) {
//            if (new_inter.type != Interaction::Type::NONE && new_inter.type != Interaction::Type::LIGHT) {
//                Eigen::Vector3f Li_cos = radiance(new_ray, depth - 1);
//                // Abs is for transmission
//                float cos_theta = abs(inter.normal.dot(-inter.wi));
//                indir_light = Li_cos.cwiseProduct(fr) * cos_theta / pdf_wi;
//            }
//        }
//    }
//    return dir_light + indir_light;
//}

Eigen::Vector3f PathTracingIntegrator::directLightning(Eigen::Vector3f inter_pos, IdealDiffusion * idealDiffusion)
{
    AreaLight * light = (AreaLight *) scene->getLight();

    Interaction light_it;
    float pdf;
    Eigen::Vector3f light_pos = light->sample(light_it, & pdf);

    Ray light_ray(inter_pos, (light_pos - inter_pos).normalized(), 1e-3);
    
    if (scene->isShadowed(light_ray)) {
        return Eigen::Vector3f(0, 0, 0);
    } else {
        return light->emission(light_pos, light_it.wo);
    }
}