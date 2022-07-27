#pragma once
#include "scene.hpp"
#include "camera.hpp"
#include "interaction.hpp"


/**
 * Base class of integrator
 */

class Integrator
{
protected:
	Scene* scene;
	Camera* camera;

public:
	Integrator(Scene* scn, Camera* cam);
	virtual ~Integrator() = default;
	
	virtual void render() = 0;
	virtual Eigen::Vector3f radiance(Ray ray, int depth, std::vector<float>& cdf) = 0;
};


/**
 * Path-tracing integrator
 */

class PathTracingIntegrator : public Integrator
{
public:
    PathTracingIntegrator(Scene* scene, Camera* camera);
    virtual void render() override;
    virtual Eigen::Vector3f radiance(Ray ray, int depth, std::vector<float>& cdf) override;
	Eigen::Vector3f directLightning(Eigen::Vector3f inter_pos, IdealDiffusion * idealDiffusion);
};	
