#pragma once
#include "Eigen/Dense"


/**
 * The Phong lighting model at the interacting point
 */

struct InteractionPhongModel
{
    Eigen::Vector3f diffusion;
    Eigen::Vector3f specular;
	Eigen::Vector3f transmission;
	Eigen::Vector3f subsurface_scatter;
    float shininess;
};


/**
 * Data structure representing interaction between objects and rays
 */

struct Interaction
{
	enum Type
	{
		NONE,
		GEOMETRY,
		LIGHT
	};

	/* Distance (in units of t) to intersection point */
	float entry_dist;
	/* Distance (in units of t) to the second intersection point(if existed) */
	float exit_dist;
	/* Position of intersection point */
	Eigen::Vector3f entry_point;
	/* Position of ray exit point, only used by subsurface scattering */
	Eigen::Vector3f exit_point;
	/* Normal of intersection point (if existed) */
	Eigen::Vector3f normal;
	/* Material at the intersected point (if existed) */
	void* material;
	/* Direction of incoming radiance */
	Eigen::Vector3f wi;
	/* Direction of outcoming radiance */
	Eigen::Vector3f wo;
	/* Type of interacting object */
	Type type;
	/* Whether it is intersecting with an bssrdf object */
	bool inter_bssrdf;
	/* Bssrdf: sampled dist */
	float r;
	/* Bssrdf: normal of intersection point */
	Eigen::Vector3f normal_wi;

	// IdealDiffusion brdf;

	Interaction() : entry_dist(-1), exit_dist(-1), material(nullptr), type(Type::NONE), inter_bssrdf(false), r(0.2) {}
};
