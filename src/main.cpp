#include <iostream>
#include <algorithm>
#include <string>
#include <map>
#include <chrono>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "geometry.hpp"
#include "brdf.hpp"
#include "scene.hpp"
#include "camera.hpp"
#include "integrator.hpp"
#include "utils.hpp"
#include "config.hpp"

#define clock_t std::chrono::time_point<std::chrono::system_clock>


Config conf;


int main(int argc, char* argv[])
{
	/////////////////////////////////////////////////////////////////////////
	// Parse input arguments
	/////////////////////////////////////////////////////////////////////////

	if (argc >= 2)
	{
		if (!conf.parseConfigFile(std::string(argv[1])))
		{
			return -1;
		}
		conf.printConfig();
	}
	else
	{
		std::cout << "Please specify a configuration file." << std::endl;
		return -1;
	}


	/////////////////////////////////////////////////////////////////////////
	// Camera settings
	/////////////////////////////////////////////////////////////////////////

	Eigen::Vector3f camera_pos;
	Eigen::Vector3f camera_look_at;
	if (conf.camera_id == 1)
	{
		// Camera setting 1
		camera_pos = Eigen::Vector3f(0, 0, 10);
		camera_look_at = Eigen::Vector3f(0, 0, 0);
	}
	else if (conf.camera_id == 2)
	{
		// Camera setting 2
		camera_pos = Eigen::Vector3f(5, 3, 12);
		camera_look_at = Eigen::Vector3f(0, 0, 0);
	}
	else if (conf.camera_id == 3)
	{
		// Camera setting 3
		camera_pos = Eigen::Vector3f (0, -4, 15);
		camera_look_at = Eigen::Vector3f(0, 0, 0);
	}
	float focal_len = 1.0f;
	float vertical_fov = 45.0f;
	Eigen::Vector2i film_res(conf.output_resolution, conf.output_resolution);
	Camera camera(camera_pos, focal_len, vertical_fov, film_res);
	camera.lookAt(camera_look_at, Eigen::Vector3f(0, 1, 0));
	

	/////////////////////////////////////////////////////////////////////////
	// Light setting
	/////////////////////////////////////////////////////////////////////////

	Eigen::Vector3f light_pos;
	Eigen::Vector3f light_radiance;
	Eigen::Vector2f light_size;

	if (conf.light_id == 1)
	{
		// Light setting 1
		light_pos = Eigen::Vector3f(0, 5.95, -6);
		light_radiance = Eigen::Vector3f(1.0, 1.0, 1.0) * 60;
		light_size = Eigen::Vector2f(2.0, 2.0);
	}
	else if (conf.light_id == 2)
	{
		// Light setting 2
		light_pos = Eigen::Vector3f(0, 5.95, -6);
		light_radiance = Eigen::Vector3f(1.0, 0.5, 0.3) * 40;
		light_size = Eigen::Vector2f(5.0, 4.0);
	}

	AreaLight light(light_pos, light_radiance, light_size);

	/////////////////////////////////////////////////////////////////////////
	// Scene setup
	/////////////////////////////////////////////////////////////////////////

	Scene scene(&light);
	

	/////////////////////////////////////////////////////////////////////////
	// Material settings
	/////////////////////////////////////////////////////////////////////////

	IdealDiffusion mat_front_wall(Eigen::Vector3f(0.7, 0.7, 0.7));
	IdealDiffusion mat_back_wall(Eigen::Vector3f(0.7, 0.7, 0.7));
	IdealDiffusion mat_left_wall(Eigen::Vector3f(0.5, 0, 0));
	IdealDiffusion mat_right_wall(Eigen::Vector3f(0, 0.5, 0));
	IdealDiffusion mat_floor_wall_diff(Eigen::Vector3f(0.7, 0.7, 0.7));
	IdealSpecular mat_floor_wall_spec(Eigen::Vector3f(0.7, 0.7, 0.7));
	IdealDiffusion mat_ceiling_wall(Eigen::Vector3f(0.7, 0.7, 0.7));

	IdealSpecular mat_sphere1_spec(Eigen::Vector3f(0.0, 0.2, 0.5));
	IdealDiffusion mat_sphere1_diff(Eigen::Vector3f(0.0, 0.2, 0.5));
	IdealTransmission mat_sphere1_trans(Eigen::Vector3f(1.0, 1.0, 1.0), IOR_GLASS);
	
	IdealSpecular mat_sphere2_spec(Eigen::Vector3f(0.9, 0.1, 0.4));
	IdealDiffusion mat_sphere2_diff(Eigen::Vector3f(0.4, 0.4, 0.4));

	IdealDiffusion mat_sphere3_diff(Eigen::Vector3f(0.0, 0.2, 0.5));
	IdealTransmission mat_sphere3_trans(Eigen::Vector3f(1.0, 1.0, 1.0), IOR_GLASS);
	
	IdealSubsurfaceScatter mat_sphere4_scatter(Eigen::Vector3f(0.0, 0.2, 0.5), IOR_GLASS, 0.05);

	IdealTransmission mat_bunny_trans(Eigen::Vector3f(1.0, 1.0, 1.0), IOR_GLASS);
	IdealDiffusion mat_bunny_diff(Eigen::Vector3f(0.3, 0.3, 0.5));
	IdealDiffusion mat_dragon(Eigen::Vector3f(0.2, 0.2, 0.2));

	//BRDF* mat_sphere1 = &mat_sphere1_diff;
	BRDF* mat_sphere1 = &mat_sphere1_trans;
	BRDF* mat_sphere2 = &mat_sphere2_spec;
	BRDF* mat_sphere3 = &mat_sphere3_diff;
	BRDF* mat_sphere4 = &mat_sphere4_scatter;
	BRDF* mat_bunny = &mat_bunny_trans;
	BRDF* mat_floor_wall = &mat_floor_wall_spec;
	if (!conf.enable_delta_mat)
	{
		mat_sphere2 = &mat_sphere2_diff;
		mat_bunny = &mat_bunny_diff;
		mat_floor_wall = &mat_floor_wall_diff;
	}


	/////////////////////////////////////////////////////////////////////////
	// Setting the Cornell box
	/////////////////////////////////////////////////////////////////////////

	Eigen::Vector3f cbox_lb(-6, -6, -10);
	Eigen::Vector3f cbox_size(12, 12, 30);

	Eigen::Vector3f p_p0(cbox_lb);
	Eigen::Vector3f p_s0(cbox_size.x(), 0, 0);
	Eigen::Vector3f p_s1(0, cbox_size.y(), 0);
	Eigen::Vector3f p_normal(0, 0, 1);
	Parallelogram back_wall(cbox_lb, p_s0, p_s1, p_normal, &mat_front_wall);

	p_p0 = Eigen::Vector3f(cbox_lb.x(), cbox_lb.y(), cbox_lb.z() + cbox_size.z());
	p_s0 = Eigen::Vector3f(cbox_size.x(), 0, 0);
	p_s1 = Eigen::Vector3f(0, cbox_size.y(), 0);
	p_normal = Eigen::Vector3f(0, 0, -1);
	Parallelogram front_wall(p_p0, p_s0 , p_s1, p_normal, &mat_back_wall);

	p_p0 = Eigen::Vector3f(cbox_lb);
	p_s0 = Eigen::Vector3f(0, 0, cbox_size.z());
	p_s1 = Eigen::Vector3f(0, cbox_size.y(), 0);
	p_normal = Eigen::Vector3f(1, 0, 0);
	Parallelogram left_wall(p_p0, p_s0, p_s1, p_normal, &mat_left_wall);

	p_p0 = Eigen::Vector3f(cbox_lb.x() + cbox_size.x(), cbox_lb.y(), cbox_lb.z());
	p_s0 = Eigen::Vector3f(0, 0, cbox_size.z());
	p_s1 = Eigen::Vector3f(0, cbox_size.y(), 0);
	p_normal = Eigen::Vector3f(-1, 0, 0);
	Parallelogram right_wall(p_p0, p_s0, p_s1, p_normal, &mat_right_wall);

	p_p0 = Eigen::Vector3f(cbox_lb);
	p_s0 = Eigen::Vector3f(cbox_size.x(), 0, 0);
	p_s1 = Eigen::Vector3f(0, 0, cbox_size.z());
	p_normal = Eigen::Vector3f(0, 1, 0);
	Parallelogram floor(p_p0, p_s0, p_s1, p_normal, mat_floor_wall);

	p_p0 = Eigen::Vector3f(cbox_lb.x(), cbox_lb.y() + cbox_size.y(), cbox_lb.z());
	p_s0 = Eigen::Vector3f(cbox_size.x(), 0, 0);
	p_s1 = Eigen::Vector3f(0, 0, cbox_size.z());
	p_normal = Eigen::Vector3f(0, -1, 0);
	Parallelogram ceiling(p_p0, p_s0, p_s1, p_normal, &mat_ceiling_wall);

	scene.addGeometry(&back_wall);
	scene.addGeometry(&front_wall);
	scene.addGeometry(&left_wall);
	scene.addGeometry(&right_wall);
	scene.addGeometry(&floor);
	scene.addGeometry(&ceiling);


	/////////////////////////////////////////////////////////////////////////
	// Setting geometries
	/////////////////////////////////////////////////////////////////////////

	Sphere sphere1(Eigen::Vector3f(3, 0, -5), 1.5f, mat_sphere1);
	Sphere sphere2(Eigen::Vector3f(-4, 2, -5), 1.8f, mat_sphere2);
	//Sphere sphere3(Eigen::Vector3f(-1.5, -4, -3), 1.2f, mat_sphere3);
	Sphere sphere4(Eigen::Vector3f(-2, -4, -6), 1.5f, mat_sphere4);

	//// TriangleMesh mesh_dragon("../resources/small_stanford_dragon.obj", &mat_dragon);
	//TriangleMesh mesh_dragon("../../resources/small_stanford_dragon.obj", &mat_dragon);
	//mesh_dragon.applyTransformation(Eigen::Translation3f(2.5, -7.0, -6.5) * Eigen::Scaling(25.0f));
	//mesh_dragon.buildUniformGrid();

	// TriangleMesh mesh_bunny("../resources/small_stanford_bunny.obj", mat_bunny);
	//TriangleMesh mesh_bunny("../../resources/stanford_bunny.obj", mat_bunny);
	//mesh_bunny.applyTransformation(Eigen::Translation3f(-2, -6.5, -6.0) * Eigen::Scaling(25.0f));
	//mesh_bunny.buildUniformGrid();

	scene.addGeometry(&sphere1);
	scene.addGeometry(&sphere2);
	//scene.addGeometry(&sphere3);
	scene.addGeometry(&sphere4);
	//scene.addGeometry(&mesh_dragon);
	//scene.addGeometry(&mesh_bunny);

	/////////////////////////////////////////////////////////////////////////
	// Setting  boxes
	/////////////////////////////////////////////////////////////////////////

	IdealTransmission mat_box1_trans(Eigen::Vector3f(1.0, 1.0, 1.0), IOR_WATER);
	//IdealSubsurfaceScatter mat_box1_trans(Eigen::Vector3f(1.0, 1.0, 1.0), 1);
	BRDF* mat_box1 = &mat_box1_trans;

	/*Eigen::Vector3f box1_lb(1, -3, -8);
	Eigen::Vector3f box1_size(3, 5.5, 1);*/

	Eigen::Vector3f box1_lb(2, -6, -6);
	Eigen::Vector3f box1_size(2, 5, 2);

	Eigen::Vector3f box1_p0(box1_lb);
	Eigen::Vector3f box1_s0(box1_size.x(), 0, 0);
	Eigen::Vector3f box1_s1(0, box1_size.y(), 0);
	Eigen::Vector3f box1_normal(0, 0, -1);
	Parallelogram box1_backwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	box1_p0 = Eigen::Vector3f(box1_lb[0], box1_lb[1], box1_lb[2] + box1_size[2]);
	box1_normal = Eigen::Vector3f(0, 0, 1);
	Parallelogram box1_frontwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	box1_p0 = box1_lb;
	box1_s0 = Eigen::Vector3f(0, box1_size.y(), 0);
	box1_s1 = Eigen::Vector3f(0, 0, box1_size.z());
	box1_normal = Eigen::Vector3f(-1, 0, 0);
	Parallelogram box1_leftwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	box1_p0 = Eigen::Vector3f(box1_lb[0] + box1_size[0], box1_lb[1], box1_lb[2]);
	box1_normal = Eigen::Vector3f(1, 0, 0);
	Parallelogram box1_rightwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	box1_p0 = box1_lb;
	box1_s0 = Eigen::Vector3f(box1_size.x(), 0, 0);
	box1_s1 = Eigen::Vector3f(0, 0, box1_size.z());
	box1_normal = Eigen::Vector3f(0, -1, 0);
	Parallelogram box1_bottomwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	box1_p0 = Eigen::Vector3f(box1_lb[0], box1_lb[1] + box1_size[1], box1_lb[2]);
	box1_normal = Eigen::Vector3f(0, 1, 0);
	Parallelogram box1_topwall(box1_p0, box1_s0, box1_s1, box1_normal, mat_box1);

	scene.addGeometry(&box1_backwall);
	scene.addGeometry(&box1_frontwall);
	scene.addGeometry(&box1_leftwall);
	scene.addGeometry(&box1_rightwall);
	scene.addGeometry(&box1_bottomwall);
	scene.addGeometry(&box1_topwall);

	/////////////////////////////////////////////////////////////////////////
	// Perform Phong lighting integrator
	/////////////////////////////////////////////////////////////////////////

	std::cout << "\n\nRendering..." << std::endl;

	clock_t start_time = std::chrono::system_clock::now();

	PathTracingIntegrator integrator(&scene, &camera);
	integrator.render();

	clock_t end_time = std::chrono::system_clock::now();
	double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "\nTime elapsed: " << time_elapsed << " ms" << std::endl;


	/////////////////////////////////////////////////////////////////////////
	// Output the image to a file
	/////////////////////////////////////////////////////////////////////////

	std::string output_path = conf.output_file;
	std::vector<unsigned char> output_data;
	output_data.reserve(film_res.x() * film_res.y() * 3);

	for (Eigen::Vector3f v : camera.getFilm().pixel_array)
	{
		for (int i = 0; i < 3; i++)
		{
			output_data.push_back(mathutils::gamma_correction(v[i]));
		}
	}
	stbi_flip_vertically_on_write(true);
	stbi_write_png(output_path.c_str(), film_res.x(), film_res.y(), 3, output_data.data(), 0);

	return 0;
}