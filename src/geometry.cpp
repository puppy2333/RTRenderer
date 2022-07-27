#include <math.h>
#include "geometry.hpp"
#include "constant.hpp"
#include "utils.hpp"
#include "objloader.hpp"

#include <Eigen/Core>

/**
 * Geometry class
 */

void Geometry::setMaterial(BRDF* new_mat)
{
    material = new_mat;
}


/**
 * Parallelogram class
 */

Parallelogram::Parallelogram(Eigen::Vector3f p0, Eigen::Vector3f s0, Eigen::Vector3f s1, Eigen::Vector3f normal, BRDF* mat)
    : p0(p0)
    , normal(normal.normalized())
{
    s0_len = s0.norm();
    s1_len = s1.norm();
    this->s0 = s0.normalized();
    this->s1 = s1.normalized();

    setMaterial(mat);
    buildBoundingBox();
}

bool Parallelogram::rayIntersection(Interaction& interaction, const Ray& ray)
{
    if (ray.direction.dot(normal) == 0)
    {
        return false;
    }
    
    float t = (p0 - ray.origin).dot(normal) / ray.direction.dot(normal);
    Eigen::Vector3f p0_p = ray.getPoint(t) - p0;
    float q0 = p0_p.dot(s0) / s0_len;
    float q1 = p0_p.dot(s1) / s1_len;
    if (q0 >= 0 && q0 <= s0.norm() && q1 >= 0 && q1 <= s1.norm() && t >= ray.range_min && t <= ray.range_max)
    {
        interaction.entry_dist = t;
        interaction.exit_dist = t;
        interaction.normal = normal;
        interaction.entry_point = ray.getPoint(t);
        interaction.wo = -ray.direction;
        // interaction.brdf = material;
        if (material != nullptr)
        {
            interaction.material = (void*)material;
        }
        interaction.type = Interaction::Type::GEOMETRY;
        return true;
    }
    return false;
}

void Parallelogram::buildBoundingBox()
{
    Eigen::Vector3f p1 = p0 + s0 + s1;
    bounding_box.lb = p0.cwiseMin(p1);
    bounding_box.ub = p0.cwiseMax(p1);
}

Eigen::Vector3f Parallelogram::getOrigin()
{
    return p0;
}

std::pair<Eigen::Vector3f, Eigen::Vector3f> Parallelogram::getDirections()
{
    return std::pair<Eigen::Vector3f, Eigen::Vector3f>(s0, s1);
}

Eigen::Vector3f Parallelogram::getNormal()
{
    return normal;
}

/**
 * Sphere class
 */

Sphere::Sphere(Eigen::Vector3f p0, float r, BRDF* mat)
    : p0(p0)
    , radius(r)
{
    setMaterial(mat);
    buildBoundingBox();
}

bool Sphere::rayIntersection(Interaction& interaction, const Ray& ray)
{
    float a = 1.0f;
    float b = 2 * ray.direction.dot(ray.origin - p0);
    float c = (ray.origin - p0).squaredNorm() - radius * radius;
    float delta = b * b - 4 * a * c;

    if (delta < 0)
    {
        return false;
    }

    float t0 = (-b - sqrt(delta)) / 2 * a;
    float t1 = (-b + sqrt(delta)) / 2 * a;
    
    if (t1 < 0)
    {
        return false;
    }
    /*else if (t0 < 0 && t1 >= 0)
    {
        t0 = t1;
    }*/
    else if (t0 < ray.range_min && t1 >= ray.range_min)
    {
        t0 = t1;
    }

    if (t0 < ray.range_min || t0 > ray.range_max)
    {
        return false;
    }

    interaction.entry_dist = t0;
    interaction.exit_dist = t1;
    interaction.entry_point = ray.getPoint(t0);
    Eigen::Vector3f r_vec = interaction.entry_point - p0;
    interaction.normal = r_vec.normalized();
    interaction.wo = -ray.direction;
    if (material != nullptr)
    {
        interaction.material = (void*)material;
        //std::cout << typeid(material).name() << "   ";
        /*if (typeid(material).name() == "IdealSubsurfaceScatter") {
            interaction.inter_bssrdf = true;
        }*/
        
    }
    interaction.type = Interaction::Type::GEOMETRY;

    return true;
}

void Sphere::buildBoundingBox()
{
    bounding_box = AABB(p0, radius);
}


/**
 * TriangleMesh class
 */

TriangleMesh::TriangleMesh(std::string file_path, BRDF* mat)
{
    setMaterial(mat);

    std::vector<Eigen::Vector2f> out_uvs;
    std::vector<int> out_vt_index;
    loadObj(file_path.c_str(), vertices, out_uvs, normals, vertices_index, out_vt_index, normals_index);
    
    num_triangles = vertices_index.size() / 3;

    has_grid = false;

    buildBoundingBox();
}

bool TriangleMesh::raySingleTriangleIntersection(Interaction& interaction, const Ray& ray, int v0_idx, int v1_idx, int v2_idx) const
{
    /**
     * TODO: Ray intersection test with single triangle
     * Note: Remember that normals are interpolated using barycentric coordinates.
     */
    Eigen::Vector3f P0 = vertices[vertices_index[v0_idx]];
    Eigen::Vector3f P1 = vertices[vertices_index[v1_idx]];
    Eigen::Vector3f P2 = vertices[vertices_index[v2_idx]];

    Eigen::Vector3f E1 = P1 - P0;
    Eigen::Vector3f E2 = P2 - P0;
    Eigen::Vector3f O  = ray.origin;
    Eigen::Vector3f S  = O - P0;
    Eigen::Vector3f D  = ray.direction;
    Eigen::Vector3f S1 = D.cross(E2);
    Eigen::Vector3f S2 = S.cross(E1);

    float t  = S2.dot(E2) / E1.dot(S1);
    float b1 = S1.dot(S)  / E1.dot(S1);
    float b2 = S2.dot(D)  / E1.dot(S1);

    if (t > 0 && b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && (b1 + b2) <= 1) {
        interaction.entry_dist = t;
        interaction.exit_dist = t;
        interaction.entry_point = ray.getPoint(t);
        interaction.normal = E1.cross(E2).normalized();
        interaction.wo = -ray.direction;
        return true;
    } else {
        return false;
    }

    UNREACHABLE;
}

bool TriangleMesh::rayIntersection(Interaction& interaction, const Ray& ray)
{
    // The closest intersection
    Interaction final_interaction;

    if (has_grid) {

        /**
         * TODO: Use uniform grid to handle triangle intersection here
         * Note: Grid traversal algorithm must be used here.
         */

    
        // Find the intersection point with the whole bounding box
        float t_in, t_out;
        if (bounding_box.rayIntersection(ray, t_in, t_out)) {

            Eigen::Vector3f inter_point = ray.getPoint(t_in);

            // std::cout << "Points: " << inter_point(0) << " " << inter_point(1) << " " << inter_point(2) << " \n"; 
            // std::cout << "Bounding box lb: " << bounding_box.lb(0) << " " << bounding_box.lb(1) << " " << bounding_box.lb(2) << " \n"; 
            // std::cout << "Bounding box ub: " << bounding_box.ub(0) << " " << bounding_box.ub(1) << " " << bounding_box.ub(2) << " \n"; 

            // Current pos of the ray
            float t = 0;

            Eigen::Vector3f dir = ray.direction;
            for (int i = 0; i < 3; i++) {
                if (dir(0) == 0) dir(0) = 1e-5;
                if (dir(1) == 0) dir(1) = 1e-5;
                if (dir(2) == 0) dir(2) = 1e-5;
            }

            // Get the ratio
            float kx = 1.0f / dir(0);
            float ky = 1.0f / dir(1);
            float kz = 1.0f / dir(2);

            if (abs(inter_point(0) - bounding_box.lb(0)) <= 1e-6) inter_point(0) += 0.00001;
            if (abs(inter_point(1) - bounding_box.lb(1)) <= 1e-6) inter_point(1) += 0.00001;
            if (abs(inter_point(2) - bounding_box.lb(2)) <= 1e-6) inter_point(2) += 0.00001;

            // std::cout << "kx ky kz: " << kx << " " << ky << " " << kz << "\n";

            int idx_x = (int) floor(inter_point(0) - bounding_box.lb(0));
            int idx_y = (int) floor(inter_point(1) - bounding_box.lb(1));
            int idx_z = (int) floor(inter_point(2) - bounding_box.lb(2));

            // std::cout << "\n" << inter_point(0) << " " << inter_point(1) << " " << inter_point(2) << " "; 
            // std::cout << bounding_box.lb(0) << " " << bounding_box.lb(1) << " " << bounding_box.lb(2) << " \n"; 
            // std::cout << "\n" << idx_x << " " << idx_y << " " << idx_z << " \n"; 

            // float dx = ceil(inter_point(0)) - inter_point(0);
            // float dy = ceil(inter_point(1)) - inter_point(1);
            // float dz = ceil(inter_point(2)) - inter_point(2);

            // float tx = (ceil(inter_point(0)) - inter_point(0)) * abs(kx);
            // float ty = (ceil(inter_point(1)) - inter_point(1)) * abs(ky);
            // float tz = (ceil(inter_point(2)) - inter_point(2)) * abs(kz);

            float dx, dy, dz;
            if (kx >= 0) {
                dx = bounding_box.lb(0) + idx_x + 1 - inter_point(0);     
            } else {
                dx = inter_point(0) - bounding_box.lb(0) - idx_x;
            }
            if (ky >= 0) {
                dy = bounding_box.lb(1) + idx_y + 1 - inter_point(1);
            } else {
                dy = inter_point(1) - bounding_box.lb(1) - idx_y;
            }
            if (kz >= 0) {
                dz = bounding_box.lb(2) + idx_z + 1 - inter_point(2);     
            } else {
                dz = inter_point(2) - bounding_box.lb(2) - idx_z;
            }

            // float dx = bounding_box.lb(0) + idx_x + 1 - inter_point(0);
            // float dy = bounding_box.lb(1) + idx_y + 1 - inter_point(1);
            // float dz = bounding_box.lb(2) + idx_z + 1 - inter_point(2);

            float tx = dx * abs(kx);
            float ty = dy * abs(ky);
            float tz = dz * abs(kz);

            // std::cout << "init dx dy dz: " << dx << " " << dy << " " << dz << "\n";
            // std::cout << "init tx ty tz: " << tx << " " << ty << " " << tz << "\n";

            int step_x, step_y, step_z;
            if (kx >= 0) step_x = 1;
            else        step_x = -1;
            if (ky >= 0) step_y = 1;
            else        step_y = -1;
            if (kz >= 0) step_z = 1;
            else        step_z = -1;


            bool intersected = false;
            while (intersected == false) {
                
                // std::cout << "\n" << inter_point(0) << " " << inter_point(1) << " " << inter_point(2) << " ";
                // std::cout << bounding_box.lb(0) << " " << bounding_box.lb(1) << " " << bounding_box.lb(2) << " \n";
                // std::cout << idx_x << " " << idx_y << " " << idx_z << " \n"; 

                // Check whether the current box has triangles
                Interaction cur_it;
                for (int i = 0; i < big_bounding_box.cellsList[idx_x][idx_y][idx_z].size(); i++) {
                    
                    int tri_idx1 = big_bounding_box.cellsList[idx_x][idx_y][idx_z][i](0);
                    int tri_idx2 = big_bounding_box.cellsList[idx_x][idx_y][idx_z][i](1);
                    int tri_idx3 = big_bounding_box.cellsList[idx_x][idx_y][idx_z][i](2);

                    // std::cout << "\n" << tri_idx1 << " " << tri_idx2 << " " << tri_idx3 << " \n"; 

                    if (raySingleTriangleIntersection(cur_it, ray, tri_idx1, tri_idx2, tri_idx3)) {
                        intersected = true;
                        // std::cout << "hhh\n";
                        if (final_interaction.entry_dist == -1 || cur_it.entry_dist < final_interaction.entry_dist) {
                            final_interaction = cur_it;
                        }
                    }
                }
                
                if (tx < ty && tx < tz) {
                    t = tx;
                    tx += abs(kx);
                    idx_x += step_x;
                } else if (ty < tz && ty < tx) {
                    t = ty;
                    ty += abs(ky);
                    idx_y += step_y;
                } else {
                    t = tz;
                    tz += abs(kz);
                    idx_z += step_z;
                }

                // std::cout << "\n" << idx_x << " " << idx_y << " " << idx_z << " \n"; 

                // std::cout << "t: " << t_out - t_in << " " << t << "\n";

                // if (t >= t_out - t_in) break;

                if (idx_x < 0 || idx_x >= big_bounding_box.dim(0) || 
                    idx_y < 0 || idx_y >= big_bounding_box.dim(1) || 
                    idx_z < 0 || idx_z >= big_bounding_box.dim(2)) 
                {
                    break;
                }
            }
        }

        
        // UNREACHABLE;

    } 

    else {
        for (int i = 0; i < num_triangles; i++) {
            Interaction cur_it;
            if (raySingleTriangleIntersection(cur_it, ray, 3 * i, 3 * i + 1, 3 * i + 2)) {
                if (final_interaction.entry_dist == -1 || cur_it.entry_dist < final_interaction.entry_dist) {
                    final_interaction = cur_it;
                }
            }
        }
    }

    // Meaning "if intersected == true"
    if (final_interaction.entry_dist != -1)
    {
        interaction = final_interaction;

        // Material is an attribute bonding to the geometry, now is asigned to the "interaction"
        if (material != nullptr)
        {
            interaction.material = material;
        }
        interaction.type = Interaction::Type::GEOMETRY;

        return true;
    }

    // If not intersected with this whole triangle mesh
    return false;
}

// Get a bounding box which covers all the triangles in the mesh
void TriangleMesh::buildBoundingBox()
{

    bounding_box.lb = vertices[0].cwiseMin(vertices[1]);
    bounding_box.ub = vertices[0].cwiseMax(vertices[1]);
    for (int i = 2; i < vertices.size(); i++)
    {
        bounding_box.lb = bounding_box.lb.cwiseMin(vertices[i]);
        bounding_box.ub = bounding_box.ub.cwiseMax(vertices[i]);
    }
}

// Transform the current mesh to a new position
void TriangleMesh::applyTransformation(const Eigen::Affine3f& t)
{
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] = t * vertices[i];
    }

    Eigen::Matrix3f t_inv_tr = t.linear().inverse().transpose();
    for (int i = 0; i < normals.size(); i++)
    {
        normals[i] = (t_inv_tr * normals[i]).normalized();
    }

    buildBoundingBox();
}

void TriangleMesh::buildUniformGrid() {
    /**
     * TODO: Build uniform grid here
     * Note: you may need to build your own data structures in the accel_struct.hpp and accel_struct.cpp
     */

    // Init a big AABB bounding box containing all the meshes
    big_bounding_box.init(bounding_box.lb, bounding_box.ub);

    for (int i = 0; i < num_triangles; i++) {

        // Idx of the current triangle
        Eigen::Vector3i tri_idx(vertices_index[3*i], vertices_index[3*i+1], vertices_index[3*i+2]);

        // An AABB box containing the current triangle
        AABB tri_bounding_box(vertices[tri_idx(0)], vertices[tri_idx(1)], vertices[tri_idx(2)]);

        // // Get lb / ub idx of the box, if tri_box_lb = 43.7, big_box_lb = 34.6, idx = (int) 9.1 = 9
        // Eigen::Vector3i tri_lb_idx;
        // tri_lb_idx(0) = (int) floor(tri_bounding_box.lb(0) - big_bounding_box.lb(0));
        // tri_lb_idx(1) = (int) floor(tri_bounding_box.lb(1) - big_bounding_box.lb(1));
        // tri_lb_idx(2) = (int) floor(tri_bounding_box.lb(2) - big_bounding_box.lb(2));

        // Eigen::Vector3i tri_ub_idx;
        // tri_ub_idx(0) = (int) floor(tri_bounding_box.ub(0) - big_bounding_box.lb(0));
        // tri_ub_idx(1) = (int) floor(tri_bounding_box.ub(1) - big_bounding_box.lb(1));
        // tri_ub_idx(2) = (int) floor(tri_bounding_box.ub(2) - big_bounding_box.lb(2));

        // // Tranverse all the bounding boxes near the triangle, put them in if intersected with
        // // the triangle's bounding box
        // for (int x = tri_lb_idx(0); x <= tri_ub_idx(0); x++) {
        //     for (int y = tri_lb_idx(1); y <= tri_ub_idx(1); y++) {
        //         for (int z = tri_lb_idx(2); z <= tri_ub_idx(2); z++) {
        //             if (big_bounding_box.aabbList[x][y][z].isOverlap(tri_bounding_box) == true) {
        //                 // std::cout << "Pushed! ";
        //                 big_bounding_box.cellsList[x][y][z].push_back(Eigen::Vector3i(3*i, 3*i+1, 3*i+2));
        //             }
        //         }
        //     }
        // }

        for (int x = 0; x < big_bounding_box.dim(0); x++) {
            for (int y = 0; y < big_bounding_box.dim(1); y++) {
                for (int z = 0; z < big_bounding_box.dim(2); z++) {
                    if (big_bounding_box.aabbList[x][y][z].isOverlap(tri_bounding_box) == true) {
                        // std::cout << "Pushed! ";
                        big_bounding_box.cellsList[x][y][z].push_back(Eigen::Vector3i(3*i, 3*i+1, 3*i+2));
                    }
                }
            }
        }
    }
    
    // UNREACHABLE;
    
    has_grid = true;
}