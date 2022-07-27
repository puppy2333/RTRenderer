#pragma once
#include "accel_struct.hpp"


/**
 * AABB class
 */

AABB::AABB()
    : lb(0, 0, 0)
    , ub(0, 0, 0)
{
}
    
AABB::AABB(float lb_x, float lb_y, float lb_z, float ub_x, float ub_y, float ub_z)
{
    lb = Eigen::Vector3f(lb_x, lb_y, lb_z);
    ub = Eigen::Vector3f(ub_x, ub_y, ub_z);
}

AABB::AABB(Eigen::Vector3f lb, Eigen::Vector3f ub)
    : lb(lb)
    , ub(ub)
{
}

/* Construct AABB for a sphere */
AABB::AABB(const Eigen::Vector3f& pos, float radius)
{
    Eigen::Vector3f r(radius, radius, radius);
    lb = pos - r;
    ub = pos + r;
}

/* Construct AABB for a triangle, finding the min & max of all three vertexes */
AABB::AABB(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, const Eigen::Vector3f& v3)
{
    lb = v1.cwiseMin(v2).cwiseMin(v3);
    ub = v1.cwiseMax(v2).cwiseMax(v3);
}

/* Construct AABB by merging two AABBs */
AABB::AABB(const AABB& a, const AABB& b)
{
    lb = Eigen::Vector3f(a.lb.cwiseMin(b.lb));
    ub = Eigen::Vector3f(a.ub.cwiseMax(b.ub));
}

Eigen::Vector3f AABB::getCenter() const
{
    return (lb + ub) / 2;
}

float AABB::getDist(int c) const
{
    return ub[c] - lb[c];
}

float AABB::getVolume() const
{
    return getDist(2) * getDist(1) * getDist(0);
}

bool AABB::isOverlap(const AABB& a) const
{
    return ((a.lb[0] >= this->lb[0] && a.lb[0] <= this->ub[0]) || (this->lb[0] >= a.lb[0] && this->lb[0] <= a.ub[0])) &&
        ((a.lb[1] >= this->lb[1] && a.lb[1] <= this->ub[1]) || (this->lb[1] >= a.lb[1] && this->lb[1] <= a.ub[1])) &&
        ((a.lb[2] >= this->lb[2] && a.lb[2] <= this->ub[2]) || (this->lb[2] >= a.lb[2] && this->lb[2] <= a.ub[2]));

}

float AABB::diagonalLength() const
{
    return (ub - lb).norm();
}

bool AABB::rayIntersection(const Ray& ray, float& t_in, float& t_out) const
{
    float dir_frac_x = (ray.direction[0] == 0.0) ? 1.0e32 : 1.0f / ray.direction[0];
    float dir_frac_y = (ray.direction[1] == 0.0) ? 1.0e32 : 1.0f / ray.direction[1];
    float dir_frac_z = (ray.direction[2] == 0.0) ? 1.0e32 : 1.0f / ray.direction[2];

    float tx1 = (lb[0] - ray.origin[0]) * dir_frac_x;
    float tx2 = (ub[0] - ray.origin[0]) * dir_frac_x;
    float ty1 = (lb[1] - ray.origin[1]) * dir_frac_y;
    float ty2 = (ub[1] - ray.origin[1]) * dir_frac_y;
    float tz1 = (lb[2] - ray.origin[2]) * dir_frac_z;
    float tz2 = (ub[2] - ray.origin[2]) * dir_frac_z;

    t_in = std::max(std::max(std::min(tx1, tx2), std::min(ty1, ty2)), std::min(tz1, tz2));
    t_out = std::min(std::min(std::max(tx1, tx2), std::max(ty1, ty2)), std::max(tz1, tz2));

    /* When t_out < 0 and the ray is intersecting with AABB, the whole AABB is behind us */
    if (t_out < 0)
    {
        return false;
    }

    return t_out >= t_in;
}

BigAABB::BigAABB()
{
}

void BigAABB::init(Eigen::Vector3f lb_in, Eigen::Vector3f ub_in)
{
    lb = lb_in;
    ub = ub_in;
    dim(0) = (int) ceil(ub(0) - lb(0));
    dim(1) = (int) ceil(ub(1) - lb(1));
    dim(2) = (int) ceil(ub(2) - lb(2));

    std::cout << "lb: " << lb(0) << " " << lb(1) << " " << lb(2) << " \n";
    std::cout << "ub: " << ub(0) << " " << ub(1) << " " << ub(2) << " \n";
    std::cout << "size: " << dim(0) << " " << dim(1) << " " << dim(2) << "\n";

    aabbList = new AABB ** [dim(0)];
    for (int x = 0; x < dim(0); x++) {
        aabbList[x] = new AABB * [dim(1)];
        for (int y = 0; y < dim(1); y++) {
            aabbList[x][y] = new AABB [dim(2)];
        }
    }
    for (int x = 0; x < dim(0); x++) {
        for (int y = 0; y < dim(1); y++) {
            for (int z = 0; z < dim(2); z++) {
                aabbList[x][y][z].lb = lb + Eigen::Vector3f(x, y, z);
                aabbList[x][y][z].ub = aabbList[x][y][z].lb + Eigen::Vector3f(1, 1, 1);
            }
        }
    }

    cellsList = new std::vector<Eigen::Vector3i> ** [dim(0)];
    for (int x = 0; x < dim(0); x++) {
        cellsList[x] = new std::vector<Eigen::Vector3i> * [dim(1)];
        for (int y = 0; y < dim(1); y++) {
            cellsList[x][y] = new std::vector<Eigen::Vector3i> [dim(2)];
        }
    }
}

// BigBox::BigBox()
// {
// }

// void BigBox::init(Eigen::Vector3f lb_in, Eigen::Vector3f ub_in)
// {
//     lb = lb_in;
//     ub = ub_in;
//     dim(0) = (int) ceil(ub(0) - lb(0));
//     dim(1) = (int) ceil(ub(1) - lb(1));
//     dim(2) = (int) ceil(ub(2) - lb(2));

//     cellsList = new std::vector<Eigen::Vector3f> ** [dim(0)];
//     for (int x = 0; x < dim(0); x++) {
//         cellsList[x] = new std::vector<Eigen::Vector3f> * [dim(1)];
//         for (int y = 0; y < dim(1); y++) {
//             cellsList[x][y] = new std::vector<Eigen::Vector3f> [dim(2)];
//         }
//     }
// }