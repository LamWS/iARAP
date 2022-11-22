#ifndef ARAP_MESH_H
#define ARAP_MESH_H

#include "Eigen/Dense"
#include <vector>
#include <set>

class Mesh {
public:
    std::vector<Eigen::Vector3d> pos;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<Eigen::Vector3d> res_pos;
    std::vector<Eigen::Vector4i> tets;


    std::vector<std::pair<uint64_t, uint64_t>> dynamic_surfEdges;
    std::vector<std::pair<uint64_t, uint64_t>> surfEdges;
    std::vector<uint64_t> surfVerts;
    std::vector<Eigen::Vector4i> surface;
    std::vector<double> mass;
    std::vector<double> volume;
    double meanMass = 0;

    bool read_tetramesh(const std::string& filename);
    Eigen::Vector3d get_maxCorner();
    Eigen::Vector3d get_minCorner();
};


#endif //ARAP_MESH_H
