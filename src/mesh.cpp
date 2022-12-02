#include <fstream>
#include <cfloat>
#include <map>

#include "mesh.h"
#include "parameters.h"

using namespace Eigen;

void split(const std::string &str, std::vector<std::string> &v, const std::string &spacer) {
    int pos1, pos2;
    int len = spacer.length();
    pos1 = 0;
    pos2 = str.find(spacer);
    while (pos2 != std::string::npos) {
        v.push_back(str.substr(pos1, pos2 - pos1));
        pos1 = pos2 + len;
        pos2 = str.find(spacer, pos1);
    }
    if (pos1 != str.length())
        v.push_back(str.substr(pos1));
}

double calculateVolum(const Vector3d &x0, const Vector3d &x1, const Vector3d &x2, const Vector3d &x3) {
    double o1x = x1[0] - x0[0];
    double o1y = x1[1] - x0[1];
    double o1z = x1[2] - x0[2];
    Vector3d OA = Vector3d(o1x, o1y, o1z);

    double o2x = x2[0] - x0[0];
    double o2y = x2[1] - x0[1];
    double o2z = x2[2] - x0[2];
    Vector3d OB = Vector3d(o2x, o2y, o2z);

    double o3x = x3[0] - x0[0];
    double o3y = x3[1] - x0[1];
    double o3z = x3[2] - x0[2];
    Vector3d OC = Vector3d(o3x, o3y, o3z);

    Vector3d heightDir = OA.cross(OB);
    double bottomArea = heightDir.norm();
    heightDir.normalize();

    double volum = bottomArea * abs(heightDir.dot(OC)) / 6;
    return volum;
}

bool Mesh::read_tetramesh(const std::string &filename) {
    surfEdges.clear();
    dynamic_surfEdges.clear();
    surfVerts.clear();
    std::ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
    }

    double x, y, z;
    int index0, index1, index2, index3;
    std::string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;
        if (line.find("$Nodes") != std::string::npos) {
            getline(ifs, line);
            nodeNumber = atoi(line.c_str());

            auto xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
            auto xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;
            for (int i = 0; i < nodeNumber; i++) {
                getline(ifs, line);
                std::vector<std::string> nodePos;
                std::string spacer = " ";
                split(line, nodePos, spacer);
                x = atof(nodePos[1].c_str());
                y = atof(nodePos[2].c_str());
                z = atof(nodePos[3].c_str());
                Vector3d d_velocity = Vector3d(0, 0, 0);
                Vector3d vertex = Vector3d(x, y, z);
                Vector3d compress_vertex = Vector3d(x, y, z);
//                if (y > 0.156) {
//                    compress_vertex.y() += (compress_vertex.y() - 0.156) * 2;
//                }
//                if (y > 0.5)compress_vertex.x() += 0.1;
//                double x1 = pow(2, 1.0 / 3.0) * ((double) rand() / (RAND_MAX));
//                double y1 = pow(2, 1.0 / 3.0) * ((double) rand() / (RAND_MAX));
//                double z1 = pow(2, 1.0 / 3.0) * ((double) rand() / (RAND_MAX));
//                compress_vertex.x() = x1;
//                compress_vertex.y() = y1;
//                compress_vertex.z() = z1;

                Matrix3d Constraint;
                Constraint.setIdentity();
                Vector3d force = Vector3d(0, 0, 0);
                Vector3d velocity = Vector3d(0, 0, 0);
                Vector3d d_pos = Vector3d(0, 0, 0);
                double mass = 0;
                pos.push_back(compress_vertex);
                res_pos.push_back(vertex);
//                forces.push_back(force);
//                velocities.push_back(velocity);
//                Constraints.push_back(Constraint);
//                isNBC.push_back(false);
//                d_velocities.push_back(d_velocity);
//                masses.push_back(mass);
//                isDelete.push_back(false);
//                d_positions.push_back(d_pos);
//                externalForce.push_back(Vector3d(0, 0, 0));

                Vector3d pos = vertex;
                if (xmin > pos[0]) xmin = pos[0];
                if (ymin > pos[1]) ymin = pos[1];
                if (zmin > pos[2]) zmin = pos[2];
                if (xmax < pos[0]) xmax = pos[0];
                if (ymax < pos[1]) ymax = pos[1];
                if (zmax < pos[2]) zmax = pos[2];
            }
//            minConer = Vector3d(xmin, ymin, zmin);
//            maxConer = Vector3d(xmax, ymax, zmax);
        }

        if (line.find("$Elements") != std::string::npos) {
            getline(ifs, line);
            elementNumber = atoi(line.c_str());
            for (int i = 0; i < elementNumber; i++) {
                getline(ifs, line);
                std::vector<std::string> elementIndexex;
                std::string spacer = " ";
                split(line, elementIndexex, spacer);
                index0 = atoi(elementIndexex[3].c_str()) - 1;
                index1 = atoi(elementIndexex[4].c_str()) - 1;
                index2 = atoi(elementIndexex[5].c_str()) - 1;
                index3 = atoi(elementIndexex[6].c_str()) - 1;

                Eigen::Vector4i tetrahedra;
                tetrahedra(0) = index0;
                tetrahedra(1) = index1;
                tetrahedra(2) = index2;
                tetrahedra(3) = index3;
                tets.push_back(tetrahedra);
            }
            break;
        }
    }
    ifs.close();
//
    std::map<std::set<int>, int> tri_map;
    std::map<std::set<int>, std::tuple<int, int, int>> tri_sort_map;
    for (auto tet: tets) {
        std::vector<int> tet_vec;
        tet_vec.emplace_back(tet.x());
        tet_vec.emplace_back(tet.y());
        tet_vec.emplace_back(tet.z());
        tet_vec.emplace_back(tet.w());
        std::sort(tet_vec.begin(), tet_vec.end());
        std::set<int> s1, s2, s3, s4;

        s1.insert(tet.x());
        s1.insert(tet.y());
        s1.insert(tet.z());

        s2.insert(tet.x());
        s2.insert(tet.y());
        s2.insert(tet.w());

        s3.insert(tet.x());
        s3.insert(tet.z());
        s3.insert(tet.w());

        s4.insert(tet.y());
        s4.insert(tet.z());
        s4.insert(tet.w());

        tri_map[s1] += 1;
        tri_sort_map[s1] = std::tuple<int, int, int>(tet.x(), tet.y(), tet.z());
        tri_map[s2] += 1;
        tri_sort_map[s2] = std::tuple<int, int, int>(tet.x(), tet.w(), tet.y());
        tri_map[s3] += 1;
        tri_sort_map[s3] = std::tuple<int, int, int>(tet.x(), tet.z(), tet.w());
        tri_map[s4] += 1;
        tri_sort_map[s4] = std::tuple<int, int, int>(tet.y(), tet.w(), tet.z());
    }

    for (auto s: tri_map) {
        if (s.second == 1) {
            auto sur = tri_sort_map[s.first];
            surface.emplace_back(std::get<0>(sur), std::get<1>(sur), std::get<2>(sur), 0);
        }
    }

    velocities.resize(pos.size(), Eigen::Vector3d::Zero());

    std::set<int> pos_set;
    for (auto s: surface) {
        pos_set.insert(s.x());
        pos_set.insert(s.y());
        pos_set.insert(s.z());
    }
    for (auto i: pos_set)surfVerts.push_back(i);
    //V_prev = vertexes;
    std::set<std::pair<int, int>> edge_set;
    for (auto tri: surface) {
        auto x = tri.x();
        auto y = tri.y();
        auto z = tri.z();
        if (x < y)edge_set.insert(std::make_pair(x, y));
        else edge_set.insert(std::make_pair(y, x));

        if (y < z)edge_set.insert(std::make_pair(y, z));
        else edge_set.insert(std::make_pair(z, y));

        if (x < z)edge_set.insert(std::make_pair(x, z));
        else edge_set.insert(std::make_pair(z, x));
    }
    for (auto p: edge_set) {
        surfEdges.emplace_back(p.first, p.second);
        dynamic_surfEdges.emplace_back(p.first, p.second);
    }
    std::vector<int> adj_v;
    std::vector<double> adj_volume;
    adj_v.resize(pos.size(), 0);
    adj_volume.resize(pos.size(), 0);
    double sum = 0;
    mass.resize(pos.size(), 0);
    for (auto tet: tets) {
        double v = calculateVolum(res_pos[tet[0]], res_pos[tet[1]], res_pos[tet[2]], res_pos[tet[3]]);
        volume.emplace_back(v);
        mass[tet[0]] += v * FEM::density / 4;
        mass[tet[1]] += v * FEM::density / 4;
        mass[tet[2]] += v * FEM::density / 4;
        mass[tet[3]] += v * FEM::density / 4;
        sum += v * FEM::density;
    }
    meanMass = sum / pos.size();
    return true;
}

Eigen::Vector3d Mesh::get_maxCorner() {
    auto maxx = DBL_MIN, maxy = DBL_MIN, maxz = DBL_MIN;
    for (auto ver: pos) {
        maxx = std::max(maxx, ver.x());
        maxy = std::max(maxy, ver.y());
        maxz = std::max(maxz, ver.z());
    }
    return Eigen::Vector3d(maxx, maxy, maxz);
}

Eigen::Vector3d Mesh::get_minCorner() {
    auto minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
    for (auto ver: pos) {
        minx = std::min(minx, ver.x());
        miny = std::min(miny, ver.y());
        minz = std::min(minz, ver.z());
    }
    return Eigen::Vector3d(minx, miny, minz);
}
