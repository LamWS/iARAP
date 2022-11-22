#ifndef ARAP_OPTIMIZER_H
#define ARAP_OPTIMIZER_H


#include "mesh.h"
#include "common.h"


class Optimizer {
public:
    Optimizer(Mesh &mesh);

    void initialize_FEM(Mesh &mesh);

    void compute_F(Mesh &mesh);

    void solve(Mesh &mesh, double dt, int total_frames);

    void compute_gradient_Hessian(Mesh &mesh, double dt, std::vector<Vector3d> &gradient, BHessian &BH);

    Eigen::Vector3d get_external_force(Mesh &mesh, int i, double dt);

    void calculateMovingDirection(Mesh &mesh, BHessian &BH, std::vector<Eigen::Vector3d> gradient,
                                  std::vector<Eigen::Vector3d> &direction);

    void line_search(Mesh &mesh, double dt, std::vector<Eigen::Vector3d> &direction, double &step);

    void step_forward(Mesh &mesh, std::vector<Eigen::Vector3d> dataV0, std::vector<Eigen::Vector3d> &direction, double &step);

    double compute_energy(Mesh &mesh, double dt);

    void export_obj(Mesh &mesh, int frame);

    std::vector<Eigen::Matrix<double, 9, 12>> tet_dfdx;
    std::vector<Eigen::Matrix<double, 3, 3>> tet_DmInv;
    std::vector<Eigen::Matrix<double, 3, 3>> tet_F;
    std::vector<Eigen::Vector3d> pos_prev;
    double lastEnergyVal = 0;

};


#endif //ARAP_OPTIMIZER_H
