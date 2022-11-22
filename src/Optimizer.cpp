#include "Optimizer.h"
#include "ARAP.h"
#include "parameters.h"
#include "fstream"
#include "cfloat"
#include <iostream>

Optimizer::Optimizer(Mesh &mesh) {
    initialize_FEM(mesh);
}

void Optimizer::initialize_FEM(Mesh &mesh) {
    tet_dfdx.resize(mesh.tets.size());
    tet_DmInv.resize(mesh.tets.size());
    tet_F.resize(mesh.tets.size());
    for (int i = 0; i < mesh.tets.size(); i++) {
        Eigen::Vector3d x0 = mesh.res_pos[mesh.tets[i][0]];
        Eigen::Vector3d x1 = mesh.res_pos[mesh.tets[i][1]];
        Eigen::Vector3d x2 = mesh.res_pos[mesh.tets[i][2]];
        Eigen::Vector3d x3 = mesh.res_pos[mesh.tets[i][3]];
        Eigen::Matrix<double, 3, 3> Dm, DmInv;
        Dm << x1 - x0, x2 - x0, x3 - x0;
        tet_DmInv[i] = DmInv = Dm.inverse();
        const double m = DmInv(0, 0);
        const double n = DmInv(0, 1);
        const double o = DmInv(0, 2);
        const double p = DmInv(1, 0);
        const double q = DmInv(1, 1);
        const double r = DmInv(1, 2);
        const double s = DmInv(2, 0);
        const double t = DmInv(2, 1);
        const double u = DmInv(2, 2);
        const double t1 = -m - p - s;
        const double t2 = -n - q - t;
        const double t3 = -o - r - u;
        Eigen::Matrix<double, 9, 12> PFPx = Eigen::Matrix<double, 9, 12>::Zero();

        PFPx(0, 0) = t1;
        PFPx(0, 3) = m;
        PFPx(0, 6) = p;
        PFPx(0, 9) = s;
        PFPx(1, 1) = t1;
        PFPx(1, 4) = m;
        PFPx(1, 7) = p;
        PFPx(1, 10) = s;
        PFPx(2, 2) = t1;
        PFPx(2, 5) = m;
        PFPx(2, 8) = p;
        PFPx(2, 11) = s;
        PFPx(3, 0) = t2;
        PFPx(3, 3) = n;
        PFPx(3, 6) = q;
        PFPx(3, 9) = t;
        PFPx(4, 1) = t2;
        PFPx(4, 4) = n;
        PFPx(4, 7) = q;
        PFPx(4, 10) = t;
        PFPx(5, 2) = t2;
        PFPx(5, 5) = n;
        PFPx(5, 8) = q;
        PFPx(5, 11) = t;
        PFPx(6, 0) = t3;
        PFPx(6, 3) = o;
        PFPx(6, 6) = r;
        PFPx(6, 9) = u;
        PFPx(7, 1) = t3;
        PFPx(7, 4) = o;
        PFPx(7, 7) = r;
        PFPx(7, 10) = u;
        PFPx(8, 2) = t3;
        PFPx(8, 5) = o;
        PFPx(8, 8) = r;
        PFPx(8, 11) = u;
        tet_dfdx[i] = PFPx;
    }
}

void Optimizer::compute_F(Mesh &mesh) {
    for (int i = 0; i < mesh.tets.size(); i++) {
        Eigen::Vector3d x0 = mesh.pos[mesh.tets[i][0]];
        Eigen::Vector3d x1 = mesh.pos[mesh.tets[i][1]];
        Eigen::Vector3d x2 = mesh.pos[mesh.tets[i][2]];
        Eigen::Vector3d x3 = mesh.pos[mesh.tets[i][3]];
        Eigen::Matrix<double, 3, 3> Ds;
        Ds << x1 - x0, x2 - x0, x3 - x0;
        tet_F[i] = Ds * tet_DmInv[i];
    }
}

void Optimizer::solve(Mesh &mesh, double dt, int total_frames) {
    pos_prev = mesh.pos;
    for (int frame = 0; frame < total_frames; frame++) {
        lastEnergyVal = compute_energy(mesh, dt);
        int k = 0;
        int iterCap = 10000;
        for (; k < iterCap; ++k) {
            std::vector<Vector3d> gradient(mesh.pos.size(), Vector3d(0, 0, 0));
            BHessian BH;
            compute_gradient_Hessian(mesh, dt, gradient, BH);
            std::vector<Eigen::Vector3d> direction(mesh.pos.size(), Vector3d(0, 0, 0));

            calculateMovingDirection(mesh, BH, gradient, direction);
            double dist_to_converge = 0;
#ifdef TEST_USE_TBB
            distToOpt_PN = parallel_reduce(
                tbb::blocked_range<int>(0, (int) direction.size()), 0.0,
                [&](const tbb::blocked_range<int> &rg, double temp_max) {
                    for (int i = rg.begin(); i != rg.end(); i++) {
                        for (int jj = 0; jj < 3; jj++) {
                            if (temp_max < abs(direction[i][jj])) {
                                temp_max = abs(direction[i][jj]);
                            }
                        }
                    }
                    return temp_max;
                },
                [&](double left, double right) {
                    return left > right ? left : right;
                }
        );
#else
            for (int ii = 0; ii < direction.size(); ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    if (dist_to_converge < abs(direction[ii][jj])) {
                        dist_to_converge = abs(direction[ii][jj]);
                    }
                }
            }
#endif
            std::cout << "dist_to_converge" << std::endl;
            std::cout << dist_to_converge << std::endl;
            double bboxDiagSize2 = (mesh.get_maxCorner() - mesh.get_minCorner()).squaredNorm();
            bool gradVanish = (dist_to_converge < sqrt(1e-6 * bboxDiagSize2 * dt * dt));
            if (k && gradVanish) {
                break;
            }
            double step_size = 1;
            line_search(mesh, dt, direction, step_size);
        }
        for (int i = 0; i < mesh.pos.size(); i++) {
            mesh.velocities[i] = (mesh.pos[i] - pos_prev[i]) / dt;
        }
        printf("frame end %d\n", frame);
        export_obj(mesh, frame);
    }
}

void Optimizer::compute_gradient_Hessian(Mesh &mesh, double dt, std::vector<Vector3d> &gradient, BHessian &BH) {
    compute_F(mesh);;
    Matrix3d massIdentity;
    massIdentity.setIdentity();
    gradient.resize(mesh.pos.size(), Eigen::Vector3d::Zero());

#ifndef NoSimulation
#ifdef TEST_USE_TBB
    tbb::parallel_for(0, vertices_cnt, 1, [&](int vI)
#else
    for (int vI = 0; vI < mesh.pos.size(); vI++)
#endif
    {
        Eigen::Vector3d xt = pos_prev[vI];
        Eigen::Vector3d v = mesh.velocities[vI];
        Eigen::Vector3d x = mesh.pos[vI];
        Eigen::Vector3d fe = get_external_force(mesh, vI, dt);
        auto tmp = (mesh.mass[vI] * massIdentity *
                    (x - (xt + dt * v + dt * dt * (1 / mesh.mass[vI]) * massIdentity * fe)));
        gradient[vI] += tmp;
        // damping
//        if (FEM::damping > 0.0) {
//            gradient[vI] += FEM::damping * mesh.mass[vI] * massIdentity * (x - xt) * dt;
//        }
    }
#ifdef TEST_USE_TBB
    );
#endif
#endif
    //internal part
    double fsum = 0;


    std::vector<Eigen::Matrix<double, 12, 1>> gradient_vec;
    gradient_vec.resize(mesh.tets.size());
    BH.H12x12.resize(mesh.tets.size(), Eigen::Matrix<double, 12, 12>::Zero());
#ifdef TEST_USE_TBB
    tbb::parallel_for(0, (int) object->get_tets().size(), 1, [&](int ii)
#else
    for (int ii = 0; ii < mesh.tets.size(); ii++)
#endif
    {

        gradient_vec[ii] =
                dt * dt * FEM::mu * mesh.volume[ii] * tet_dfdx[ii].transpose() * compute_ARAP_gradient(tet_F[ii]);
        BH.H12x12[ii] =
                dt * dt * mesh.volume[ii] * tet_dfdx[ii].transpose() * compute_ARAP_Hessian(tet_F[ii]) * tet_dfdx[ii];
    }
#ifdef TEST_USE_TBB
    );
#endif
    BH.D4Index = mesh.tets;
    for (auto tetIndex = 0; tetIndex < BH.D4Index.size(); tetIndex++) {
        auto tet = BH.D4Index[tetIndex];
        for (int i = 0; i < 3; i++) {
            gradient[tet.x()][i] += gradient_vec[tetIndex](i);
            gradient[tet.y()][i] += gradient_vec[tetIndex](i + 3);
            gradient[tet.z()][i] += gradient_vec[tetIndex](i + 6);
            gradient[tet.w()][i] += gradient_vec[tetIndex](i + 9);
        }
    }
}

Eigen::Vector3d Optimizer::get_external_force(Mesh &mesh, int i, double dt) {
    return Vector3d(0, 0 * mesh.mass[i], 0);
}

void Optimizer::calculateMovingDirection(Mesh &mesh, BHessian &BH, std::vector<Eigen::Vector3d> gradient,
                                         std::vector<Eigen::Vector3d> &direction) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    std::vector<Triplet<double>> triplets = BH.toTriplets(mesh.pos.size());
    int offset = triplets.size();
    triplets.resize(offset + 3 * mesh.pos.size());
#ifdef TEST_USE_TBB
    tbb::parallel_for(0, 3 * mesh.pos.size(), 1, [&](int i)
#else
    for (int i = 0; i < 3 * mesh.pos.size(); i++)
#endif
    {

        triplets[offset + i] = Triplet<double>(i, i, mesh.mass[i / 3]);
    }
#ifdef TEST_USE_TBB
    );
#endif
    Eigen::SparseMatrix<double> sparseFEM(3 * mesh.pos.size(), 3 * mesh.pos.size());
    sparseFEM.setFromTriplets(triplets.begin(), triplets.end());
//    std::cout << sparseFEM << std::endl;

    auto vec_gradient = Eigen::VectorXd(gradient.size() * 3);
    for (int i = 0; i < gradient.size(); i++) {
        vec_gradient.block<3, 1>(3 * i, 0) = gradient[i];
    }
    auto result = Eigen::VectorXd(gradient.size());
    solver.compute(sparseFEM);
    result = solver.solve(vec_gradient);
    for (int i = 0; i < gradient.size(); i++) {
        direction[i] = result.block<3, 1>(3 * i, 0);
    }

}

void Optimizer::line_search(Mesh &mesh, double dt, std::vector<Eigen::Vector3d> &direction, double &step) {
    bool stopped = false;

    double c1m = 0.0;
    double armijoParam = 0;
    if (armijoParam > 0.0) {
#ifdef TEST_USE_TBB
        c1m = parallel_reduce(
                tbb::blocked_range<int>(0, mesh.pos.size()), 0.0,
                [&](const tbb::blocked_range<int> &rg, double temp_deltaE) {
                    for (int i = rg.begin(); i != rg.end(); i++) {
                        temp_deltaE += direction[i].dot(direction[i]);
                    }
                    return temp_deltaE;
                },
                [&](double left, double right) {
                    return left + right;
                }
        );
#else
        for (int i = 0; i < mesh.pos.size(); i++) {
            c1m += direction[i].dot(direction[i]);
        }
#endif
    }

    std::vector<Eigen::Vector3d> resultV0 = mesh.pos;
    step_forward(mesh, resultV0, direction, step);

    bool rehash = true;

//    sh->build(object);
//    std::cout << "step spatial hash:" << std::endl;
//    while (isIntersected(resultV0)) {
//        std::cout << "intersected" << std::endl;
//        profile::ccd_error += 1;
//        step /= 2.0;
//        stepForward(resultV0, direction, step);
//        sh->build(object);
//    }
//    std::cout << "step before line search:" << std::endl;
//    std::cout << step << std::endl;
    double alpha_before = step;
    rehash = false;
    double testingE;


//    buildCollisionSets();

    testingE = compute_energy(mesh, dt);
    int numOfLineSearch = 0;
    double LFStepSize = step;
    while ((testingE > lastEnergyVal + step * c1m) &&
           (step > 1e-8)) {
        ++numOfLineSearch;

        if (step == 0.0) {
            stopped = true;
            break;
        }
        step_forward(mesh, resultV0, direction, step);
        step /= 2.0;
//        buildCollisionSets();
        testingE = compute_energy(mesh, dt);
//        std::cout << "step during line search:" << std::endl;
//        std::cout << step << std::endl;

    }
    lastEnergyVal = testingE;
}

double Optimizer::compute_energy(Mesh &mesh, double dt) {
    compute_F(mesh);
    double energyVal = 0;
    for (int i = 0; i < mesh.tets.size(); i++) {
        energyVal += dt * dt * compute_ARAP_energy(tet_F[i]);
    }
    double deltaE = 0;
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    for (int i = 0; i < mesh.pos.size(); i++) {
        Eigen::Vector3d xt = pos_prev[i];
        Eigen::Vector3d v = mesh.velocities[i];
        Eigen::Vector3d x = mesh.pos[i];
        Eigen::Vector3d fe = get_external_force(mesh, i, dt);
        double tep = (x - (xt + dt * v + dt * dt * (1 / mesh.mass[i]) * identity * fe)).squaredNorm() *
                     mesh.mass[i] / 2.0;
        deltaE += tep;
    }
    energyVal += deltaE;
    return energyVal;
}

void Optimizer::step_forward(Mesh &mesh, std::vector<Eigen::Vector3d> dataV0, std::vector<Eigen::Vector3d> &direction,
                             double &step) {
    for (int vI = 0; vI < mesh.pos.size(); vI++) {
        mesh.pos[vI] = dataV0[vI] - step * direction[vI];
    }
}

void Optimizer::export_obj(Mesh &mesh, int frame) {
    std::ofstream cloth_stream("../data/output/surface_" + std::to_string(frame) + ".obj");
    cloth_stream << "# Generated by lam"
                 << "\n";
    cloth_stream << std::fixed << std::setprecision(6) << "cloth\n";
    auto cloth_vertices = mesh.pos;
    auto cloth_triangles = mesh.surface;

    for (auto vec: cloth_vertices) {
        cloth_stream << "v " << vec(0) << " " << vec(1) << " "
                     << vec(2) << "\n";
    }
    cloth_stream << "s 1\n";
    for (auto tri: cloth_triangles) {
        cloth_stream << "f " << tri(0) + 1 << " " << tri(1) + 1
                     << " " << tri(2) + 1 << "\n";
    }
    cloth_stream.close();
}

