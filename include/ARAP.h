//
// Created by 林焕承 on 22/11/2022.
//

#ifndef ARAP_ARAP_H
#define ARAP_ARAP_H

#include <Eigen/Dense>
#include "quartic.h"

using namespace Eigen;


void ComputeEigenvector0(double a00, double a01, double a02, double a11, double a12, double a22,
                         double eval0, Vector3d &evec0) {
    Vector3d row0 = {a00 - eval0, a01, a02};
    Vector3d row1 = {a01, a11 - eval0, a12};
    Vector3d row2 = {a02, a12, a22 - eval0};
    Vector3d r0xr1 = row0.cross(row1);
    Vector3d r0xr2 = row0.cross(row2);
    Vector3d r1xr2 = row1.cross(row2);
    double d0 = r0xr1.dot(r0xr1);
    double d1 = r0xr2.dot(r0xr2);
    double d2 = r1xr2.dot(r1xr2);

    double dmax = d0;
    int32_t imax = 0;
    if (d1 > dmax) {
        dmax = d1;
        imax = 1;
    }
    if (d2 > dmax) {
        imax = 2;
    }

    if (imax == 0) {
        evec0 = r0xr1 / std::sqrt(d0);
    } else if (imax == 1) {
        evec0 = r0xr2 / std::sqrt(d1);
    } else {
        evec0 = r1xr2 / std::sqrt(d2);
    }
}

void ComputeOrthogonalComplement(Vector3d const &W,
                                 Vector3d &U, Vector3d &V) {
    double invLength;
    if (std::fabs(W[0]) > std::fabs(W[1])) {
        invLength = (double) 1 / std::sqrt(W[0] * W[0] + W[2] * W[2]);
        U = {-W[2] * invLength, (double) 0, +W[0] * invLength};
    } else {
        invLength = (double) 1 / std::sqrt(W[1] * W[1] + W[2] * W[2]);
        U = {(double) 0, +W[2] * invLength, -W[1] * invLength};
    }
    V = W.cross(U);
}

void ComputeEigenvector1(double a00, double a01, double a02, double a11, double a12, double a22,
                         Vector3d const &evec0, double eval1, Vector3d &evec1) {
    Vector3d U, V;
    ComputeOrthogonalComplement(evec0, U, V);

    Vector3d AU =
            {
                    a00 * U[0] + a01 * U[1] + a02 * U[2],
                    a01 * U[0] + a11 * U[1] + a12 * U[2],
                    a02 * U[0] + a12 * U[1] + a22 * U[2]
            };

    Vector3d AV =
            {
                    a00 * V[0] + a01 * V[1] + a02 * V[2],
                    a01 * V[0] + a11 * V[1] + a12 * V[2],
                    a02 * V[0] + a12 * V[1] + a22 * V[2]
            };

    double m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
    double m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
    double m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eval1;

    double absM00 = std::fabs(m00);
    double absM01 = std::fabs(m01);
    double absM11 = std::fabs(m11);
    double maxAbsComp;
    if (absM00 >= absM11) {
        maxAbsComp = std::max(absM00, absM01);
        if (maxAbsComp > (double) 0) {
            if (absM00 >= absM01) {
                m01 /= m00;
                m00 = (double) 1 / std::sqrt((double) 1 + m01 * m01);
                m01 *= m00;
            } else {
                m00 /= m01;
                m01 = (double) 1 / std::sqrt((double) 1 + m00 * m00);
                m00 *= m01;
            }
            evec1 = m01 * U - m00 * V;
        } else {
            evec1 = U;
        }
    } else {
        maxAbsComp = std::max(absM11, absM01);
        if (maxAbsComp > (double) 0) {
            if (absM11 >= absM01) {
                m01 /= m11;
                m11 = (double) 1 / std::sqrt((double) 1 + m01 * m01);
                m01 *= m11;
            } else {
                m11 /= m01;
                m01 = (double) 1 / std::sqrt((double) 1 + m11 * m11);
                m11 *= m01;
            }
            evec1 = m11 * U - m01 * V;
        } else {
            evec1 = U;
        }
    }
}

double compute_ARAP_energy(const Matrix3d &F) {
    double i1 = F.squaredNorm();
    double i2 = (F.transpose() * F).squaredNorm();
    double i3 = (F.transpose() * F).determinant();
    double J = F.determinant();
    double a = 0;
    double b = -2 * i1;
    double c = -8 * J;
    double d = i1 * i1 - 2 * (i1 * i1 - i2);
//        std::cout << a << " " << b << " " << c << " " << d << std::endl;
    std::complex<double> *solutions = solve_quartic(a, b, c, d);
    // in this quartic solver, the solutions can be reduced to the singular values as follows:
    // solutions[0] = s1 - s2 - s3
    // solutions[1] = s2 - s1 - s3
    // solutions[2] = s1 + s2 + s3
    // solutions[3] = s3 - s1 - s2

    double f = solutions[2].real();
    delete[] solutions;
    return 0.5 * (i1 - 2 * f + 3);
}

Eigen::Matrix<double, 9, 1> compute_ARAP_gradient(const Matrix3d &F) {
    Eigen::Matrix<double, 9, 1> g1;
    g1.block<3, 1>(0, 0) = 2 * F.col(0);
    g1.block<3, 1>(3, 0) = 2 * F.col(1);
    g1.block<3, 1>(6, 0) = 2 * F.col(2);
    Eigen::Matrix3d mat_g2 = 4 * F * F.transpose() * F;
    Eigen::Matrix<double, 9, 1> g2;
    g2.block<3, 1>(0, 0) = mat_g2.col(0);
    g2.block<3, 1>(3, 0) = mat_g2.col(1);
    g2.block<3, 1>(6, 0) = mat_g2.col(2);
    double J = F.determinant();
    Eigen::Matrix<double, 9, 1> gJ;
    gJ.block<3, 1>(0, 0) = F.col(1).cross(F.col(2));
    gJ.block<3, 1>(3, 0) = F.col(2).cross(F.col(0));
    gJ.block<3, 1>(6, 0) = F.col(0).cross(F.col(1));
    Eigen::Matrix<double, 9, 1> g3 = 2 * J * gJ;
    double i1 = F.squaredNorm();
    double i2 = (F.transpose() * F).squaredNorm();
    double i3 = (F.transpose() * F).determinant();
    double a = 0;
    double b = -2 * i1;
    double c = -8 * J;
    double d = i1 * i1 - 2 * (i1 * i1 - i2);

    std::complex<double> *solutions = solve_quartic(a, b, c, d);
    double f = solutions[2].real();
    double f1 = (2 * f * f + 2 * i1) / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
    double f2 = -2 / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
    double fJ = (8 * f) / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
    delete[] solutions;
    return 0.5 * (g1 - 2 * (f1 * g1 + f2 * g2 + fJ * gJ));
}

Matrix<double, 9, 9> compute_ARAP_Hessian(const Matrix3d &F) {
    double i1 = F.squaredNorm();
    double i2 = (F.transpose() * F).squaredNorm();
    double i3 = (F.transpose() * F).determinant();
    double J = F.determinant();
    double a = 0;
    double b = -2 * i1;
    double c = -8 * J;
    double d = i1 * i1 - 2 * (i1 * i1 - i2);

    std::complex<double> *solutions = solve_quartic(a, b, c, d);
    double x1 = solutions[0].real();
    double x2 = solutions[1].real();
    double x3 = solutions[3].real();
    double x4 = solutions[2].real();
    double sig1 = (x1 + x4) / 2;
    double sig2 = (x2 + x4) / 2;
    double sig3 = (x3 + x4) / 2;
    double f = solutions[2].real();
    // in this quartic solver, the solutions can be reduced to the singular values as follows:
    // solutions[0] = s1 - s2 - s3
    // solutions[1] = s2 - s1 - s3
    // solutions[2] = s1 + s2 + s3
    // solutions[3] = s3 - s1 - s2

    delete[] solutions;
    if (sig1 > sig2) {
        std::swap(sig1, sig2);
    }
    if (sig1 > sig3) {
        std::swap(sig1, sig3);
    }
    if (sig2 > sig3) {
        std::swap(sig2, sig3);
    }
    double f1 = (2 * f * f + 2 * i1) / (4 * f * f * f - 4 * i1 * f - 8 * J);
    double f2 = -2 / (4 * f * f * f - 4 * i1 * f - 8 * J);
    double fJ = (8 * f) / (4 * f * f * f - 4 * i1 * f - 8 * J);
    Eigen::Matrix3d g1 = 2 * F;
    Eigen::Matrix3d g2 = 4 * F * F.transpose() * F;
    Eigen::Matrix3d gJ;
    gJ.col(0) = F.col(1).cross(F.col(2));
    gJ.col(1) = F.col(2).cross(F.col(0));
    gJ.col(2) = F.col(0).cross(F.col(1));
    Matrix3d R = f1 * g1 + f2 * g2 + fJ * gJ;
    Matrix3d S = R.transpose() * F;
//                gte::NISymmetricEigensolver3x3<float> ses;
//                std::array<float, 3> eval{};
//                std::array<std::array<T, 3>, 3> evec{};
//                ses(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), 0, eval, evec);
//                cout << "eval:" << endl;
//                cout << eval[0] << " " << eval[1] << " " << eval[2] << endl;

//                cout << S << endl;
    Vector3d V0, V1, V2;

    Matrix3d V;
    double norm = S(0, 1) * S(0, 1) + S(0, 2) * S(0, 2) + S(1, 2) * S(1, 2);
    if (norm > 0) {
        double q = (S(0, 0) + S(1, 1) + S(2, 2)) / (double) 3;
        double b00 = S(0, 0) - q;
        double b11 = S(1, 1) - q;
        double b22 = S(2, 2) - q;
        double p = std::sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * (double) 2) / (double) 6);
        double c00 = b11 * b22 - S(1, 2) * S(1, 2);
        double c01 = S(0, 1) * b22 - S(1, 2) * S(0, 2);
        double c02 = S(0, 1) * S(1, 2) - b11 * S(0, 2);
        double det = (b00 * c00 - S(0, 1) * c01 + S(0, 2) * c02) / (p * p * p);

        double halfDet = det * (double) 0.5;
        halfDet = std::min(std::max(halfDet, (double) -1), (double) 1);

        if (halfDet >= (double) 0) {
            ComputeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sig3, V2);
            ComputeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V2, sig2, V1);
            V0 = V1.cross(V2);
        } else {
            ComputeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sig1, V0);
            ComputeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V0, sig2, V1);
            V2 = V0.cross(V1);
        }
        V.col(2) = V0;
        V.col(1) = V1;
        V.col(0) = V2;
    } else {
        // The matrix is diagonal.
        V = Matrix3d::Identity();
    }
    Matrix3d Sigma = Matrix3d::Identity();
    Sigma(0, 0) = sig3;
    Sigma(1, 1) = sig2;
    Sigma(2, 2) = sig1;

    Matrix3d U = R * V;

    Eigen::Matrix3d T0{
            {0, -1, 0},
            {1, 0,  0},
            {0, 0,  0}
    };
    T0 = (1 / sqrt(2)) * U * T0 * V.transpose();
    Eigen::Matrix3d T1{
            {0, 0,  0},
            {0, 0,  1},
            {0, -1, 0}
    };
    T1 = (1 / sqrt(2)) * U * T1 * V.transpose();
    Eigen::Matrix3d T2{
            {0,  0, 1},
            {0,  0, 0},
            {-1, 0, 0}
    };
    T2 = (1 / sqrt(2)) * U * T2 * V.transpose();
    Eigen::Matrix<double, 9, 1> t0, t1, t2;
    t0.block<3, 1>(0, 0) = T0.col(0);
    t0.block<3, 1>(3, 0) = T0.col(1);
    t0.block<3, 1>(6, 0) = T0.col(2);
    t1.block<3, 1>(0, 0) = T1.col(0);
    t1.block<3, 1>(3, 0) = T1.col(1);
    t1.block<3, 1>(6, 0) = T1.col(2);
    t2.block<3, 1>(0, 0) = T2.col(0);
    t2.block<3, 1>(3, 0) = T2.col(1);
    t2.block<3, 1>(6, 0) = T2.col(2);
    double sx = Sigma(0, 0);
    double sy = Sigma(1, 1);
    double sz = Sigma(2, 2);
    double lambda0 = 2 / (sx + sy);
    double lambda1 = 2 / (sz + sy);
    double lambda2 = 2 / (sx + sz);

    if (sx + sy < 2)lambda0 = 1;
    if (sz + sy < 2)lambda1 = 1;
    if (sx + sz < 2)lambda2 = 1;
//
    Eigen::Matrix<double, 9, 9> SH = Eigen::Matrix<double, 9, 9>::Identity();
//////
    SH -= lambda0 * (t0 * t0.transpose());
    SH -= lambda1 * (t1 * t1.transpose());
    SH -= lambda2 * (t2 * t2.transpose());
    return SH;
}

#endif //ARAP_ARAP_H
