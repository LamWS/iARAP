//
// Created by 林焕承 on 22/11/2022.
//

#ifndef ARAP_COMMON_H
#define ARAP_COMMON_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

class BHessian {
public:
    std::vector<int> D1Index;//pIndex, DpeIndex, DptIndex;
    std::vector<Vector3i> D3Index;
    std::vector<Vector4i> D4Index;
    std::vector<Vector2i> D2Index;//, DempIndex;
    std::vector<Eigen::Matrix<double, 12, 12>> H12x12;//, HDee, HDpt;
    std::vector<Eigen::Matrix<double, 3, 3>> H3x3;//, HDpm, HDpm, HDpm;
    std::vector<Eigen::Matrix<double, 6, 6>> H6x6;//, HDeme;
    std::vector<Eigen::Matrix<double, 9, 9>> H9x9;

    std::vector<Eigen::Triplet<double>> toTriplets(int n) {
        Eigen::Matrix3d identity3 = Eigen::Matrix3d::Identity();
        std::vector<Triplet < double>>
        coefficients;
        coefficients.resize(9 * (D1Index.size() + D2Index.size() * 4 + D3Index.size() * 9 + D4Index.size() * 16),
                            Triplet<double>(0, 0, 0));
        int offset = 0;

#ifdef TEST_USE_TBB
        tbb::parallel_for(0, (int) D1Index.size(), 1, [&](int i)
#else
        for (int i = 0; i < D1Index.size(); i++)
#endif
        {
            if (D1Index[i] >= 0) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        coefficients[offset + i * 9 + j * 3 + k] = Triplet<double>(3 * D1Index[i] + j,
                                                                           3 * D1Index[i] + k,
                                                                           H3x3[i](j, k));
                    }
                }
            }
        }
#ifdef TEST_USE_TBB
        );
#endif
        offset = D1Index.size() * 9;

#ifdef TEST_USE_TBB
        tbb::parallel_for(0, (int) D2Index.size(), 1, [&](int i)
#else
        for (int i = 0; i < D2Index.size(); i++)
#endif
        {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (D2Index[i][0] >= 0)
                        coefficients[offset + i * 9 * 4 + j * 3 + k] = Triplet<double>(3 * D2Index[i][0] + j,
                                                                               3 * D2Index[i][0] + k,
                                                                               H6x6[i](j, k));
                    if (D2Index[i][0] >= 0 && D2Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 4 + 9 + j * 3 + k] = Triplet<double>(3 * D2Index[i][0] + j,
                                                                                   3 * D2Index[i][1] + k,
                                                                                   H6x6[i](j, k + 3));
                    if (D2Index[i][0] >= 0 && D2Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 4 + 18 + j * 3 + k] = Triplet<double>(3 * D2Index[i][1] + j,
                                                                                    3 * D2Index[i][0] + k,
                                                                                    H6x6[i](j + 3, k));
                    if (D2Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 4 + 27 + j * 3 + k] = Triplet<double>(3 * D2Index[i][1] + j,
                                                                                    3 * D2Index[i][1] + k,
                                                                                    H6x6[i](j + 3, k + 3));
                }
            }
        }
#ifdef TEST_USE_TBB
        );
#endif
        offset += D2Index.size() * 36;

#ifdef TEST_USE_TBB
        tbb::parallel_for(0, (int) D3Index.size(), 1, [&](int i)
#else
        for (int i = 0; i < D3Index.size(); i++)
#endif
        {
//                              IglUtils::makePD(H9x9[i]);
//            std::cout << D3Index[i].transpose() << std::endl;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (D3Index[i][0] >= 0)
                        coefficients[offset + i * 9 * 9 + j * 3 + k] = Triplet<double>(3 * D3Index[i][0] + j,
                                                                               3 * D3Index[i][0] + k,
                                                                               H9x9[i](j, k));
                    if (D3Index[i][0] >= 0 && D3Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 9 + 9 + j * 3 + k] = Triplet<double>(3 * D3Index[i][0] + j,
                                                                                   3 * D3Index[i][1] + k,
                                                                                   H9x9[i](j, k + 3));
                    if (D3Index[i][0] >= 0 && D3Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 9 + 18 + j * 3 + k] = Triplet<double>(3 * D3Index[i][0] + j,
                                                                                    3 * D3Index[i][2] + k,
                                                                                    H9x9[i](j, k + 6));
                    if (D3Index[i][0] >= 0 && D3Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 9 + 27 + j * 3 + k] = Triplet<double>(3 * D3Index[i][1] + j,
                                                                                    3 * D3Index[i][0] + k,
                                                                                    H9x9[i](j + 3, k));
                    if (D3Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 9 + 36 + j * 3 + k] = Triplet<double>(3 * D3Index[i][1] + j,
                                                                                    3 * D3Index[i][1] + k,
                                                                                    H9x9[i](j + 3, k + 3));

                    if (D3Index[i][1] >= 0 && D3Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 9 + 45 + j * 3 + k] = Triplet<double>(3 * D3Index[i][1] + j,
                                                                                    3 * D3Index[i][2] + k,
                                                                                    H9x9[i](j + 3, k + 6));
                    if (D3Index[i][2] >= 0 && D3Index[i][0] >= 0)
                        coefficients[offset + i * 9 * 9 + 54 + j * 3 + k] = Triplet<double>(3 * D3Index[i][2] + j,
                                                                                    3 * D3Index[i][0] + k,
                                                                                    H9x9[i](j + 6, k));
                    if (D3Index[i][2] >= 0 && D3Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 9 + 63 + j * 3 + k] = Triplet<double>(3 * D3Index[i][2] + j,
                                                                                    3 * D3Index[i][1] + k,
                                                                                    H9x9[i](j + 6, k + 3));
                    if (D3Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 9 + 72 + j * 3 + k] = Triplet<double>(3 * D3Index[i][2] + j,
                                                                                    3 * D3Index[i][2] + k,
                                                                                    H9x9[i](j + 6, k + 6));
                }
            }
        }
#ifdef TEST_USE_TBB
        );
#endif
        offset += D3Index.size() * 81;
#ifdef TEST_USE_TBB
        tbb::parallel_for(0, (int) D4Index.size(), 1, [&](int i)
#else
        for (int i = 0; i < D4Index.size(); i++)
#endif
        {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (D4Index[i][0] >= 0)
                        coefficients[offset + i * 9 * 16 + j * 3 + k] = Triplet<double>(3 * D4Index[i][0] + j,
                                                                                3 * D4Index[i][0] + k,
                                                                                H12x12[i](j, k));
                    if (D4Index[i][0] >= 0 && D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 9 + j * 3 + k] = Triplet<double>(3 * D4Index[i][0] + j,
                                                                                    3 * D4Index[i][1] + k,
                                                                                    H12x12[i](j, k + 3));
                    if (D4Index[i][0] >= 0 && D4Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 16 + 18 + j * 3 + k] = Triplet<double>(3 * D4Index[i][0] + j,
                                                                                     3 * D4Index[i][2] + k,
                                                                                     H12x12[i](j, k + 6));
                    if (D4Index[i][0] >= 0 && D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 27 + j * 3 + k] = Triplet<double>(3 * D4Index[i][0] + j,
                                                                                     3 * D4Index[i][3] + k,
                                                                                     H12x12[i](j, k + 9));
                    if (D4Index[i][0] >= 0 && D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 36 + j * 3 + k] = Triplet<double>(3 * D4Index[i][1] + j,
                                                                                     3 * D4Index[i][0] + k,
                                                                                     H12x12[i](j + 3, k));
                    if (D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 45 + j * 3 + k] = Triplet<double>(3 * D4Index[i][1] + j,
                                                                                     3 * D4Index[i][1] + k,
                                                                                     H12x12[i](j + 3, k + 3));
                    if (D4Index[i][2] >= 0 && D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 54 + j * 3 + k] = Triplet<double>(3 * D4Index[i][1] + j,
                                                                                     3 * D4Index[i][2] + k,
                                                                                     H12x12[i](j + 3, k + 6));
                    if (D4Index[i][3] >= 0 && D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 63 + j * 3 + k] = Triplet<double>(3 * D4Index[i][1] + j,
                                                                                     3 * D4Index[i][3] + k,
                                                                                     H12x12[i](j + 3, k + 9));
                    if (D4Index[i][0] >= 0 && D4Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 16 + 72 + j * 3 + k] = Triplet<double>(3 * D4Index[i][2] + j,
                                                                                     3 * D4Index[i][0] + k,
                                                                                     H12x12[i](j + 6, k));
                    if (D4Index[i][2] >= 0 && D4Index[i][1] >= 0)
                        coefficients[offset + i * 9 * 16 + 81 + j * 3 + k] = Triplet<double>(3 * D4Index[i][2] + j,
                                                                                     3 * D4Index[i][1] + k,
                                                                                     H12x12[i](j + 6, k + 3));
                    if (D4Index[i][2] >= 0)
                        coefficients[offset + i * 9 * 16 + 90 + j * 3 + k] = Triplet<double>(3 * D4Index[i][2] + j,
                                                                                     3 * D4Index[i][2] + k,
                                                                                     H12x12[i](j + 6, k + 6));
                    if (D4Index[i][2] >= 0 && D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 99 + j * 3 + k] = Triplet<double>(3 * D4Index[i][2] + j,
                                                                                     3 * D4Index[i][3] + k,
                                                                                     H12x12[i](j + 6, k + 9));
                    if (D4Index[i][0] >= 0 && D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 108 + j * 3 + k] = Triplet<double>(3 * D4Index[i][3] + j,
                                                                                      3 * D4Index[i][0] + k,
                                                                                      H12x12[i](j + 9, k));
                    if (D4Index[i][1] >= 0 && D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 117 + j * 3 + k] = Triplet<double>(3 * D4Index[i][3] + j,
                                                                                      3 * D4Index[i][1] + k,
                                                                                      H12x12[i](j + 9, k + 3));
                    if (D4Index[i][2] >= 0 && D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 126 + j * 3 + k] = Triplet<double>(3 * D4Index[i][3] + j,
                                                                                      3 * D4Index[i][2] + k,
                                                                                      H12x12[i](j + 9, k + 6));
                    if (D4Index[i][3] >= 0)
                        coefficients[offset + i * 9 * 16 + 135 + j * 3 + k] = Triplet<double>(3 * D4Index[i][3] + j,
                                                                                      3 * D4Index[i][3] + k,
                                                                                      H12x12[i](j + 9, k + 9));
                }
            }
        }
#ifdef TEST_USE_TBB
        );
#endif
        return coefficients;
    }
};

#endif //ARAP_COMMON_H
