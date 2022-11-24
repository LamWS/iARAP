//
// Created by 林焕承 on 24/11/2022.
//

#ifndef ARAP_GROUND_H
#define ARAP_GROUND_H

#include "Eigen/Dense"
#include "mesh.h"

using namespace Eigen;

class Ground {
public: // data
    Vector3d normal;
    double D;
public:
    Ground();

    void init(const Vector3d &m_normal, const double &d);

    double calculateGapFromObj(Mesh &mesh, const int &vId) const;

    void calculateActivateSet(Mesh &mesh,
                              std::vector<int> &Environment_ActiveSet, double d_hat) const;
    //void calculateConstraintVal();
};


#endif //ARAP_GROUND_H
