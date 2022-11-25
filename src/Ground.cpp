

#include "Ground.h"

Ground::Ground() {
    init(Vector3d(0, 1, 0), -1);
}

void Ground::init(const Vector3d &m_normal, const double &d) {
    normal = m_normal;
    D = d;
}


double Ground::calculateGapFromObj(Mesh &mesh, const int &vId) const {
    double dist = normal.dot(mesh.pos[vId]) - D;
    return dist * dist;
}

//double Ground::calculateConstraintVal(const mesh3D& mesh, const int& vId) {
//    double dist = normal.dot(mesh.vertexes[vId]) - D;
//    return dist * dist;
//}

void Ground::calculateActivateSet(Mesh &mesh,
                                  std::vector<int> &Environment_ActiveSet, double d_hat) const {
    Environment_ActiveSet.resize(0);
#ifdef TEST_USE_TBB
    tbb::spin_mutex Mutex;
    tbb::parallel_for(0, (int) surfVerts.size(), 1, [&](int i)
#else
    for (int i = 0; i < mesh.surfVerts.size(); i++)
#endif

    {
        double dis = calculateGapFromObj(mesh, mesh.surfVerts[i]);
        if (dis < d_hat) {
#ifdef TEST_USE_TBB
            Mutex.lock();
#endif
            Environment_ActiveSet.push_back(mesh.surfVerts[i]);
#ifdef TEST_USE_TBB
            Mutex.unlock();
#endif
        }
    }
#ifdef TEST_USE_TBB
    );
#endif
}
