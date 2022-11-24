
#ifndef ARAP_IPC_H
#define ARAP_IPC_H

#include <cmath>

namespace IPC {
    void compute_b(double d, double dHat, double &b) {
        b = -(d - dHat) * (d - dHat) * log(d / dHat);
    }

    void compute_g_b(double d, double dHat, double &g) {
        double t = d - dHat;
        g = t * std::log(d / dHat) * -2.0 - (t * t) / d;
    }

    void compute_H_b(double d, double dHat, double &H) {
        double t = d - dHat;
        H = (std::log(d / dHat) * -2.0 - t * 4.0 / d) + 1.0 / (d * d) * (t * t);
    }
}

#endif //ARAP_IPC_H
