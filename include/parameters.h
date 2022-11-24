//
// Created by 林焕承 on 22/11/2022.
//

#ifndef ARAP_PARAMETERS_H
#define ARAP_PARAMETERS_H
#define YoungModulus 1500
#define PoissonRate 0.3
namespace FEM {
    const static double mu = YoungModulus / (2 * (1 + PoissonRate));
    const static double lambda = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
}
namespace IPC {
    const static double Kappa = 1e15;
    const static double d_hat = 1e-3;
};
#endif //ARAP_PARAMETERS_H
