# Isotropic-ARAP-energy-using-Cauchy-Green-invariants

Demo code for volumetric dynamic simulation in paper Isotropic-ARAP-energy-using-Cauchy-Green-invariants. Given a volumetric mesh as the input, this code output the simulation of dropping the mesh to the ground.

![](https://github.com/LamWS/iARAP/raw/main/file/teaser.png)

[Paper], [supplementary], [video], and [code] are now available.

## Usage

```shell
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./arap ../data/<cube5.msh>
```

The result sequence will be output to `build/output` by default.

Please check the code in `main.cpp` and `parameters.h` for more configuration.

## Abstract

Isotropic As-Rigid-As-Possible (ARAP) energy has been popular for shape editing, mesh parametrisation and soft-body
simulation for almost two decades.
However, a formulation using Cauchy-Green (CG) invariants has always been unclear, due to a rotation-polluted trace term that cannot be directly expressed using these invariants.
We show how this incongruent trace term can be understood via an implicit relationship to the CG invariants.
Our analysis reveals this relationship to be a polynomial where the roots equate to the trace term, and where the derivatives also give rise to closed-form expressions of the Hessian to guarantee positive semi-definiteness for a fast and concise Newton-type implicit time integration.
A consequence of this analysis is a novel analytical formulation to compute rotations and singular values of
deformation-gradient tensors without explicit/numerical factorization, which is significant, resulting in up to 3.5x speedup and benefits energy function evaluation for reducing solver time.
We validate our energy formulation by experiments and comparison, demonstrating that our resulting eigendecomposition using the CG invariants is equivalent to existing ARAP formulations.
We thus reveal isotropic ARAP energy to be a member of the "Cauchy-Green club", meaning that it can indeed be defined using CG invariants and therefore that the closed-form expressions of the resulting Hessian are shared with other energies written in their terms.


[Paper]: https://drive.google.com/file/d/1iJ7XS7T8d9ViS-nuCF0BX-M146Zd2vG4/view?usp=share_link

[video]: https://drive.google.com/file/d/1pP31gomGFFMi9U8qM5d9JXPKNB-uih8p/view?usp=share_link

[code]: https://github.com/LamWS/iARAP

[supplementary]: https://drive.google.com/file/d/1-xbcU3GJ1DcDBt3g9MF8tfdFL61Y3BUA/view?usp=share_link