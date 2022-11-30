# Isotropic-ARAP-energy-using-Cauchy-Green-invariants
Paper Isotropic ARAP energy using Cauchy-Green invariants

![](https://github.com/LamWS/iARAP/raw/main/file/teaser.png)

## Abstract

Many strategies exist for optimizing non-linear distortion energies in geometry and physics applications, but devising an approach that achieves the convergence promised by Newton-type methods remains challenging. In order to guarantee the positive semi-definiteness required by these methods, a numerical eigendecomposition or approximate regularization is usually needed. In this article, we present analytic expressions for the eigensystems at each quadrature point of a wide range of isotropic distortion energies. These systems can then be used to project energy Hessians to positive semi-definiteness analytically. Unlike previous attempts, our formulation provides compact expressions that are valid both in 2D and 3D, and does not introduce spurious degeneracies. At its core, our approach utilizes the invariants of the stretch tensor that arises from the polar decomposition of the deformation gradient. We provide closed-form expressions for the eigensystems for all these invariants, and use them to systematically derive the eigensystems of any isotropic energy. Our results are suitable for geometry optimization over flat surfaces or volumes, and agnostic to both the choice of discretization and basis function. To demonstrate the efficiency of our approach, we include comparisons against existing methods on common graphics tasks such as surface parameterization and volume deformation.

[Paper], [Video], and [Code] are now available.

[Paper]: https://drive.google.com/file/d/1iJ7XS7T8d9ViS-nuCF0BX-M146Zd2vG4/view?usp=share_link
[Video]: https://drive.google.com/file/d/1pP31gomGFFMi9U8qM5d9JXPKNB-uih8p/view?usp=share_link
[Code]: https://github.com/LamWS/iARAP