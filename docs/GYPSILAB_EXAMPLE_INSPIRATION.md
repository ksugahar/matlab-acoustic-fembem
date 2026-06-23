# Gypsilab Example Inspiration

The first 100 examples in this repository are organized around the example
families that made Gypsilab useful as readable MATLAB source.

This repository does not copy Gypsilab examples verbatim. It absorbs the
teaching structure and rewrites examples around the local `.vol` mesh lane,
radia-ngsolve validation, and readable first-order FEM/BEM code.

## Source Families

- `openMsh`: mesh construction, cleanup, boundary extraction, refinement, IO
- `openFem`: Lagrange, Nedelec, RWG, DOF maps, FEM operators
- `openOpr`: Green kernels, BEM operator builders, ACA
- `openHmx`: H-matrix blocks, low-rank compression, algebra, builder examples
- `openEbd`: 2D kernels, Helmholtz/Laplace scalar products, radial quadrature
- `miscellaneous`: `sphereHelmholtz`, `diskHelmholtz`, `sphereMaxwell`,
  vibro-acoustic sketches
- `nonRegressionTest/finiteElement`: FEM convergence, wave, Dirichlet examples
- `nonRegressionTest/operators`: Helmholtz, Maxwell, Stokes operator examples
- `nonRegressionTest/scattering2d` and `scattering3d`: acoustic and EM
  scattering examples, including H-matrix variants
- `nonRegressionTest/radiationImpedances`: baffled/unbaffled acoustic disk
  radiation impedance
- `nonRegressionTest/vibroAcoustic`: FEM/BEM vibro-acoustic coupling examples
- `nonRegressionTest/femBemDielectrique`: RWG/Nedelec-style FEM/BEM dielectric
  coupling examples
- `nonRegressionTest/rayTracing`: geometric acoustics/ray-tracing examples

## Local Rule

Each local example should be short enough for a student to open first. The
script may delegate to a readable implementation, but the case id, topic, and
Gypsilab inspiration must be visible immediately.

Status after initial example pass:

- all 100 case scripts exist under `examples`
- 10 mesh/topology scripts are verified
- 90 scripts are planned teaching cards awaiting radia-ngsolve promotion
