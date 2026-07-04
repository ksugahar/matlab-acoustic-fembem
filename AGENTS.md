# Repository Policy

## Purpose

This is a public MATLAB education repository for readable FEM/BEM prototypes
and small optimization gates. It prioritizes readability, cross validation,
and student understanding over performance.

## Solver Scope

- Keep the API short and inspectable, in the spirit of Gypsilab.
- Follow Gypsilab's readable class style when introducing MATLAB classes:
  properties should expose the mathematical objects, methods should be short,
  and performance caches should not hide the algorithm.
- Use Cubit/Coreform Netgen `.vol` as the mesh intake.
- Accept only first-order triangle surface elements and tetrahedron volume
  elements.
- First-order spaces only for now: H1 P1, HCurl/Nedelec0, scalar BEM P1, and
  RWG0.
- Include acoustics as a first-class teaching lane.
- Include low-frequency-stable Helmholtz kernels by splitting Laplace and
  Helmholtz correction terms.
- Do not optimize prematurely. NGSolve/NGSolve.BEM is the performance target;
  this repo is for understanding the method.
- Prefer readable class/source organization over vectorized cleverness or
  production-style performance plumbing.

## Optimization Scope

- Learn MATLAB optimization alongside Gypsilab-style FEM/BEM.
- Prefer small, inspectable optimization examples with explicit design
  variables, objective functions, constraints, sensitivities when available,
  and reproducible validation gates.
- Do not maintain Optuna in this repository. When Optuna operation is needed,
  use the official `optuna-mcp-server`.
- Couple optimization examples to readable FEM/BEM operators where useful:
  shape/mesh parameters, material coefficients, acoustic objectives, and
  scalar/vector field targets are good teaching cases.
- Do not hide optimization logic behind opaque solver wrappers. Students
  should be able to open the `.m` file and see the algorithmic loop.

## License Boundary

- This repository is GPL-3.0-or-later.
- Do not vendor third-party source trees. Optional comparisons against
  separately installed tools should remain explicit verification helpers.
- Keep examples reproducible from committed fixtures, analytic references, or
  openly available numerical packages.
