# Repository Policy

## Purpose

This is a local, non-public MATLAB education repository for Sugawara Lab,
which is treated as a CAE-AI Lab. The repository teaches Gypsilab-style
FEM/BEM prototypes and MATLAB optimization methods. It prioritizes readability,
cross validation, and student understanding over performance.

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
- Treat COMSOL acoustic FEM/BEM as an important internal reference for acoustic
  coupling, sign conventions, radiation boundaries, and low-frequency behavior.
  Use only the existing COMSOL LiveLink MATLAB session.
- Do not optimize prematurely. NGSolve/NGSolve.BEM is the performance target;
  this repo is for understanding the method.
- Prefer readable class/source organization over vectorized cleverness or
  production-style performance plumbing.

## Optimization Scope

- Learn MATLAB optimization alongside Gypsilab-style FEM/BEM.
- Prefer small, inspectable optimization examples with explicit design
  variables, objective functions, constraints, sensitivities when available,
  and reproducible validation gates.
- Treat `W:\00_CAE\MATLAB\30_Optimization` as the main learning source for
  optimization methods, but distill only private-safe, license-safe teaching
  code into this repository.
- Compare MATLAB optimization examples with Optuna when useful. Use the
  installed official public `optuna/optuna-mcp` or plain Optuna APIs; do not
  create a lab-specific Optuna MCP server.
- Couple optimization examples to readable FEM/BEM operators where useful:
  shape/mesh parameters, material coefficients, acoustic objectives, and
  scalar/vector field targets are good teaching cases.
- Do not hide optimization logic behind opaque solver wrappers. Students
  should be able to open the `.m` file and see the algorithmic loop.

## MATLAB / COMSOL Session

- On INTEL11, use the existing COMSOL LiveLink MATLAB session for MATLAB work
  when available.
- Do not start or kill COMSOL server, COMSOL batch, or the existing MATLAB
  session.
- Use `S:\COMSOL\mcp-server\livelink\scripts\eval_shared_matlab.py` if the
  MATLAB MCP cannot discover the running session.

## Public Boundary

- This repository is private/local.
- Do not push it to public GitHub or PyPI.
- Public-safe generic helpers may be mirrored into radia-mcp only when they do
  not mention private MATLAB/Gypsilab/Lukas provenance or commercial tools.
