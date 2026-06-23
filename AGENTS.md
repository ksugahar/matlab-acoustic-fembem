# Repository Policy

## Purpose

This is a local, non-public MATLAB education repository for Gypsilab-style
FEM/BEM prototypes. It prioritizes readability over performance.

## Solver Scope

- Keep the API short and inspectable, in the spirit of Gypsilab.
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

