# COMSOL Acoustic FEM/BEM Internal Reference

COMSOL's acoustic FEM/BEM workflow is an important internal reference for this
MATLAB education layer.

The primary public-safe validation target remains radia-ngsolve / NGSolve.BEM.
COMSOL is used privately to understand practical FEM/BEM modeling choices:

- acoustic pressure FEM domains
- exterior BEM radiation treatment
- FEM/BEM coupling signs and normal conventions
- low-frequency behavior near the Laplace limit
- scattering and radiation benchmark setup
- how production software presents acoustic open-boundary workflows

## Policy

- Use only the already-running COMSOL LiveLink MATLAB session.
- Do not start `comsolmphserver`, `comsolbatch`, or a second MATLAB session.
- Do not kill MATLAB or COMSOL processes.
- Keep COMSOL-derived values and model files private.
- Do not copy COMSOL implementation details into public repos.

## Role in the 100-Case Campaign

For acoustic examples, a case can have:

- required reference: radia-ngsolve / NGSolve.BEM
- optional internal reference: COMSOL acoustic FEM/BEM

The COMSOL reference is especially valuable for:

- Helmholtz exterior scattering
- radiation/open-boundary comparisons
- FEM pressure domain coupled to exterior BEM
- low-frequency stabilization sanity checks
- source-normal derivative sign conventions

The example is not considered verified unless the radia-ngsolve comparison
passes. A COMSOL run strengthens the case, but does not replace the
radia-ngsolve gate.

