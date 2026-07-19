# Acoustic FEM/BEM MCP Layer

This directory exposes the MATLAB acoustic FEM/BEM teaching code through the
official MathWorks MATLAB MCP Server.

The solver, validation fixtures, and MCP entry points live in the same
repository:

- MATLAB package: `+acoustic_fembem`
- Extension file: `mcp/extensions/acoustic-fembem-tools.json`
- Verifier: `mcp/tools/verify_acoustic_fembem_mcp.ps1`

## Tools

| Tool | MATLAB function | Purpose |
|------|-----------------|---------|
| `acoustic_fembem_check_result_manifest_file` | `acoustic_fembem.check_result_manifest_file` | Validate script/MCP-ready result manifests |
| `acoustic_fembem_acoustic_gate` | `acoustic_fembem.check_fembem_acoustic_gate` | Run acoustic FEM/BEM gates against analytic references |
| `acoustic_fembem_crossval_gate` | `acoustic_fembem.check_fembem_crossval_gate` | Run `.vol`-backed cross-validation against radia-ngsolve/NGSolve references |
| `acoustic_fembem_knowledge` | `acoustic_fembem.fembem_knowledge_tool` | Serve compact acoustic FEM/BEM teaching knowledge |
| `acoustic_fembem_vol_mesh_summary` | `acoustic_fembem.check_vol_mesh_summary` | Summarize `.vol` meshes and viewer guidance |
| `acoustic_fembem_cq_time_grid` | `acoustic_fembem.check_cq_time_grid` | Validate CQ time and Laplace grids |
| `acoustic_fembem_cq_response_reality` | `acoustic_fembem.check_cq_response_reality` | Check CQ conjugate symmetry and real response |
| `acoustic_fembem_soft_sphere_cq_causality` | `acoustic_fembem.check_soft_sphere_cq_causality` | Gate causal soft-sphere CQ behavior |
| `acoustic_fembem_adjoint_scaling` | `acoustic_fembem.check_adjoint_scaling` | Gate acoustic adjoint scaling |
| `acoustic_fembem_hmatrix_scaling` | `acoustic_fembem.check_hmatrix_scaling` | Run readable ACA+ scaling and gate dense-reference error, rank, and storage growth |

The acoustic gate includes adjoint optimization cases:

- `focus_adjoint`: gradient of focused pressure intensity vs finite difference.
- `radiation_force`: acoustic radiation-force postprocess from the BEM field.
- `thrust_adjoint`: elastic-bead force quadratic form and Wirtinger gradient.

The knowledge tool includes `matlab_execution_policy`, `vol_visualization`, and
`pde_vol_bridge`.  The policy is normal `.m` functions/scripts plus MCP JSON,
Netgen for native interactive `.vol` inspection, MATLAB `plotVolMesh` for
figure previews, and `acoustic_fembem_vol_mesh_summary` for LLM/headless
preflight.

## Official Server

Install the official MathWorks MATLAB MCP Server separately, then pass:

```powershell
--extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

For an already shared R2026a session, preflight that exactly one shared
session is visible, then use the official server's capture route:

```powershell
--matlab-session-mode=auto --matlab-display-mode=nodesktop `
  --extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

The official CLI does not allow `existing` together with `nodesktop`.  The
`auto` route is acceptable here only after the single-session preflight.  The
MCP log must say `Attaching to existing session`, and the MATLAB process set
must be unchanged before and after the call.  This keeps command-window output
on the capture `FEval` path without closing or replacing the shared desktop
session.

The official server supplies the MCP runtime.  This repository supplies the
domain tools and validation gates.

## MATLAB Agentic Toolkit

Use the MathWorks MATLAB Agentic Toolkit as the setup and skills reference
layer.  It can install/update the official MATLAB MCP Server and register
agent skills for MATLAB workflows.  This repository does not replace those
pieces; it adds only the acoustic FEM/BEM/CQ domain extension.

Recommended split:

- official MATLAB MCP Server: MATLAB session/runtime, code execution, tests,
  Code Analyzer checks, and toolbox detection;
- MATLAB Agentic Toolkit: skills and setup guidance, selected narrowly for the
  project;
- `matlab-acoustic-fembem`: `.vol` tri/tet intake, P1 FEM/BEM, Johnson-Nedelec
  coupling, Lubich CQ, Gmsh artifacts, and NGSolve.BEM/radia cross-validation
  gates.

When an external solver bridge is involved, attach to the already shared MATLAB
session through the official server's existing-session workflow. Pure MATLAB
acoustic FEM/BEM checks may run separately when they do not touch that bridge.

## Repository Boundary

The acoustic FEM/BEM MCP layer stays inside this repository.  Do not maintain a
separate lab `matlab-mcp` repository for these tools: the extension file is the
public contract, while the official MathWorks server remains the runtime.
Internal cross-validation provenance, local solver logs, and private MATLAB
automation notes stay outside the public package unless rewritten as scrubbed
domain knowledge.
