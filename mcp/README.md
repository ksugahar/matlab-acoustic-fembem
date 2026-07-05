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
| `acoustic_fembem_repository_health` | `acoustic_fembem.check_repository_health` | Pre-push health check for the integrated repo and MCP extension |

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

The official server supplies the MCP runtime.  This repository supplies the
domain tools and validation gates.

## Repository Boundary

The acoustic FEM/BEM MCP layer stays inside this repository.  Do not maintain a
separate lab `matlab-mcp` repository for these tools: the extension file is the
public contract, while the official MathWorks server remains the runtime.
Internal cross-validation provenance, local solver logs, and private MATLAB
automation notes stay outside the public package unless rewritten as scrubbed
domain knowledge.
