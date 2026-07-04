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
| `acoustic_fembem_check_result_manifest_file` | `acoustic_fembem.check_result_manifest_file` | Validate notebook-ready result manifests |
| `acoustic_fembem_acoustic_gate` | `acoustic_fembem.check_fembem_acoustic_gate` | Run acoustic FEM/BEM gates against analytic references |
| `acoustic_fembem_crossval_gate` | `acoustic_fembem.check_fembem_crossval_gate` | Run `.vol`-backed cross-validation against radia-ngsolve/NGSolve references |
| `acoustic_fembem_knowledge` | `acoustic_fembem.fembem_knowledge_tool` | Serve compact acoustic FEM/BEM teaching knowledge |
| `acoustic_fembem_repository_health` | `acoustic_fembem.check_repository_health` | Pre-push health check for the integrated repo and MCP extension |

The acoustic gate includes adjoint optimization cases:

- `focus_adjoint`: gradient of focused pressure intensity vs finite difference.
- `radiation_force`: acoustic radiation-force postprocess from the BEM field.
- `thrust_adjoint`: elastic-bead force quadratic form and Wirtinger gradient.

The knowledge tool includes `pde_vol_bridge`, which records the MATLAB PDE
Toolbox `generateMesh` -> first-order Netgen `.vol` path for simple geometries.

## Official Server

Install the official MathWorks MATLAB MCP Server separately, then pass:

```powershell
--extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

The official server supplies the MCP runtime.  This repository supplies the
domain tools and validation gates.
