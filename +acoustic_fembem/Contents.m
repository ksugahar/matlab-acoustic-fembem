% ACOUSTIC_FEMBEM  MCP-facing helpers for MATLAB acoustic FEM/BEM.
%
%   check_result_manifest_file  MCP custom-tool entry point for result manifests.
%   result_manifest_gate  Check notebook-ready result JSON/struct metadata,
%                         timing, versions, parameter sets, and objectives.
%
%   Integrated acoustic FEM/BEM (Gypsilab-readable lane):
%   repository_root          Resolve this repository root.
%   fembem_knowledge          Queryable FEM/BEM knowledge (topic dispatcher).
%   fembem_knowledge_tool     MCP custom-tool entry point: print one topic.
%   repository_health         Lightweight repository and MCP extension health.
%   check_repository_health   MCP custom-tool entry point: print JSON verdict.
%   vol_mesh_summary          Summarize .vol mesh metadata and viewer guidance.
%   check_vol_mesh_summary    MCP custom-tool entry point: print JSON summary.
%   fembem_acoustic_gate      Run a Gypsilab acoustic solve vs the analytic
%                             series; return a struct verdict.
%   check_fembem_acoustic_gate  MCP custom-tool entry point: print JSON verdict.
%   fembem_crossval_gate      Run Gypsilab FEM/BEM vs radia-ngsolve/NGSolve.
%   check_fembem_crossval_gate  MCP custom-tool entry point: print JSON verdict.
