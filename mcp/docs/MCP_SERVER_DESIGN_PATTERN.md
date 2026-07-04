# MCP Design Pattern

This repository uses the official MathWorks MATLAB MCP Server as the
runtime/transport model and keeps acoustic FEM/BEM behavior as tested MATLAB
functions plus extension JSON.

## Pattern

- Prefer official custom-tool extension files before writing another server.
- Tool JSON contains the MCP-facing name, description, schema, and annotations.
- MATLAB package functions in `+acoustic_fembem` contain the domain behavior.
- Complex inputs are passed by file path when the official custom-tool
  extension needs scalar JSON argument types.
- Verification may use an existing MATLAB Engine session or `matlab -batch`.

## Current Tools

- `acoustic_fembem_check_result_manifest_file`
- `acoustic_fembem_acoustic_gate`
- `acoustic_fembem_crossval_gate`
- `acoustic_fembem_knowledge`
- `acoustic_fembem_repository_health`

## Add A Tool

1. Add a MATLAB package function under `+acoustic_fembem`.
2. Add or update a test under `tests`.
3. Register the entry point in `mcp/extensions/acoustic-fembem-tools.json`.
4. Run `mcp/tools/verify_acoustic_fembem_mcp.ps1`.
5. Keep file-writing helpers, such as PDE Toolbox mesh-to-.vol export, as normal
   MATLAB APIs unless the MCP tool explicitly needs write access.
