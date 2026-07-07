# MathWorks MATLAB MCP Notes

Use the official MathWorks MATLAB MCP Server as the runtime.  Do not vendor or
patch the server binary in this repository.

Official resources:

- `https://github.com/matlab/matlab-mcp-server`
- `https://github.com/matlab/matlab-mcp-server/releases`
- `https://github.com/matlab/matlab-agentic-toolkit`
- `https://www.mathworks.com/products/matlab-mcp-server.html`

This repository contributes only the acoustic FEM/BEM domain layer:

- extension metadata in `mcp/extensions/acoustic-fembem-tools.json`
- MATLAB entry points in `+acoustic_fembem`
- focused tests in `tests/testMcpAcousticFembemTools.m`

Keep MCP functions as stable entry points around normal MATLAB functions.  The
engineering behavior should remain testable without an MCP client.

## Layering

The official MATLAB MCP Server is the runtime.  It owns MATLAB session
management, code execution, test execution, Code Analyzer checks, and toolbox
detection.

The MATLAB Agentic Toolkit is the setup and skills layer.  Use it to install or
update the official server and to provide MATLAB workflow guidance to the agent.
Do not copy the server binary or toolkit implementation into this repository.

This repository is the domain extension.  It provides only the acoustic
FEM/BEM/CQ entry points, validation gates, and extension metadata needed by the
official server:

- `.vol` tri/tet preflight and summaries
- P1 volume FEM and P1 boundary BEM gates
- Johnson-Nedelec / Calderon coupling checks
- Lubich CQ acoustic transient teaching gates
- NGSolve.BEM / radia-ngsolve cross-validation handoff metadata

When COMSOL LiveLink is involved, use the official server's existing-session
workflow against the MATLAB session that is already connected to COMSOL.  Pure
MATLAB acoustic FEM/BEM checks may run separately when they do not touch
COMSOL.
