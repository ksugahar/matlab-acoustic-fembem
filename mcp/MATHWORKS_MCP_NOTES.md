# MathWorks MATLAB MCP Notes

Use the official MathWorks MATLAB MCP Server as the runtime.  Do not vendor or
patch the server binary in this repository.

Official resources:

- `https://github.com/matlab/matlab-mcp-server`
- `https://github.com/matlab/matlab-mcp-server/releases`
- `https://www.mathworks.com/products/matlab-mcp-server.html`

This repository contributes only the acoustic FEM/BEM domain layer:

- extension metadata in `mcp/extensions/acoustic-fembem-tools.json`
- MATLAB entry points in `+acoustic_fembem`
- focused tests in `tests/testMcpAcousticFembemTools.m`

Keep MCP functions as stable entry points around normal MATLAB functions.  The
engineering behavior should remain testable without an MCP client.
