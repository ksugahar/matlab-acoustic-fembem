# MCP Requirements

This repository is an acoustic FEM/BEM MATLAB package with an optional MCP
surface for the official MathWorks MATLAB MCP Server.

## Required

- MATLAB R2021a or later.
- Official MathWorks MATLAB MCP Server:
  `https://github.com/matlab/matlab-mcp-server`.
- An MCP client that can launch a stdio server.

## Recommended

- MATLAB Agentic Toolkit:
  `https://github.com/matlab/matlab-agentic-toolkit`.

The toolkit is the preferred way to install/update the official server and add
MATLAB workflow skills.  Install only the skill groups that match this project;
for this repository the useful groups are MATLAB Core, MATLAB Programming,
MATLAB Software Development, and Math/Optimization/PDE-style workflows.

## Installing The Official Server

Use the official release binary for Windows/Linux, or build from source with
Go:

```powershell
go install github.com/matlab/matlab-mcp-server/cmd/matlab-mcp-server@latest
```

For an existing MATLAB session workflow, install the MATLAB-side toolbox with
the official server and run `shareMATLABSession()` in the MATLAB session you
want the MCP server to attach to.

## Extension File

After installing the official server, pass this repository's extension file:

```powershell
matlab-mcp-server-windows-x64.exe --extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

Add `--matlab-session-mode=existing` when you want the server to attach to a
shared MATLAB session instead of starting a new one.

The acoustic FEM/BEM extension is intentionally thin.  The official MathWorks
server owns process/session management; this repository owns only domain MATLAB
functions, validation gates, and the extension-file contract.

## Existing MATLAB Session Policy

For COMSOL LiveLink work, attach the official server to the already-running
MATLAB session that shares the COMSOL server.  Do not start a second MATLAB for
COMSOL-bound checks.  Batch MATLAB is acceptable for pure MATLAB acoustic
FEM/BEM tests that do not touch COMSOL.
