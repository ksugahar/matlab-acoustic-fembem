function root = repository_root()
%REPOSITORY_ROOT Resolve this MATLAB acoustic FEM/BEM repository.
%
%   root = acoustic_fembem.repository_root();
%
% The MCP layer is shipped inside this repository.  This resolver keeps the
% MCP entry points independent of the current working directory.
%
% Resolution order (fail-loud, no silent fallback):
%   1. the ACOUSTIC_FEMBEM_ROOT environment variable, if set
%   2. this package's parent repository
% A directory qualifies only if it contains matlab_api/; otherwise error.

candidates = strings(1, 0);
envRoot = getenv("ACOUSTIC_FEMBEM_ROOT");
if strlength(envRoot) > 0
    candidates(end + 1) = string(envRoot);
end
repoRoot = fileparts(fileparts(mfilename("fullpath")));   % +acoustic_fembem -> repo root
candidates(end + 1) = string(repoRoot);

for c = candidates
    if isfolder(fullfile(c, "matlab_api"))
        root = c;
        return
    end
end

error("acoustic_fembem:repositoryRootNotFound", ...
    "MATLAB acoustic FEM/BEM repo not found (needs matlab_api/). Set " + ...
    "ACOUSTIC_FEMBEM_ROOT or run from the integrated repository. Tried: %s", ...
    strjoin(candidates, "; "));
end
