function repoRoot = gypsilabRepoRoot()
%GYPSILABREPOROOT Return this repository's root directory.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
end
