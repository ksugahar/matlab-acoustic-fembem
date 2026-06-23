function writeExampleScripts()
%WRITEEXAMPLESCRIPTS Generate thin example scripts for all catalog entries.
%
% The generated scripts are intentionally tiny. Students open the script to see
% the case id, then follow `runExampleById` into the readable implementation.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "examples"));
cases = validationCatalog();

for k = 1:numel(cases)
    path = fullfile(repoRoot, cases(k).examplePath);
    folder = fileparts(path);
    if ~isfolder(folder)
        mkdir(folder);
    end
    if isfile(path) && cases(k).status == "verified"
        continue
    end
    lines = [
        "% " + cases(k).id + " " + cases(k).title
        "% Gypsilab inspiration: " + cases(k).gypsilabInspiration
        "repoRoot = fileparts(fileparts(fileparts(mfilename(""fullpath""))));"
        "addpath(fullfile(repoRoot, ""examples""));"
        "result = runExampleById(""" + cases(k).id + """);"
    ];
    fid = fopen(path, "w");
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "%s\n", lines);
    clear cleanup
end
end
