function report = repository_health(options)
%REPOSITORY_HEALTH Check the public acoustic FEM/BEM MCP checkout.
%
%   report = acoustic_fembem.repository_health()
%
% This is a lightweight pre-push and MCP smoke gate.  It verifies that the
% integrated repository has the expected acoustic name, package paths,
% extension-tool registrations, .vol fixtures, and 100-case teaching catalog.

arguments
    options.ScanOldNames (1,1) logical = true
end

root = acoustic_fembem.repository_root();
addpath(fullfile(root, "examples"));

extensionPath = fullfile(root, "mcp", "extensions", "acoustic-fembem-tools.json");
toolNames = strings(1, 0);
extensionParsed = false;
if isfile(extensionPath)
    try
        extension = jsondecode(fileread(extensionPath));
        toolNames = string({extension.tools.name});
        extensionParsed = true;
    catch
        toolNames = strings(1, 0);
    end
end

requiredTools = [
    "acoustic_fembem_check_result_manifest_file"
    "acoustic_fembem_acoustic_gate"
    "acoustic_fembem_crossval_gate"
    "acoustic_fembem_knowledge"
    "acoustic_fembem_repository_health"
    "acoustic_fembem_vol_mesh_summary"
];

catalogCount = NaN;
verifiedCount = NaN;
try
    cases = validationCatalog();
    catalogCount = numel(cases);
    verifiedCount = nnz(string({cases.status}) == "verified");
catch
end

volFixtures = dir(fullfile(root, "fixtures", "mesh_topology", "*.vol"));
oldNameHits = struct("file", {}, "pattern", {});
if options.ScanOldNames
    oldNameHits = scanOldNames(root);
end
mcpUiPolicyHits = scanMcpUiPolicy(root);

requiredPaths = struct( ...
    "matlab_api", isfolder(fullfile(root, "matlab_api")), ...
    "fixtures", isfolder(fullfile(root, "fixtures", "mesh_topology")), ...
    "validation", isfolder(fullfile(root, "validation")), ...
    "mcp", isfolder(fullfile(root, "mcp")), ...
    "package", isfolder(fullfile(root, "+acoustic_fembem")), ...
    "extension_file", isfile(extensionPath));

checks = struct();
checks.required_paths_present = all(structfun(@(x) logical(x), requiredPaths));
checks.extension_json_parsed = extensionParsed;
checks.required_tools_registered = all(ismember(requiredTools, toolNames));
checks.catalog_has_100_cases = catalogCount == 100;
checks.catalog_verified = verifiedCount == 100;
checks.vol_fixture_count_sufficient = numel(volFixtures) >= 10;
checks.pde_vol_bridge_present = isfile(fullfile(root, "matlab_api", "mesh", "writePdeMeshVol.m"));
checks.vol_plot_preview_present = isfile(fullfile(root, "matlab_api", "mesh", "plotVolMesh.m"));
checks.old_repository_names_absent = isempty(oldNameHits);
checks.mcp_surface_avoids_live_document_policy = isempty(mcpUiPolicyHits);

failed = failedCheckNames(checks);

report = struct();
report.tool = "acoustic_fembem_repository_health";
report.repository_name = "matlab-acoustic-fembem";
report.repository_url = "https://github.com/ksugahar/matlab-acoustic-fembem";
report.root = string(root);
report.extension_file = string(extensionPath);
report.status = "needs_attention";
report.pass = isempty(failed);
report.checks = checks;
report.failed_checks = failed;
report.required_paths = requiredPaths;
report.required_tools = requiredTools;
report.registered_tools = toolNames(:);
report.num_validation_cases = catalogCount;
report.num_verified_cases = verifiedCount;
report.num_vol_fixtures = numel(volFixtures);
report.old_name_hits = oldNameHits;
report.mcp_ui_policy_hits = mcpUiPolicyHits;

if report.pass
    report.status = "ok";
end
end


function hits = scanMcpUiPolicy(root)
patterns = [
    "note" + "book"
    "ip" + "ynb"
    "m" + "lx"
    "Live" + " Script"
];

files = [
    string(fullfile(root, "README.md"))
    textFilesUnder(fullfile(root, "+acoustic_fembem"))
    textFilesUnder(fullfile(root, "mcp"))
];

hits = scanTextPatterns(root, files, patterns, ...
    [string(fullfile(root, "+acoustic_fembem", "repository_health.m"))]);
end


function hits = scanOldNames(root)
patterns = [
    "caeai-" + "matlab-mcp"
    "caeai-" + "matlab-fembem"
    "caeai" + "_check_"
    "caeai" + "_fembem_"
    "cae-ai" + "-lab"
    "GYPSILAB" + "_ROOT"
    "gypsilab" + "_root"
];

files = [
    string(fullfile(root, "README.md"))
    string(fullfile(root, "run_tests.m"))
    textFilesUnder(fullfile(root, "+acoustic_fembem"))
    textFilesUnder(fullfile(root, "mcp"))
    textFilesUnder(fullfile(root, "docs"))
    textFilesUnder(fullfile(root, "matlab_api"))
    textFilesUnder(fullfile(root, "tests"))
];

hits = scanTextPatterns(root, files, patterns, ...
    [string(fullfile(root, "+acoustic_fembem", "repository_health.m"))]);
end


function hits = scanTextPatterns(root, files, patterns, excludedFiles)
hits = struct("file", {}, "pattern", {});
excludedFiles = string(excludedFiles);
for k = 1:numel(files)
    file = files(k);
    if ~isfile(file)
        continue
    end
    if any(strcmpi(file, excludedFiles))
        continue
    end
    try
        body = string(fileread(file));
    catch
        continue
    end
    for p = 1:numel(patterns)
        pattern = patterns(p);
        if contains(body, pattern)
            hits(end + 1).file = relativePath(root, file); %#ok<AGROW>
            hits(end).pattern = pattern;
        end
    end
end
end


function files = textFilesUnder(folder)
allowedExtensions = [".m", ".md", ".json", ".ps1"];
entries = dir(fullfile(folder, "**", "*"));
files = strings(0, 1);
for k = 1:numel(entries)
    entry = entries(k);
    if entry.isdir
        continue
    end
    [~, ~, ext] = fileparts(entry.name);
    if any(strcmpi(ext, allowedExtensions))
        files(end + 1, 1) = string(fullfile(entry.folder, entry.name)); %#ok<AGROW>
    end
end
end


function rel = relativePath(root, file)
prefix = string(root) + filesep;
rel = erase(string(file), prefix);
end


function names = failedCheckNames(checks)
names = strings(1, 0);
checkNames = string(fieldnames(checks));
for k = 1:numel(checkNames)
    checkName = checkNames(k);
    if ~checks.(checkName)
        names(end + 1) = checkName; %#ok<AGROW>
    end
end
end
