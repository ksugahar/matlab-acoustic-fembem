function cases = validationCatalog()
%VALIDATIONCATALOG Planned 100-case radia-ngsolve validation catalog.
%
% Cases start as planned. Promote to verified only after a reproducible
% MATLAB example and radia-ngsolve reference agree within the stated tolerance.

groups = [
    group("01_mesh_topology", "mesh", 1e-14, [
        "unit tetra vol topology"
        "four tet interior node trace"
        "closed tetra surface manifold"
        "two material tetra labels"
        "boundary name propagation"
        "reversed surface orientation signs"
        "tet edge orientation signs"
        "quad surface rejection"
        "hex volume rejection"
        "real cubit cylinder vol intake"
    ])
    group("02_h1_scalar_fem", "h1", 1e-10, [
        "unit tetra p1 stiffness"
        "unit tetra p1 mass"
        "constant gradient patch"
        "linear manufactured potential"
        "dirichlet elimination"
        "neumann face load"
        "two tet patch continuity"
        "material coefficient scaling"
        "volume residual balance"
        "ngsolve p1 scalar comparison"
    ])
    group("03_hcurl_edge_fem", "hcurl", 1e-10, [
        "single tet nedelec edge count"
        "single tet edge signs"
        "two tet shared edge continuity"
        "curl basis orientation"
        "mass matrix symmetry"
        "curl curl matrix symmetry"
        "gradient nullspace trace"
        "boundary edge extraction"
        "material mu scaling"
        "ngsolve hcurl comparison"
    ])
    group("04_laplace_bem_dense", "laplace_dense", 1e-10, [
        "single triangle centroid kernel"
        "two separated triangles"
        "laplace symmetry weighted"
        "near field diagonal policy"
        "sphere coarse capacitance trend"
        "cube coarse capacitance trend"
        "constant density potential"
        "normal derivative sign"
        "dense matvec reproducibility"
        "ngsolve bem laplace comparison"
    ])
    group("05_laplace_hmatrix", "laplace_hmatrix", 1e-8, [
        "cluster tree split"
        "admissibility threshold"
        "low rank far block"
        "dense near block"
        "hmatrix matvec vs dense"
        "rank tolerance sweep"
        "leaf size sweep"
        "storage ratio report"
        "sphere hmatrix laplace"
        "ngsolve bem hmatrix comparison"
    ])
    group("06_acoustic_low_frequency", "acoustic_low_frequency", 1e-9, [
        "helmholtz k zero laplace limit"
        "single layer expm1 correction"
        "single layer taylor correction"
        "double layer source normal correction"
        "small kr sweep"
        "complex symmetry weighted"
        "low frequency matvec"
        "low frequency sphere monopole"
        "low frequency hmatrix candidate"
        "ngsolve bem low frequency comparison"
    ])
    group("07_acoustic_helmholtz", "acoustic_helmholtz", 1e-8, [
        "two point helmholtz kernel"
        "plane wave boundary data"
        "pulsating sphere low ka"
        "pulsating sphere mid ka"
        "rigid sphere scattering small ka"
        "dirichlet sphere scattering"
        "frequency sweep phase"
        "complex matvec reproducibility"
        "helmholtz hmatrix candidate"
        "ngsolve bem helmholtz comparison"
    ])
    group("08_scalar_fem_bem_coupling", "scalar_coupling", 1e-8, [
        "trace matrix unit tetra"
        "trace matrix interior node"
        "dirichlet to neumann toy"
        "laplace open boundary toy"
        "flux balance check"
        "single material coupling"
        "two material coupling"
        "sphere exterior scalar potential"
        "scalar fem bem hmatrix"
        "ngsolve coupled scalar comparison"
    ])
    group("09_rwg_hcurl_trace", "rwg_hcurl", 1e-10, [
        "rwg closed manifold dofs"
        "rwg open edge rejection"
        "rwg opposite vertex map"
        "rwg to hcurl edge ids"
        "boundary edge orientation"
        "two tet shared boundary"
        "surface current toy"
        "magnetic vector trace toy"
        "edge trace sparse map"
        "ngsolve hcurl trace comparison"
    ])
    group("10_ngsolve_bem_reference", "ngsolve_bem", 1e-8, [
        "ngsolve mesh import smoke"
        "ngsolve p1 scalar smoke"
        "ngsolve hcurl smoke"
        "ngsolve laplace bem smoke"
        "ngsolve helmholtz bem smoke"
        "ngsolve low frequency smoke"
        "ngsolve sphere analytic laplace"
        "ngsolve sphere analytic acoustic"
        "ngsolve fem bem coupling smoke"
        "full pipeline capstone"
    ])
];

cases = repmat(emptyCase(), 100, 1);
caseIndex = 0;
for g = 1:numel(groups)
    for k = 1:numel(groups(g).titles)
        caseIndex = caseIndex + 1;
        title = groups(g).titles(k);
        id = sprintf("GYP-%03d", caseIndex);
        cases(caseIndex, 1) = struct( ...
            "id", string(id), ...
            "category", groups(g).category, ...
            "shortName", groups(g).shortName + "_" + sprintf("%02d", k), ...
            "title", title, ...
            "status", "planned", ...
            "reference", "radia-ngsolve", ...
            "secondaryReference", secondaryReferenceFor(groups(g).category), ...
            "gypsilabInspiration", gypsilabInspirationFor(groups(g).category, title), ...
            "tolerance", groups(g).tolerance, ...
            "examplePath", fullfile("examples", groups(g).category, lower(id) + "_" + slug(title) + ".m"), ...
            "validationLog", "");
    end
end
cases = markMeshTopologyVerified(cases);
end


function g = group(category, shortName, tolerance, titles)
g = struct();
g.category = string(category);
g.shortName = string(shortName);
g.tolerance = tolerance;
g.titles = string(titles);
end


function c = emptyCase()
c = struct( ...
    "id", "", ...
    "category", "", ...
    "shortName", "", ...
    "title", "", ...
    "status", "", ...
    "reference", "", ...
    "secondaryReference", "", ...
    "gypsilabInspiration", "", ...
    "tolerance", NaN, ...
    "examplePath", "", ...
    "validationLog", "");
end


function reference = secondaryReferenceFor(category)
if ismember(category, ["06_acoustic_low_frequency", "07_acoustic_helmholtz", "10_ngsolve_bem_reference"])
    reference = "COMSOL acoustic FEM/BEM internal";
else
    reference = "";
end
end


function source = gypsilabInspirationFor(category, title)
switch category
    case "01_mesh_topology"
        source = "openMsh + nonRegressionTest/meshManagement";
    case "02_h1_scalar_fem"
        source = "openFem + nonRegressionTest/finiteElement";
    case "03_hcurl_edge_fem"
        source = "openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m";
    case "04_laplace_bem_dense"
        source = "openOpr + nonRegressionTest/operators + openEbd scalar products";
    case "05_laplace_hmatrix"
        source = "openHmx + nonRegressionTest/hierarchicalMatrix";
    case "06_acoustic_low_frequency"
        source = "openEbd Helmholtz kernels + radiationImpedances + acoustic papers";
    case "07_acoustic_helmholtz"
        source = "miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d";
    case "08_scalar_fem_bem_coupling"
        source = "doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic";
    case "09_rwg_hcurl_trace"
        source = "openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique";
    case "10_ngsolve_bem_reference"
        source = "Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM";
    otherwise
        source = "Gypsilab readable FEM/BEM examples";
end
source = source + " / " + title;
end


function cases = markMeshTopologyVerified(cases)
logPath = "S:\MATLAB\_crossval\gypsilab_mesh_topology_10of100_20260624.md";
for k = 1:10
    cases(k).status = "verified";
    cases(k).validationLog = logPath;
end
end


function text = slug(title)
text = lower(regexprep(title, "[^a-zA-Z0-9]+", "_"));
text = regexprep(text, "^_+|_+$", "");
end
