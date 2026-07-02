function tests = testLaplaceDirichletSolve
%TESTLAPLACEDIRICHLETSOLVE Interior Dirichlet BVP rung of the ladder.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testLinearPatchTestIsExact(testCase)
% P1 reproduces a linear potential exactly: the classic patch test.
model = FemBemModel(fixturePath("four_tet_interior_node.vol"));
verifyNotEmpty(testCase, setdiff((1:size(model.mesh.vtx, 1)).', ...
    model.mesh.traceNodeIds), "fixture must have an interior node");

a0 = 1.0; a = [2.0; 3.0; 4.0];
uExact = a0 + model.mesh.vtx * a;
g = uExact(model.mesh.traceNodeIds);

sol = laplaceDirichletSolve(model, g);

verifyEqual(testCase, sol.status, "ok");
verifyEqual(testCase, sol.kind, "laplace_dirichlet_interior_solve");
verifyEqual(testCase, sol.u, uExact, "AbsTol", 1e-12);
verifyLessThan(testCase, sol.interiorResidualNorm, 1e-12);
verifyTrue(testCase, sol.checks.boundaryValuesImposedExactly);
verifyTrue(testCase, sol.checks.interiorEquationsSatisfied);
verifyTrue(testCase, sol.checks.reactionBalancesToZero);
end


function testUniformCoefficientDoesNotChangeThePatchTest(testCase)
% -div(c grad u) = 0 with constant c keeps linear fields exact.
model = FemBemModel(fixturePath("four_tet_interior_node.vol"));
uExact = -0.5 + model.mesh.vtx * [1.0; -2.0; 0.5];
g = uExact(model.mesh.traceNodeIds);

sol = laplaceDirichletSolve(model, g, "MaterialCoef", 2.5);

verifyEqual(testCase, sol.status, "ok");
verifyEqual(testCase, sol.u, uExact, "AbsTol", 1e-12);
end


function testEnergyMatchesQuadraticForm(testCase)
model = FemBemModel(fixturePath("four_tet_interior_node.vol"));
g = sin((1:numel(model.mesh.traceNodeIds)).');

sol = laplaceDirichletSolve(model, g);

[K, ~] = model.h1.stiffness();
verifyEqual(testCase, sol.energy, 0.5 * (sol.u.' * K * sol.u), "AbsTol", 1e-12);
verifyGreaterThanOrEqual(testCase, sol.energy, 0);
end


function testBoundaryOnlyMeshSkipsTheInteriorSolve(testCase)
% The unit tetra has no interior node: the solution is the boundary data.
model = FemBemModel(fixturePath("unit_tetra.vol"));
g = [1; 2; 3; 4];

sol = laplaceDirichletSolve(model, g);

verifyEqual(testCase, sol.solver, "boundary_only_no_interior_unknowns");
verifyEqual(testCase, sol.u, g);
verifyEmpty(testCase, sol.interiorNodeIds);
verifyEqual(testCase, sol.status, "ok");
end


function testPerTetCoefficientVectorIsAccepted(testCase)
model = FemBemModel(fixturePath("two_material_tetra_labels.vol"));
coef = 1 + double(model.mesh.tetMat(:));
g = ones(numel(model.mesh.traceNodeIds), 1);

sol = laplaceDirichletSolve(model, g, "MaterialCoef", coef);

verifyEqual(testCase, sol.status, "ok");
verifyEqual(testCase, sol.materialCoef, coef);
% constant boundary data keeps the constant field for any coefficient
verifyEqual(testCase, sol.u, ones(size(model.mesh.vtx, 1), 1), "AbsTol", 1e-12);
end


function testRejectsWrongBoundaryLength(testCase)
model = FemBemModel(fixturePath("unit_tetra.vol"));
verifyError(testCase, @() laplaceDirichletSolve(model, [1; 2]), ...
    "laplaceDirichletSolve:boundary");
end


function path = fixturePath(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
path = string(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end
