function tests = testDrumRollGmsh3dArtifact
%TESTDRUMROLLGMSH3DARTIFACT 3D Gmsh scene of the cylinder drum roll.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumRollGmsh3dSceneWritesCylinderDrumAndField(testCase)
outBase = string(tempname(tempdir));
testCase.addTeardown(@() cleanupArtifacts(outBase));

artifact = writeDrumRollGmsh3dArtifact("OutputBase", outBase, ...
    "NumBeats", 2, "TimeStep", 0.25, "NumGrid", 36);

verifyEqual(testCase, artifact.schema, "matlab-acoustic-fembem.drum-roll-gmsh3d.v1");
verifyEqual(testCase, artifact.status, "ok");
verifyEqual(testCase, artifact.drum_shape, "cylinder");
verifyEqual(testCase, cellstr(artifact.beat_spots(:)).', {'A', 'B'});

% the scene is a CYLINDER drum struck on its top drumhead at two spots
verifyTrue(testCase, artifact.checks.cylinder_not_sphere);
verifyTrue(testCase, artifact.checks.drumhead_two_spots_top_face);
verifyEqual(testCase, artifact.strike_spot_a(3), artifact.strike_spot_b(3));  % same top-face z
verifyGreaterThan(testCase, artifact.strike_spot_a(1), 0);
verifyLessThan(testCase, artifact.strike_spot_b(1), 0);

% both Gmsh views written: x-z pressure scalar + deforming drumhead vector
verifyTrue(testCase, artifact.checks.pressure_xz_plane_written);
verifyTrue(testCase, artifact.checks.deforming_drum_written);
verifyTrue(testCase, artifact.checks.alternating_two_spot);
verifyTrue(testCase, artifact.checks.cq_residual_small);
verifyTrue(testCase, artifact.checks.gmsh_v41);

verifyTrue(testCase, isfile(artifact.gmsh_msh));
verifyTrue(testCase, isfile(artifact.gmsh_geo));
verifyTrue(testCase, isfile(artifact.gmsh_geo_opt));
verifyTrue(testCase, contains(firstLines(artifact.gmsh_msh, 2), "4.1"));  % v4.1 header
end


function head = firstLines(path, n)
fid = fopen(path, "r");
cleanup = onCleanup(@() fclose(fid));
lines = strings(1, 0);
for i = 1:n
    lines(end + 1) = string(fgetl(fid)); %#ok<AGROW>
end
head = join(lines, " ");
end


function cleanupArtifacts(outBase)
for ext = [".msh", ".geo", ".geo.opt", ".opt", ".result.json"]
    p = char(outBase + ext);
    if isfile(p)
        delete(p);
    end
end
end
