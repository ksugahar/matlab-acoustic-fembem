function tests = testDrumScatterGmsh3dArtifact
%TESTDRUMSCATTERGMSH3DARTIFACT 3D Gmsh scene of the drum radiating into a scatterer.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumScatterSceneWritesTwoBodyScatteringField(testCase)
outBase = string(tempname(tempdir));
testCase.addTeardown(@() cleanupArtifacts(outBase));

artifact = writeDrumScatterGmsh3dArtifact("OutputBase", outBase, ...
    "NumBeats", 2, "TimeStep", 0.25, "NumGrid", 36);

verifyEqual(testCase, artifact.schema, "matlab-acoustic-fembem.drum-scatter-gmsh3d.v1");
verifyEqual(testCase, artifact.status, "ok");
verifyEqual(testCase, artifact.drum_shape, "cylinder");
verifyEqual(testCase, artifact.scatterer_shape, "sphere");
verifyEqual(testCase, cellstr(artifact.beat_spots(:)).', {'A', 'B'});

% two disjoint bodies: a cylinder drum struck on top + a sphere scatterer above
verifyTrue(testCase, artifact.checks.two_bodies_detected);
verifyEqual(testCase, artifact.num_bodies, 2);
verifyTrue(testCase, artifact.checks.cylinder_not_sphere);
verifyTrue(testCase, artifact.checks.scatterer_present);
verifyGreaterThan(testCase, artifact.scatterer_nodes, 0);

% the scatterer floats above the drumhead, directly over the rim
verifyTrue(testCase, artifact.checks.scatterer_above_drum);
verifyGreaterThan(testCase, artifact.scatterer_bottom_z, artifact.drum_top_z);
verifyTrue(testCase, artifact.checks.scatterer_over_drum_rim);
verifyLessThan(testCase, abs(artifact.scatterer_offaxis_radius - artifact.drum_radius), ...
    0.35 * artifact.drum_radius);

% the drum is struck at +x / -x on the top face
verifyGreaterThan(testCase, artifact.strike_spot_a(1), 0);
verifyLessThan(testCase, artifact.strike_spot_b(1), 0);
verifyEqual(testCase, artifact.strike_spot_a(3), artifact.strike_spot_b(3));

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
