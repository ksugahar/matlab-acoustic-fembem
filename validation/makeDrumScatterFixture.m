function makeDrumScatterFixture(outFile)
%MAKEDRUMSCATTERFIXTURE Generate the drum + scatterer mesh fixture (drum_scatter.vol).
%
%   makeDrumScatterFixture() writes
%   fixtures/mesh_topology/drum_scatter.vol: a CYLINDER DRUM (radius 1.0,
%   height 0.6, so the drumhead is the top face at z = 0.6 and the bottom
%   face at z = 0) UNION a SPHERE SCATTERER (radius 0.5) placed directly
%   above the drum RIM at center [1.0 0 1.4] -- its off-axis radius (1.0)
%   equals the drum radius, and its bottom (z = 0.9) floats 0.3 above the
%   drumhead.  The two bodies are disjoint, so the .vol carries two closed
%   boundary surfaces (the same two-body class drumScatterField
%   and the CQ single-layer solver expect).
%
%   The mesh is generated with the PDE Toolbox (multicylinder + multisphere
%   + union + generateMesh, linear order, Hmax = 0.28) and written via
%   writePdeMeshVol.  Runtime code only READS the committed .vol -- this
%   generator exists so the fixture is reproducible, not so it runs live.
%
%   makeDrumScatterFixture(outFile) writes to a custom path.

arguments
    outFile (1,1) string = fullfile(gypsilabRepoRoot(), ...
        "fixtures", "mesh_topology", "drum_scatter.vol")
end

drumRadius = 1.0; drumHeight = 0.6;
scRadius = 0.5;   scCenter = [1.0 0 1.4];
Hmax = 0.28;

g = union(multicylinder(drumRadius, drumHeight), ...
          translate(multisphere(scRadius), scCenter));
gm = generateMesh(g, "Hmax", Hmax, "GeometricOrder", "linear");
writePdeMeshVol(gm.Mesh, outFile, MaterialName="domain", BoundaryName="outer");

% ---------- report the two-body split so the fixture can be sanity-checked ----------
surface = VolMesh(outFile).boundary();
X = surface.vtx; tri = surface.tri; nB = size(X, 1);
E = unique(sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])], 2), "rows");
comp = conncomp(graph(E(:, 1), E(:, 2), [], nB)).';
loZ = accumarray(comp, X(:, 3), [max(comp) 1], @min);
[~, drumBody] = min(loZ);
nDrum = nnz(comp == drumBody); nScat = nB - nDrum;
zTop = max(X(comp == drumBody, 3));
fprintf("wrote %s\n", outFile);
fprintf("  surface nodes = %d, bodies = %d (drum = %d, scatterer = %d)\n", ...
    nB, max(comp), nDrum, nScat);
fprintf("  drum top z = %.3f, scatterer center = [%.2f %.2f %.2f], radius ~ %.3f\n", ...
    zTop, scCenter, scRadius);
end
