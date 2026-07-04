function text = fembem_knowledge(topic)
%FEMBEM_KNOWLEDGE Queryable FEM/BEM knowledge for the Sugahara-lab MATLAB MCP.
%
%   text = acoustic_fembem.fembem_knowledge();            % overview (default)
%   text = acoustic_fembem.fembem_knowledge("acoustic");  % one topic
%
% Sugahara-lab is a FEM/BEM lab; this makes the lab MATLAB MCP server
% *know* FEM/BEM the way it already knows optimization (the Touchstone /
% Tikhonov / L-curve / Morozov lessons in CONVERSION_LESSONS.md). The
% underlying implementation is the integrated Gypsilab readable lane
% (acoustic_fembem.repository_root), exercised by the runnable acoustic_fembem.fembem_acoustic_gate.
%
% Topics: overview, spaces, galerkin_bem, coupled_fem_bem, multiphysics,
%         acoustic, sonic_crystal, adjoint_ad, matlab_execution_policy,
%         vol_visualization, pde_vol_bridge,
%         radia_ngsolve_crossval, validation_discipline, optimization_link,
%         all.

arguments
    topic (1,1) string = "overview"
end

t = lower(strip(topic));
switch t
    case {"overview", "intro", ""}
        text = OVERVIEW;
    case {"adjoint_ad", "adjoint", "ad", "automatic_differentiation", ...
            "sensitivity", "wavefront", "inverse_design"}
        text = ADJOINT_AD;
    case {"spaces", "first_order_spaces", "h1", "hcurl", "rwg"}
        text = SPACES;
    case {"galerkin_bem", "bem", "single_layer", "double_layer"}
        text = GALERKIN_BEM;
    case {"coupled_fem_bem", "coupled", "johnson_nedelec", "transmission"}
        text = COUPLED_FEM_BEM;
    case {"acoustic", "helmholtz", "scattering"}
        text = ACOUSTIC;
    case {"sonic_crystal", "band_gap", "bloch", "metamaterial", "duct"}
        text = SONIC_CRYSTAL;
    case {"validation_discipline", "validation", "gates", "cross_check"}
        text = VALIDATION_DISCIPLINE;
    case {"radia_ngsolve_crossval", "radia-ngsolve", "ngsolve_crossval", ...
            "vol_crossval", "vol"}
        text = RADIA_NGSOLVE_CROSSVAL;
    case {"pde_vol_bridge", "pde_toolbox", "generate_mesh", "matlab_mesh"}
        text = PDE_VOL_BRIDGE;
    case {"vol_visualization", "visualization", "netgen_viewer", "vol_viewer"}
        text = VOL_VISUALIZATION;
    case {"matlab_execution_policy", "execution_policy", "no_live_documents", ...
            "scripts", "mcp_json"}
        text = MATLAB_EXECUTION_POLICY;
    case {"optimization_link", "optimization", "inverse", "design"}
        text = OPTIMIZATION_LINK;
    case {"multiphysics", "interface", "coupling_difficulty", ...
            "nonconforming", "mortar", "mesh_coupling"}
        text = MULTIPHYSICS;
    case "all"
        text = strjoin([OVERVIEW, SPACES, GALERKIN_BEM, COUPLED_FEM_BEM, ...
            MULTIPHYSICS, ACOUSTIC, SONIC_CRYSTAL, ADJOINT_AD, ...
            MATLAB_EXECUTION_POLICY, VOL_VISUALIZATION, PDE_VOL_BRIDGE, ...
            RADIA_NGSOLVE_CROSSVAL, VALIDATION_DISCIPLINE, ...
            OPTIMIZATION_LINK], [newline newline]);
    otherwise
        text = "Unknown topic '" + topic + "'. Available: overview, " + ...
            "spaces, galerkin_bem, coupled_fem_bem, multiphysics, acoustic, " + ...
            "sonic_crystal, adjoint_ad, matlab_execution_policy, " + ...
            "vol_visualization, pde_vol_bridge, " + ...
            "radia_ngsolve_crossval, validation_discipline, optimization_link, all.";
end
end


function s = OVERVIEW()
s = strjoin([
    "# Sugahara-lab FEM/BEM (integrated Gypsilab readable lane)"
    ""
    "One user-facing path as short as Gypsilab, with the mathematics"
    "visible in each class (see the Gypsilab READABLE_CLASS_STYLE)."
    "The lab MATLAB MCP drives it through two runnable gates"
    "(acoustic_fembem.fembem_acoustic_gate) and this knowledge dispatcher, so"
    "FEM/BEM is a first-class lab capability alongside optimization."
    ""
    "Canonical mesh intake: Netgen .vol (surface triangles + volume"
    "tetrahedra only; curved/quad/hex/wedge rejected fail-loud)."
    ""
    "Cross-validation ladder (landed and tested):"
    "  3  interior Dirichlet FEM (P1 patch test 1e-12)"
    "  4  exterior Galerkin single layer (sphere capacitance; 3-way vs"
    "     the real Gypsilab AND ngsolve.bem)"
    "  5  Johnson-Nedelec coupled FEM/BEM (analytic unit ball)"
    "  6  H(curl)/RWG vector trace coupling (magnetized-sphere gate)"
    "  7  acoustic Helmholtz single layer (3-way validated)"
    "  8  sonic-crystal chain (4-leg; free space has NO Bragg gap)"
    "  9  duct band gap (Bloch cell + finite-crystal transmission)"
    " 10  coupled acoustic transmission (invisibility + Anderson sphere)"
    " 11  rigid scattering + irregular frequencies + CHIEF"
    ""
    "Topics: spaces, galerkin_bem, coupled_fem_bem, multiphysics,"
    "acoustic, sonic_crystal, matlab_execution_policy, vol_visualization,"
    "pde_vol_bridge, radia_ngsolve_crossval,"
    "validation_discipline, optimization_link."
    ], newline);
end


function s = SPACES()
s = strjoin([
    "# First-order spaces (one class per mathematical object)"
    ""
    "VolMesh / SurfaceMesh  - parsed .vol + compact boundary, orientation"
    "                         signs to outward (fail-loud if unknown)."
    "H1Space                - P1 tetrahedra: stiffness(c), mass(c)."
    "Nedelec0Space          - HCurl volume edges: matrices/mass/curlCurl."
    "SurfaceP1Space         - boundary P1 scalar trace + exact mass."
    "RwgSpace               - boundary RWG0; rotatedTraceMap gives the RWG"
    "                         coefficients of n x u|_Gamma EXACTLY from a"
    "                         volume Nedelec field (FEEC trace identity,"
    "                         machine precision)."
    "TraceOperator          - one-hot H1 volume -> boundary P1."
    ""
    "Lesson (locked): a boundary-adjacency identity RECORDED in a package"
    "is not VERIFIED until a consumer uses it - a swapped ismember in the"
    ".vol reader silently returned tet 1 for every boundary triangle and"
    "passed every convex-mesh gate until the RWG rung consumed the value."
    ], newline);
end


function s = GALERKIN_BEM()
s = strjoin([
    "# Galerkin BEM (Laplace + Helmholtz, the Gypsilab split)"
    ""
    "GalerkinSingleLayer / GalerkinDoubleLayer, kernel"
    "  G_k(r) = 1/(4 pi r) + (exp(i k r) - 1)/(4 pi r)"
    "split like Gypsilab integral + regularize: the singular Laplace part"
    "is integrated ANALYTICALLY per source triangle (Wilton closed forms,"
    "machine precision vs subdivision/polar), the smooth"
    "low-frequency-stable correction by quadrature. k -> 0 is exact by"
    "construction. Double layer = outward-normal principal value; the"
    "+-1/2 jump lives in the BIE, not the matrix."
    ""
    "Sphere spectral gates pin the conventions NUMERICALLY:"
    "  V[Y_l] = i k j_l(k) h_l(k)"
    "  K[Y_l] = 1/2 + i k^2 j_l(k) h_l'(k)   (Laplace: -1/(2(2l+1)))"
    "The imaginary signs distinguish e^{+ikr} from e^{-ikr}; an operator"
    "must match a trusted reference far better than its conjugate."
    ""
    "Scaling: the all-dense correction is nGauss x nGauss - gss 7 on a"
    "634-triangle multi-body mesh is a ~315 MB transient; the H-matrix"
    "teaching path (HMatrix) is where compression belongs."
    ], newline);
end


function s = COUPLED_FEM_BEM()
s = strjoin([
    "# Johnson-Nedelec coupled FEM/BEM (open boundary, NO absorbing BC)"
    ""
    "The BIE row IS the exact radiation condition - the coupled solve"
    "needs no ABC (a volume-FEM leg with a first-order Sommerfeld ABC"
    "sits ~4-8% off at exterior probes; the coupling has zero such error)."
    ""
    "Laplace:  A u - T' M lambda = F ;  (1/2 M - K) T u + V lambda = 0."
    "Helmholtz (acoustic transmission, medium 1 inside / medium 0 out):"
    "  (1/rho1)(A - k1^2 Mv) u - T' M sigma = F + T' G_inc"
    "  (1/2 M - K_k) T u + V_k sigma = (1/2 M - K_k) g_inc"
    "exterior scattered u_s = -S_k[sigma] + D_k[u_s|_Gamma]."
    ""
    "Gates: unit ball f=1 -> u = 1/2 - r^2/6 (Laplace); acoustic"
    "INVISIBILITY (k1=k0, rho1=rho0) is the patch test - the scatterer"
    "must vanish, the residual is the discretization null and must"
    "CONVERGE (4.1e-2 -> 1.3e-2, fine -> finer ball); real contrast is"
    "gated by the Anderson (1950) fluid-sphere series."
    ""
    "WHERE COUPLING TRULY EARNS ITS KEEP: a non-trivial bounded interior"
    "(inhomogeneous / nonlinear / ELASTIC) with an unbounded homogeneous"
    "exterior. A rigid/soft scatterer needs NO interior PDE (pure BEM is"
    "right); the coupling pays off when the interior physics is real -"
    "the FSI (fluid-structure) case: an ELASTIC bead has internal"
    "compressional + shear resonances that rigid/soft MISS. Analytic"
    "reference landed: elasticSphereScattering (Faran 1951), validated"
    "vs Anderson (shear->0, 1.3e-10) AND rigid (stiff limit); its"
    "radiation force Y_p(kR) peaks at a genuine internal resonance"
    "(lucite-like: 2.83 at kR=3 vs the rigid ~0.75). The FSI COUPLED"
    "SOLVE LANDED (fsiCoupledSolve): vector P1 elasticity FEM interior +"
    "acoustic BEM + interface (pressure<->normal traction, normal"
    "velocity<->normal displacement); unknowns (u, p_s, q_s), the scalar"
    "Johnson-Nedelec pattern with a VECTOR interior. Validated vs the"
    "elastic sphere: stiff limit -> rigid sphere 5e-4 (formulation gate,"
    "coupling EXACT), elastic field CONVERGES under refinement (25% coarse"
    "-> 7% fine; residual is the P1 INTERIOR elasticity, not the"
    "coupling). Run it: acoustic_fembem.fembem_acoustic_gate(""fsi"")."
    ""
    "FAST EXTERIOR (Kelvin where Kelvin applies): for a SPHERE truncation,"
    "ExteriorMethod=""dtn"" swaps the dense Galerkin single/double layer for"
    "the EXACT spherical Helmholtz DtN (sphericalDtnOperator, Lambda_n ="
    "k h_n'(kR)/h_n(kR) - the Kelvin operator on the sphere / its radiating"
    "extension, Sugahara IEICE 2024). The scattered field is a spherical-"
    "harmonic expansion p_s = Phi c, so the exterior reduces to (N+1)^2"
    "coefficients and the coupled system stays FULL RANK (a nodal low-rank"
    "DtN leaves p_s under-determined - the bug that made the naive nodal swap"
    "574% wrong; the harmonic-coefficient reduction is the fix). DtN operator"
    "exact per multipole (independent point-source D->N check 2.6e-5 at"
    "degree 10), coupled field matches the rigid sphere to 3.8e-3 (BEM leg"
    "4.2e-3; DtN-vs-BEM 7e-4), singular N^2 assembly (~85 s) skipped ->"
    "sub-second solve. FAIL-LOUD on a non-sphere. The acoustic instance of"
    "the lab exact-operator (NOT absorber) open-boundary paradigm. Run it:"
    "acoustic_fembem.fembem_acoustic_gate(""fsi_dtn"")."
    ""
    "This SCALAR solve is the precursor. Lab EM analogues where coupling is the"
    "workhorse: nonlinear iron + open boundary, eddy currents + open"
    "boundary - the Kelvin transform / Radia integral method are the same"
    "problem class."
    ], newline);
end


function s = ACOUSTIC()
s = strjoin([
    "# Acoustic scattering (Helmholtz, e^{+ikr})"
    ""
    "Exterior Dirichlet (sound-soft): first-kind V_k q = M g. Rigid"
    "(sound-hard): total-field second-kind (1/2 M - K_k) t = M g_inc,"
    "u = u_inc + D_k[t] (single-layer term drops since dp/dn = 0)."
    ""
    "IRREGULAR FREQUENCIES: both fail at the interior Dirichlet"
    "eigenvalues of the surface (unit sphere: first at kR = pi). The"
    "sharpest lesson in the lab suite: at kR = pi the solution is ~96%"
    "wrong while the DISCRETE condition number stays benign (~29, the"
    "faceted eigenvalue shifts off pi) - condition monitoring cannot"
    "detect it, only analytic gates do. CHIEF (Schenck 1968; interior"
    "null-field rows at jittered points, least squares) restores"
    "regular-class accuracy with zero cost at regular k; Burton-Miller"
    "(1971; combined field, hypersingular) is the production alternative."
    ""
    "Analytic references: soft/rigid single-sphere partial-wave series,"
    "the interior point source (EXACT exterior-reproduction gate), Foldy"
    "monopole multiple scattering (LOW-k, with a MEASURED validity"
    "envelope), Anderson fluid-sphere transmission (log-derivative-stable"
    "per-mode solve; size the series on the FARTHEST probe k*r_max)."
    "Run any of these via acoustic_fembem.fembem_acoustic_gate."
    ], newline);
end


function s = SONIC_CRYSTAL()
s = strjoin([
    "# Sonic crystals / acoustic band gaps: confinement is not optional"
    ""
    "A sparse FREE-SPACE chain of strong scatterers develops NO Bragg"
    "stop band at k d = pi (verified four ways: Foldy, two BEM codes, a"
    "volume FEM). Insertion loss is a flat broadband plateau - a"
    "sound-soft sphere keeps scattering length ~ R down to k -> 0, and"
    "in open 3D the energy leaks sideways so the 1D interference never"
    "accumulates."
    ""
    "CONFINE the same chain in a rigid duct (single-mode below the"
    "transverse cutoff) and the gap appears: Bloch band structure of the"
    "unit cell (Floquet phase; empty-lattice analytic gate) + finite"
    "N-cell transmission (one-way ports; empty duct transparent) show"
    "stop / pass / stop aligned quantitatively (contrast > 100)."
    ""
    "Bug class the gates catch: OCC face classification by AREA"
    "thresholds silently turns the cell into a Dirichlet cavity and"
    "misses the inclusions - classify GEOMETRICALLY (face center on a"
    "wall plane). The empty-lattice gate also caught a wrong analytic"
    "reference table (forgotten transverse mode families) - gates cut"
    "both ways. Public write-up: radia_mcp.metamaterial"
    "acoustic_sonic_crystal."
    ], newline);
end


function s = ADJOINT_AD()
s = strjoin([
    "# Adjoint automatic differentiation through the BEM/FEM solve"
    ""
    "Yes, AD works here - the readable, correct way: differentiate the"
    "EQUATION (implicit-function / adjoint rule), never the LU algorithm."
    "MATLAB's black-box AD (dlarray) is a dead end for this stack (complex"
    "+ sparse + a linear solve are unsupported/impractical); complex-step"
    "is blocked because the imaginary axis is occupied by e^{ikr}."
    ""
    "For a linear system A(q) u = b(q) and objective J(u,q):"
    "  A' lambda = (dJ/du)'   (one transpose solve)"
    "  dJ/dq = dJ/dq|_explicit - lambda' (dA/dq u - db/dq)"
    "Reverse mode: ONE extra solve for the whole gradient, ANY number of"
    "design variables (vs one forward re-solve per variable for finite"
    "differences). AD helps only in the small assembly derivatives"
    "dA/dq, db/dq; the solve is handled by the adjoint identity."
    ""
    "Landed (Gypsilab acousticFocusAdjoint): phased-array WAVEFRONT"
    "SYNTHESIS - design the complex source amplitudes p to focus"
    "J = |u(target)|^2 behind a rigid scatterer. The field is affine in p"
    "(u = w p, w = S0 + lambda' M S), so w is EXACT (residual 4e-18) and"
    "the gradient matches central finite differences to 1.7e-10; a"
    "gradient ascent focuses the energy. Run it: acoustic_fembem.fembem_acoustic_gate"
    "(""focus_adjoint"") returns the adjoint-vs-FD relative error."
    ""
    "Difficulty map (what is easy vs hard):"
    "  parameter derivatives (k, material/density contrast, source"
    "    amplitude/phase): EASY - smooth, non-singular dependence;"
    "  shape derivatives (moving/resizing a scatterer): HARD - singular"
    "    BEM shape sensitivity;"
    "  non-holomorphic objectives (|u|^2, radiation force): use Wirtinger"
    "    calculus. Sign trap (recorded): the steepest-ASCENT direction is"
    "    dJ/dconj(p) = 2 u conj(w), NOT 2 conj(u) w - a near-zero field"
    "    hides the flip, so test a gradient's SIGN away from the minimum."
    ""
    "Radiation-force / thrust (LANDED - acousticRadiationForce): the net"
    "time-averaged force is the control-sphere integral of the Brillouin"
    "radiation-stress tensor of the LINEAR field, a quadratic form in the"
    "surface trace, so the SAME one-adjoint-solve reverse mode gives"
    "dF/d(phases) for wavefront-synthesised thrust. Run it:"
    "acoustic_fembem.fembem_acoustic_gate(""radiation_force""). Air @ 40 kHz"
    "(lambda 8.58 mm): a kR=2 (2.73 mm) sphere feels ~5 uN at 140 dB,"
    "~0.5 mN at 160 dB (levitation regime; F ~ W/c ~ 2.9 mN per radiated"
    "watt in air)."
    ""
    "CAPSTONE (LANDED - elasticThrustAdjoint): the SAME reverse-mode idea one"
    "level up - steer the radiation FORCE on a solid ELASTIC bead (through the"
    "FSI coupled solve) by the array phases. The field is LINEAR in the"
    "amplitudes (one FSI solve per source, potentials shared over sources), so"
    "the force is an exact QUADRATIC FORM F_i(p) = p^H Q_i p (Q_i Hermitian,"
    "from the Brillouin stress); the Wirtinger ascent is dF_i/dconj(p) = 2 Q_i"
    "p, so force AND full gradient are closed forms once Q_i is built - a design"
    "loop REUSES forceForm with NO re-solve. Validated: the form reproduces a"
    "vectorised direct Brillouin integral to ~1e-15 (INDEPENDENT assembly, which"
    "caught a conj(z) z.' vs z z' outer-product conjugation error that only"
    "bites for COMPLEX amplitudes), matches the golden acousticRadiationForce,"
    "is control-radius independent, and the Wirtinger gradient matches FD to"
    "~1e-13. Run it: acoustic_fembem.fembem_acoustic_gate(""thrust_adjoint"")."
    ""
    "This is why AD lives in a FEM/BEM lab's MCP: the"
    "solver and the design loop that drives it, together."
    ], newline);
end


function s = PDE_VOL_BRIDGE()
s = strjoin([
    "# PDE Toolbox to .vol bridge"
    ""
    "Goal: simple geometries should not require an external mesher."
    "MATLAB PDE Toolbox can generate a linear tetrahedral mesh, then"
    "writePdeMeshVol turns that mesh into the same Netgen .vol contract"
    "used by VolMesh, FemBemModel, and the BEM boundary views."
    ""
    "Minimal path:"
    "  model = createpde();"
    "  model.Geometry = multicuboid(1, 1, 1);"
    "  generateMesh(model, ""GeometricOrder"", ""linear"");"
    "  writePdeMeshVol(model.Mesh, ""box.vol"");"
    "  mesh = VolMesh(""box.vol"");"
    ""
    "Policy: no implicit element conversion. The bridge accepts only"
    "linear tetrahedral PDE meshes; boundary faces are derived as"
    "outward-oriented triangles, and the round trip is checked by"
    "readVolTriTet. For very complex CAD, Cubit/Netgen remain better,"
    "but for teaching cubes, boxes, balls, and parameter sweeps, MATLAB"
    "alone can now create the solver-facing .vol input."
    ], newline);
end


function s = VOL_VISUALIZATION()
s = strjoin([
    "# .vol visualization policy"
    ""
    "Best native GUI viewer: Netgen. The .vol file is Netgen's own mesh"
    "handoff format, so Netgen is the most faithful place to inspect"
    "surface triangles, volume tetrahedra, material labels, and boundary"
    "names interactively."
    ""
    "MATLAB role: quick figure preview and solver preflight."
    "  plotVolMesh(""mesh.vol"")                 -> boundary triangle preview"
    "  acoustic_fembem.vol_mesh_summary(""mesh.vol"") -> counts, bbox, labels"
    ""
    "Practical split:"
    "  - GUI/user inspection: Netgen native .vol viewer."
    "  - MATLAB script/report: plotVolMesh plus the summary JSON."
    "  - LLM/headless preflight: acoustic_fembem_vol_mesh_summary."
    ""
    "Do not turn visualization into a hidden mesh conversion step. The"
    "solver-facing contract remains first-order tri/tet .vol, and the"
    "reader rejects quad/hex/wedge/pyramid/curved records fail-loud."
    ], newline);
end


function s = MATLAB_EXECUTION_POLICY()
s = strjoin([
    "# MATLAB execution policy"
    ""
    "Do not make MATLAB Live Editor documents the default product surface"
    "for this repository.  The durable interfaces are normal MATLAB functions,"
    "small .m scripts, JSON result manifests, and MCP custom tools."
    ""
    "Preferred split:"
    "  - human mesh inspection: Netgen native .vol viewer"
    "  - MATLAB local work: .m functions/scripts plus figures from plotVolMesh"
    "  - LLM/headless work: MCP tools that print compact JSON"
    "  - reproducible artifacts: JSON manifests with versions, run dates,"
    "    timing breakdowns, schema ids, and convention ids"
    ""
    "This keeps the solver readable for students while avoiding another UI"
    "format to maintain.  If a report is needed, generate it from scripts or"
    "ordinary documentation after the manifest and gates have passed."
    ], newline);
end


function s = VALIDATION_DISCIPLINE()
s = strjoin([
    "# BEM/FEM validation discipline (the lab's standing rules)"
    ""
    "1. Conventions are MEASURED, never assumed: sphere eigenvalues pin"
    "   V_k, K_k and the e^{+ikr} sign; an operator must beat its"
    "   conjugate by >~40x. The armchair Green-identity derivation had a"
    "   sign slip (BIE constant 3/2); the Gauss check + two spherical"
    "   modes fixed it."
    "2. Prefer an EXACT gate: an interior point source must be reproduced"
    "   at every exterior point (uniqueness + radiation), zero truncation"
    "   - it generalizes verbatim to multi-body."
    "3. MULTI-code, MULTI-method: the two BEM codes agree 10-30x tighter"
    "   with each other than either matches the true sphere (=> the"
    "   deviation is shared faceting, not a bug); add a volume-FEM+ABC"
    "   leg so a shared BEM formulation blind spot cannot hide."
    "4. Approximate references get a MEASURED envelope (Foldy: ~4% at"
    "   kR=0.18 -> ~47% at kR=0.9) and lock the degradation itself."
    "5. Two-fixture CONVERGENCE assertions turn P1 (k h)^2 resolution"
    "   limits into locked tests, not loose bands."
    "6. Committed reference artifacts store their own intorder"
    "   self-convergence so a stale .mat fails loudly."
    ""
    "Runnable now: acoustic_fembem.fembem_acoustic_gate(kind) returns a JSON"
    "verdict (relative error vs the analytic series, pass within a"
    "per-kind band)."
    "Cross-code gate: acoustic_fembem.fembem_crossval_gate(kind) returns a JSON"
    "verdict against the radia-ngsolve/NGSolve validation ladder using"
    "the same Netgen .vol tri/tet fixtures."
    ], newline);
end


function s = RADIA_NGSOLVE_CROSSVAL()
s = strjoin([
    "# radia-ngsolve / NGSolve cross validation"
    ""
    "Both sides use the SAME Netgen .vol input. MATLAB/Gypsilab reads the"
    "volume tetrahedra as the FEM view and the boundary triangles as the"
    "BEM view; radia-ngsolve/NGSolve reads the same .vol fixture as the"
    "reference-side mesh. This makes node order, boundary triangles, and"
    "tetrahedra part of the validation contract, not incidental setup."
    ""
    "Accepted topology in this teaching lane:"
    "  - surface triangles"
    "  - volume tetrahedra"
    "  - first-order only"
    "Quad, hex, wedge, pyramid, curved, or hidden splitting are rejected"
    "on the MATLAB side. NGSolve can read broader .vol families, but this"
    "cross-validation lane deliberately uses tri/tet fixtures."
    ""
    "Runnable MCP gate:"
    "  acoustic_fembem.fembem_crossval_gate(""mesh_topology"")      -> live .vol intake"
    "  acoustic_fembem.fembem_crossval_gate(""galerkin_ngsolve"")   -> Laplace P1 BEM"
    "  acoustic_fembem.fembem_crossval_gate(""helmholtz_ngsolve"")  -> Helmholtz BEM"
    "  acoustic_fembem.fembem_crossval_gate(""catalog_100"")        -> 100-case catalog"
    ""
    "Interpretation discipline:"
    "  mesh_topology checks file-format identity; Galerkin checks BEM"
    "  operator conventions and errors; Helmholtz checks e^{+ikr} by"
    "  requiring the operator to match the reference much better than its"
    "  conjugate. A passing artifact can be promoted into MCP knowledge"
    "  only after the gate, artifact metadata, and focused verifier pass."
    ], newline);
end


function s = OPTIMIZATION_LINK()
s = strjoin([
    "# FEM/BEM <-> optimization (why both live in this lab MCP)"
    ""
    "The optimization knowledge already in the lab MCP"
    "(CONVERSION_LESSONS.md: Touchstone port/match/preflight, Tikhonov"
    "path, L-curve corner, Morozov discrepancy, box-constrained and"
    "trace least squares) is the DESIGN side of the same FEM/BEM story:"
    ""
    "- FEM/BEM trace least squares fits a boundary trace to a target"
    "  with T, boundary mass M, Tikhonov weight, and a finite-difference"
    "  gradient check all visible (femBemTraceLeastSquares)."
    "- inverse acoustic/EM design (sonic-crystal band-gap tuning, coil"
    "  shaping) is a forward FEM/BEM solve wrapped in a regularized"
    "  least-squares / L-curve / Morozov loop - the same homogenization"
    "  cell-problem mathematics as effective-medium metamaterial design."
    "- result manifests from either side pass acoustic_fembem.result_manifest_gate"
    "  so scripts, MCP tools, and reports can ingest FEM/BEM and optimization"
    "  runs uniformly."
    ""
    "That is the point of integrating Gypsilab here: a FEM/BEM lab's"
    "MATLAB MCP should know the solver AND the design loop that drives it."
    ], newline);
end


function s = MULTIPHYSICS()
s = strjoin([
    "# Multiphysics coupling: the interface is the battlefield"
    ""
    "The individual physics are each well-understood; the DIFFICULTY"
    "concentrates at the INTERFACE where they meet. A coupled solver's bugs,"
    "convention mismatches, and accuracy bottlenecks live in the COUPLING,"
    "not the components (all five hit in the acoustic FSI build)."
    ""
    "1. INTERFACE OPERATOR = the whole coupling content. Vector displacement"
    "   <-> scalar pressure bridged by G_ij = int_Gamma mu_i (n . phi_j) - a"
    "   geometric pairing (a scalar trace with a vector normal component);"
    "   sign / normal-direction / convention bugs nest here."
    "2. DIFFERENT MATHEMATICAL HOMES. Acoustics = de Rham 0-form; elasticity"
    "   = the ELASTICITY (BGG) complex (rank-4 stress). The premetric picture"
    "   unifies the SETUP (d = topology metric-free; Star = material+metric),"
    "   but the elasticity complex is genuinely bigger than de Rham - the"
    "   coupling bridges two different function spaces."
    "3. CONVENTION RECONCILIATION. Outward normal, pressure-load sign,"
    "   e^{+ikr} vs e^{-ikr}, PV jump: each physics carries its own"
    "   conventions; the coupling works only when they are reconciled at"
    "   Gamma (the DtN outward normal had to match the interface-G normal)."
    "4. THE WEAKER PHYSICS IS THE BOTTLENECK. Coupled accuracy is limited by"
    "   the less-resolved side, not the coupling: the FSI residual ~10% is"
    "   the P1 INTERIOR elasticity (shear wavelength ~ sphere size at kR=2);"
    "   the coupling itself is exact (stiff->rigid 3.8e-3, DtN-vs-BEM 7e-4)."
    "5. BUGS COMPOUND AND HIDE. The 574%-wrong stiff limit came from the"
    "   INTERACTION of the fast DtN's low rank with the coupled block system"
    "   (p_s under-determined in the mesh null space) - neither the DtN"
    "   (exact 2.6e-5) nor the elasticity was wrong ALONE. And it HID: the"
    "   elastic case passed by min-norm luck while the stiff case failed."
    "   Fix = represent the exterior in its harmonic basis (full rank)."
    "   LESSON: a multiphysics bug lives in the coupling and shows only in a"
    "   specific regime - probe multiple regimes (stiff AND soft)."
    ""
    "NON-CONFORMING MESH COUPLING (hard even for the SAME physics):"
    "The interface has a SECOND difficulty independent of the physics - if"
    "the two sides of Gamma carry DIFFERENT meshes (non-matching nodes), the"
    "coupling int_Gamma (.)(.) dS is no longer a same-triangulation product."
    "The acoustic FSI here DODGES it by construction: the BEM surface IS the"
    "boundary of the FEM volume mesh (surface.volNodeIds), so mu_i and phi_j"
    "share ONE triangulation and G is exact. Real production FSI (a"
    "structural mesh + an independent fluid mesh) does NOT get to dodge it."
    "TWO WEAK-COUPLING routes (impose the interface condition weakly):"
    "  - MORTAR (Bernardi-Maday-Patera 1993; named for the binding agent -"
    "    it 'mortars' non-matching subdomains together): a Lagrange"
    "    multiplier enforces weak continuity - conservative, but a"
    "    SADDLE-POINT system whose multiplier space needs an inf-sup"
    "    condition."
    "  - NITSCHE (Joachim Nitsche 1971 - the mathematician, NOT the"
    "    philosopher Nietzsche): PRIMAL, NO multiplier - a consistency term"
    "    + a symmetric adjoint term + a penalty term. Variationally"
    "    consistent (unlike a raw penalty), SPD, optimal; the penalty must"
    "    be big enough for coercivity (too big -> ill-conditioned). The"
    "    backbone of DG and CutFEM/unfitted coupling, the modern FSI"
    "    favorite (no saddle point)."
    "TWO TRANSFER routes (move a trace across):"
    "  - SUPERMESH / COMMON REFINEMENT: integrate int_Gamma on the two"
    "    meshes' intersection mesh - conservative + consistent but"
    "    geometrically expensive (polygon clipping)."
    "  - POINTWISE INTERPOLATION: cheap but NON-CONSERVATIVE - can violate"
    "    the interface energy/flux balance, the quiet source of a slow drift."
    "A CONFORMING interface (shared boundary nodes, our case) sidesteps ALL"
    "of it - which is why the FEM-volume-boundary = BEM-surface is the clean"
    "teaching setup."
    ""
    "MANAGEMENT (why it stays tractable):"
    "- GEOMETRY makes the coupling principled: the interface is a geometric"
    "  pairing (de Rham 0-form trace <-> elasticity-complex traction), not"
    "  ad-hoc glue."
    "- VALIDATE THE COUPLING, not just the parts: gate against a reference"
    "  that genuinely couples (the Faran elastic-sphere resonance exists ONLY"
    "  because the elastic interior couples to the acoustic exterior), with"
    "  an INDEPENDENT reference (the earlier sign bug passed the"
    "  self-consistent limit, failed Anderson) and MULTIPLE regimes (stiff"
    "  caught the 574% the soft case hid)."
    "- RIGHT TOOL PER PART + WEAK COUPLING at the interface: FEM for the"
    "  complex bounded interior, exact operator / BEM for the open exterior"
    "  - do NOT force one method (a giant PML air box) to do everything."
    "  Lab EM analogues share the pattern AND the management: IH"
    "  electromagnetic+thermal, motor electromagnetic+mechanical, maglev"
    "  Radia-IEM + NGSolve-FEM weak coupling (no re-meshing)."
    ], newline);
end
