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
%         acoustic, convolution_quadrature, public_acoustic_blog_lessons,
%         public_acoustic_nonboundary_10, sonic_crystal, adjoint_ad,
%         matlab_execution_policy, mathworks_agentic_toolkit,
%         vol_visualization, pde_vol_bridge, gmsh_artifact,
%         radia_ngsolve_crossval, ngsolve_bem_50, catalog_100,
%         vibroacoustic_drum, curved_vol_geometry,
%         validation_discipline, optimization_link,
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
    case {"convolution_quadrature", "cq", "time_domain", "time_domain_bem", ...
            "lubich", "retarded", "td_bem"}
        text = CONVOLUTION_QUADRATURE;
    case {"public_acoustic_blog_lessons", "acoustic_blog", ...
            "method_selection", "public_acoustic_modeling_lessons"}
        text = PUBLIC_ACOUSTIC_BLOG_LESSONS;
    case {"public_acoustic_nonboundary_10", "nonboundary_acoustics", ...
            "acoustic_nonboundary_10", "modeling_nonboundary_10"}
        text = PUBLIC_ACOUSTIC_NONBOUNDARY_10;
    case {"sonic_crystal", "band_gap", "bloch", "metamaterial", "duct"}
        text = SONIC_CRYSTAL;
    case {"validation_discipline", "validation", "gates", "cross_check"}
        text = VALIDATION_DISCIPLINE;
    case {"radia_ngsolve_crossval", "radia-ngsolve", "ngsolve_crossval", ...
            "vol_crossval", "vol"}
        text = RADIA_NGSOLVE_CROSSVAL;
    case {"ngsolve_bem_50", "bem_50", "ngsolve_bem", "bem_crossval_50", ...
            "fifty"}
        text = NGSOLVE_BEM_50;
    case {"catalog_100", "hundred", "100", "gyp_100"}
        text = CATALOG_100;
    case {"vibroacoustic_drum", "drum", "taiko", "radiating_membrane", ...
            "baffled_drum"}
        text = VIBROACOUSTIC_DRUM;
    case {"curved_vol_geometry", "curved_vol", "high_order_curve", ...
            "geometry_order"}
        text = CURVED_VOL_GEOMETRY;
    case {"pde_vol_bridge", "pde_toolbox", "generate_mesh", "matlab_mesh"}
        text = PDE_VOL_BRIDGE;
    case {"vol_visualization", "visualization", "netgen_viewer", "vol_viewer"}
        text = VOL_VISUALIZATION;
    case {"gmsh_artifact", "gmsh", "node_data", "vol_fembem_gmsh", ...
            "artifact_contract"}
        text = GMSH_ARTIFACT;
    case {"matlab_execution_policy", "execution_policy", "no_live_documents", ...
            "scripts", "mcp_json"}
        text = MATLAB_EXECUTION_POLICY;
    case {"mathworks_agentic_toolkit", "agentic_toolkit", ...
            "official_matlab_mcp", "matlab_mcp_server", "mcp_runtime"}
        text = MATHWORKS_AGENTIC_TOOLKIT;
    case {"optimization_link", "optimization", "inverse", "design"}
        text = OPTIMIZATION_LINK;
    case {"multiphysics", "interface", "coupling_difficulty", ...
            "nonconforming", "mortar", "mesh_coupling"}
        text = MULTIPHYSICS;
    case "all"
        text = strjoin([OVERVIEW, SPACES, GALERKIN_BEM, COUPLED_FEM_BEM, ...
            MULTIPHYSICS, ACOUSTIC, CONVOLUTION_QUADRATURE, ...
            PUBLIC_ACOUSTIC_BLOG_LESSONS, ...
            PUBLIC_ACOUSTIC_NONBOUNDARY_10, ...
            SONIC_CRYSTAL, ADJOINT_AD, ...
            MATLAB_EXECUTION_POLICY, VOL_VISUALIZATION, PDE_VOL_BRIDGE, ...
            MATHWORKS_AGENTIC_TOOLKIT, ...
            GMSH_ARTIFACT, ...
            RADIA_NGSOLVE_CROSSVAL, NGSOLVE_BEM_50, CATALOG_100, ...
            VIBROACOUSTIC_DRUM, CURVED_VOL_GEOMETRY, ...
            VALIDATION_DISCIPLINE, ...
            OPTIMIZATION_LINK], [newline newline]);
    otherwise
        text = "Unknown topic '" + topic + "'. Available: overview, " + ...
            "spaces, galerkin_bem, coupled_fem_bem, multiphysics, acoustic, " + ...
            "convolution_quadrature, " + ...
            "public_acoustic_blog_lessons, public_acoustic_nonboundary_10, " + ...
            "sonic_crystal, adjoint_ad, " + ...
            "matlab_execution_policy, mathworks_agentic_toolkit, " + ...
            "vol_visualization, pde_vol_bridge, gmsh_artifact, " + ...
            "radia_ngsolve_crossval, ngsolve_bem_50, catalog_100, " + ...
            "vibroacoustic_drum, curved_vol_geometry, " + ...
            "validation_discipline, optimization_link, all.";
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
    "acoustic, public_acoustic_blog_lessons, public_acoustic_nonboundary_10,"
    "sonic_crystal,"
    "matlab_execution_policy, vol_visualization,"
    "pde_vol_bridge, gmsh_artifact, radia_ngsolve_crossval, ngsolve_bem_50,"
    "catalog_100, vibroacoustic_drum, curved_vol_geometry,"
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
    "ARBITRARY SHAPE (a box, not just the analytic sphere): the SAME FSI runs on a"
    "non-separable RECTANGULAR box via the BEM exterior (no truncation - the only"
    "mesh is the scatterer surface).  A box has no analytic reference, so its"
    "goldens are SHAPE-INDEPENDENT physics: far-field reciprocity"
    "f(x_hat; d) = f(-d; -x_hat) (to discretisation error ~1.5%) and the Sommerfeld"
    "1/r radiation decay (testElasticBoxScattering).  Box .vol comes from"
    "structuredBoxVol (no PDE Toolbox) or writePdeGeometryVol (any PDE geometry)."
    "DECISION RULE (FEM the scatterer?): FEM the interior ONLY when it carries"
    "physics BEM cannot - elastic (this FSI) or inhomogeneous fluid.  Rigid/soft ="
    "pure BEM (a surface BC; FEM would mesh + truncate the exterior, what BEM avoids)."
    ""
    "FAST EXTERIOR (spherical DtN, not a Kelvin acoustic boundary): for a SPHERE truncation,"
    "ExteriorMethod=""dtn"" swaps the dense Galerkin single/double layer for"
    "the EXACT spherical Helmholtz DtN (sphericalDtnOperator, Lambda_n ="
    "k h_n'(kR)/h_n(kR)). The scattered field is a spherical-"
    "harmonic expansion p_s = Phi c, so the exterior reduces to (N+1)^2"
    "coefficients and the coupled system stays FULL RANK (a nodal low-rank"
    "DtN leaves p_s under-determined - the bug that made the naive nodal swap"
    "574% wrong; the harmonic-coefficient reduction is the fix). DtN operator"
    "exact per multipole (independent point-source D->N check 2.6e-5 at"
    "degree 10), coupled field matches the rigid sphere to 3.8e-3 (BEM leg"
    "4.2e-3; DtN-vs-BEM 7e-4), singular N^2 assembly (~85 s) skipped ->"
    "sub-second solve. FAIL-LOUD on a non-sphere. The acoustic instance of"
    "the lab exact-operator (NOT absorber, NOT Kelvin boundary) open-boundary"
    "paradigm. Run it:"
    "acoustic_fembem.fembem_acoustic_gate(""fsi_dtn"")."
    ""
    "This SCALAR solve is the precursor. Lab EM analogues where coupling is the"
    "workhorse: nonlinear iron + open boundary, eddy currents + open"
    "boundary - but the acoustic lane should say BEM, Sommerfeld, DtN, or"
    "high-order impedance boundary, never Kelvin boundary."
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


function s = CONVOLUTION_QUADRATURE()
s = strjoin([
    "# Convolution quadrature (Lubich CQ) - the time-domain acoustic BEM"
    ""
    "CQ is the canonical time-domain lane, and for acoustics it is the STRONGEST"
    "route: it inherits the Laplace-domain kernel's A-stability, so it stays"
    "causal and late-time stable where direct marching-on-in-time (MOT) is"
    "notoriously unstable.  The retarded single/double layer are never formed as"
    "time kernels - the BDF generating function delta(zeta) is sampled on a"
    "contour, the Laplace-domain BEM V(s), K(s) are evaluated at s = delta(zeta)/dt,"
    "and an FFT recovers the time sequence."
    ""
    "THE COMPACTNESS IS STRUCTURAL: CQ = (frequency-domain operator) x (a thin,"
    "generic FFT wrapper).  The hard physics (the BEM kernel, the FEM/BEM coupling"
    "block) is the SAME frequency-domain code, evaluated at complex frequencies."
    "So you validate ONCE in the frequency domain and the time-domain solver"
    "inherits its correctness plus CQ's stability.  volTdBemConvolutionQuadrature"
    "(exterior single-layer Dirichlet) and volFemBemCoupledConvolutionQuadrature"
    "(volume-FEM/BEM coupled) share ONE goldened operator, laplaceSingleLayerGalerkin."
    ""
    "THE GOLDEN (the CQ correctness anchor): the CQ core is the Laplace-domain"
    "single layer V(s), kernel exp(-s r/c)/(4 pi r).  On the IMAGINARY axis"
    "s = -1i c k it becomes exp(1i k r) = the Helmholtz kernel, so"
    "  laplaceSingleLayerGalerkin(surface, -1i c k, c, q)"
    "     == GalerkinSingleLayer(surface, k, q).matrix + Delta"
    "to machine precision (measured relerr ~3e-17, all quad orders and k).  This"
    "pins the s/c scaling, the exponent sign, and the 1/(4 pi) to the analytically-"
    "validated (Faran/Anderson/soft-sphere) frequency-domain single layer - far"
    "stronger than a self-consistent CQ residual (all the old smoke checks tested)."
    ""
    "THE COINCIDENT-NODE CONVENTION (a real subtlety, verified): the smooth"
    "correction (exp(-alpha r)-1)/(4 pi r) is a REGULAR integrand with the finite"
    "limit -alpha/(4 pi) at coincident quadrature points.  Keeping that limit is"
    "the correct product-Gauss sample and the MORE accurate choice;"
    "GalerkinSingleLayer/HelmholtzKernel instead ZERO the smooth diagonal, so the"
    "two differ by the KNOWN term Delta_ij = (-alpha/4 pi) sum_g w_g^2 phi_i phi_j."
    "Do NOT bend the CQ to zero it just to make the golden exact (that moves the CQ"
    "~1e-2..3e-1 away from the validated operator) - keep -alpha and add Delta on"
    "the reference side.  Golden: testLaplaceSingleLayerGalerkin."
    ], newline);
end


function s = PUBLIC_ACOUSTIC_BLOG_LESSONS()
s = strjoin([
    "# Public acoustic modeling lessons -> readable FEM/BEM method selection"
    ""
    "This topic distills public acoustic modeling material into a"
    "solver-independent teaching checklist.  It is not a solver benchmark"
    "and it does not import proprietary models or numerical values.  The"
    "lesson is how to choose the modeling family before writing a Gypsilab-"
    "style MATLAB script."
    ""
    "1. Unbounded exterior radiation: use BEM or an exact DtN/open-boundary"
    "   operator in the frequency domain when the exterior is homogeneous."
    "   This repository's policy for wave boundaries is high-order surface"
    "   impedance Zs; do not use PML as the default validation route for"
    "   acoustic or electromagnetic waves."
    "   In this repository, the readable route is .vol boundary triangles ->"
    "   P1 BEM; the production reference is NGSolve.BEM/radia-ngsolve."
    "   Do not claim a generic time-domain BEM just because the frequency-"
    "   domain BEM exists.  Our time-domain lane is explicitly CQ:"
    "   volTdBemConvolutionQuadrature / volFemBemCoupledConvolutionQuadrature."
    ""
    "2. Acoustic-structure interaction: the interface must carry both"
    "   structural displacement/normal velocity and acoustic pressure/normal"
    "   traction.  A drum, bell, transducer, membrane, or elastic bead is a"
    "   FEM interior plus acoustic exterior problem, not a painted pressure"
    "   animation.  The teaching gate is two-way coupling plus a result"
    "   artifact schema for the interface convention."
    ""
    "3. Impedance lumping: a detailed thermo/viscous or resonator submodel"
    "   can be replaced by a frequency-dependent impedance only if the"
    "   impedance definition travels with the power balance.  Record whether"
    "   it is specific impedance p=Z_s v or acoustic impedance p=Z Q, and"
    "   keep incident/reflected/transmitted/dissipated power explicit."
    ""
    "4. Absorbing boundaries: local reaction and extended reaction are"
    "   different modeling assumptions, but they are metadata below the lab"
    "   boundary policy.  For acoustic and electromagnetic waves, record"
    "   high-order Zs as the boundary model and explicitly record PML=false."
    "   Extended-reaction porous or interior-impedance physics can inform the"
    "   fitted Zs, but a scalar absorption coefficient alone is not a reusable"
    "   validation artifact."
    ""
    "5. Room acoustics: choose by acoustic/geometric scale.  Below the"
    "   Schroeder-frequency region, wave/FEM or eigenmode analysis teaches"
    "   modal behavior.  At high frequency, ray or diffusion models are"
    "   efficient but no longer resolve wave effects such as diffraction and"
    "   standing waves.  A hybrid example should record the split frequency."
    ""
    "MCP/radia bridge: use radia-mcp acoustic_method_selection_manifest_gate"
    "before promoting any public acoustic literature lesson into a"
    "Gypsilab example, report, or radia-ngsolve validation candidate."
    ], newline);
end


function s = PUBLIC_ACOUSTIC_NONBOUNDARY_10()
s = strjoin([
    "# 10 public acoustic non-boundary problems for Gypsilab/radia-acoustic"
    ""
    "This catalog deliberately excludes boundary-condition topics such as"
    "absorbing boundaries, PML, port boundaries, and impedance-boundary-only"
    "examples.  When an open wave closure is later needed, this project's policy"
    "still records high-order Zs and PML=false, but these ten lessons are"
    "about the acoustic physics inside the model."
    ""
    "1. Acoustic trap with Gor'kov potential, streaming, and particles:"
    "   Gypsilab task = P1 pressure plus a particle-force post map."
    "   radia-acoustic task = manifest gate for pressure extrema, force"
    "   direction, and particle equilibrium."
    ""
    "2. Surface-acoustic-wave droplet streaming:"
    "   Gypsilab task = SAW source field -> quadratic streaming forcing."
    "   radia-acoustic task = reduced streaming-force artifact with wave"
    "   direction, steady-flow sign, and mixing observable."
    ""
    "3. Thermoviscous acoustic radiation force:"
    "   Gypsilab task = inviscid Gor'kov fast gate, with thermoviscous"
    "   corrections marked as the heavy lane."
    "   radia-acoustic task = gate first-order-field provenance, particle"
    "   contrast, perturbation order, and force vector."
    ""
    "4. Small microphone with thermoviscous and electromechanical losses:"
    "   Gypsilab task = lumped diaphragm plus thermoviscous duct/cavity."
    "   radia-acoustic task = frequency-response artifact with acoustic"
    "   compliance, mechanical resonance, and loss-channel identity."
    ""
    "5. Thermoacoustic engine:"
    "   Gypsilab task = readable 1D transfer-matrix or weak-form gate."
    "   radia-acoustic task = heat-work sign, acoustic power flux, and"
    "   stack-location metadata."
    ""
    "6. Acoustic topology optimization with density/bulk-modulus"
    "   interpolation:"
    "   Gypsilab task = tiny Helmholtz design field with objective,"
    "   constraint, filter, and material interpolation metadata."
    "   radia-acoustic task = design-variable bounds and objective-region"
    "   gate."
    ""
    "7. Microacoustic topology optimization with thermoviscous losses:"
    "   Gypsilab task = lossy equivalent-fluid objective before full"
    "   thermoviscous elements."
    "   radia-acoustic task = loss-model identity and dissipated-power"
    "   columns in the result artifact."
    ""
    "8. Room response split into wave and statistical/high-frequency"
    "   methods:"
    "   Gypsilab task = Schroeder-split teaching manifest."
    "   radia-acoustic task = low-band method, high-band method, transition"
    "   frequency, and merge-rule gate."
    ""
    "9. Small-speaker room impulse response with FEM-to-ray source handoff:"
    "   Gypsilab task = save source directivity/near-field data before"
    "   room propagation."
    "   radia-acoustic task = receiver list, source-map id, and time-axis"
    "   convention gate."
    ""
    "10. Ultrasonic pipe pulse-echo reduced replay:"
    "    Gypsilab task = symmetry-reduced pulse/echo manifest with material"
    "    wave speeds."
    "    radia-acoustic task = pulse definition, echo windows, and arrival-"
    "    ordering gate."
    ""
    "Reusable gate: radia-mcp public_acoustic_nonboundary_problem_catalog"
    "returns the same ten case families, and"
    "acoustic_nonboundary_problem_catalog_manifest_gate rejects boundary"
    "families, PML use, missing Gypsilab/radia tasks, and missing observables."
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
    "GENERAL geometry (any PDE Toolbox shape, not just a box):"
    "  writePdeGeometryVol(""sphere.vol"", multisphere(0.5), Hmax=0.12)"
    "  writePdeGeometryVol(""cyl.vol"",    multicylinder(0.4, 1.0), Hmax=0.15)"
    "  writePdeGeometryVol(""part.vol"",   importGeometry(""part.step""), Hmax=0.05)"
    "writePdeBoxVol is a thin box wrapper over it.  In R2026a the geometry builders"
    "are BUILT-INS (exist == 5), so an availability check must use exist(name)==0,"
    "NOT exist(...,""file"")==2 (which wrongly rejects an installed toolbox)."
    ""
    "NO-TOOLBOX box: structuredBoxVol writes a box .vol with a structured"
    "6-tet-per-hex (Kuhn) mesh and NO PDE Toolbox at all - the box scatterer test"
    "uses it so it runs on a minimal checkout."
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
    "  - GUI/user inspection: Netgen native .vol viewer, but do not rely"
    "    on a raw Windows netgen.exe file association."
    "  - MATLAB script/report: plotVolMesh plus the summary JSON."
    "  - LLM/headless preflight: acoustic_fembem_vol_mesh_summary."
    ""
    "Windows double-click rule:"
    "  Raw .vol -> netgen.exe can open a blank GUI because the startup path"
    "  does not necessarily call the native Tcl mesh loader. Use the Radia"
    "  helper association (`radia-vol-viewer --register`) or an equivalent"
    "  Netgen startup hook that runs:"
    "    Ng_LoadMesh ""mesh.vol"""
    "    set selectvisual mesh"
    "    Ng_SetVisParameters; redraw; Ng_ReadStatus"
    "  .sol files are mesh-free GridFunction coefficient dumps, so double-click"
    "  handling must locate or receive the companion .vol and the matching FES"
    "  order before selecting solution visual mode."
    ""
    "Do not turn visualization into a hidden mesh conversion step. The"
    "solver-facing contract remains first-order tri/tet .vol, and the"
    "reader rejects quad/hex/wedge/pyramid/curved records fail-loud."
    ], newline);
end


function s = GMSH_ARTIFACT()
s = strjoin([
    "# Acoustic field movies: native MATLAB here; gmsh is radia-acoustic"
    ""
    "This MATLAB lane is gmsh-free.  Time-domain acoustic FIELD MOVIES are made"
    "INSIDE MATLAB as indexed-image animated GIFs (headless, no gmsh):"
    "  field = softSphereScatterField();  writeSoftSphereScatterGif(field, gif);"
    "  field = drumRollField();           writeSoftSphereScatterGif(field, gif);"
    "  field = drumScatterField();        writeSoftSphereScatterGif(field, gif);"
    "  field = drumStepTimeField();       writeDrumStepTimeGif(field, gif);"
    "The GIF writer maps the pressure array directly to an indexed image, so it"
    "works in batch/headless MATLAB and does not require Gmsh."
    ""
    "The gmsh route for acoustic field movies belongs to the radia-acoustic"
    "(Python) side, where GmshPostExport already writes .msh v4.1 NodeData time"
    "series (scalar pressure / vector displacement) for the gmsh animation"
    "player.  Division of labor: MATLAB teaching repo -> native MATLAB figures /"
    "GIF only; radia-acoustic -> gmsh via GmshPostExport (.msh v4.1).  Do not"
    "reintroduce a gmsh writer into this MATLAB lane; see the radia_mcp"
    "acoustic_fembem_cross_learnings ""Visualization"" section for the same rule."
    ], newline);
end


function s = MATHWORKS_AGENTIC_TOOLKIT()
s = strjoin([
    "# MathWorks MATLAB MCP Server and Agentic Toolkit policy"
    ""
    "Use the official MathWorks MATLAB MCP Server as the execution runtime."
    "It owns MATLAB process/session management, code evaluation, test running,"
    "Code Analyzer checks, and toolbox detection.  This repository must not"
    "vendor, fork, or wrap a second general MATLAB MCP server."
    ""
    "Use the MATLAB Agentic Toolkit as the setup and skills reference layer."
    "The toolkit can install/update the official MCP server and register"
    "agent skills.  Its role is agent guidance and configuration, not domain"
    "solver logic.  Keep only the skill groups that are relevant to this"
    "acoustic FEM/BEM project so tool selection stays reliable."
    ""
    "This repository contributes the thin domain extension:"
    "  - mcp/extensions/acoustic-fembem-tools.json is the MCP contract;"
    "  - +acoustic_fembem functions are ordinary MATLAB entry points;"
    "  - matlab_api contains the readable P1 FEM/BEM/CQ implementation;"
    "  - tests verify the extension without requiring an MCP client."
    ""
    "Existing-session policy: when sharing a MATLAB session with COMSOL"
    "LiveLink, attach to that already-running MATLAB session through the"
    "official server's existing-session workflow.  Do not start a second"
    "MATLAB only to reach COMSOL.  Pure MATLAB acoustic FEM/BEM checks may"
    "run in batch when they do not touch COMSOL."
    ""
    "Domain division: MathWorks runtime and skills teach MATLAB mechanics;"
    "the acoustic_fembem extension teaches .vol tri/tet intake, P1 volume FEM,"
    "P1 boundary BEM, Johnson-Nedelec coupling, Lubich CQ, native GIF"
    "visualization, and NGSolve.BEM/radia-ngsolve cross-validation gates."
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


function s = NGSOLVE_BEM_50()
s = strjoin([
    "# NGSolve.BEM 50-case comparison lane"
    ""
    "Yes: 50 examples are realistic before the full 100-case catalog. Use the"
    "BEM-heavy half of catalog_100 as the first milestone, all driven by the"
    "same first-order Netgen .vol tri/tet fixtures:"
    ""
    "  GYP-031..040  Laplace P1 BEM, dense assembly"
    "  GYP-041..050  Laplace P1 BEM, H-matrix/compression checks"
    "  GYP-051..060  acoustic low-frequency stability checks"
    "  GYP-061..070  acoustic Helmholtz P1 BEM"
    "  GYP-091..100  NGSolve.BEM reference and convention checks"
    ""
    "The .vol visualization skill for this lane is:"
    "  1. open the solver-facing .vol in Netgen for human inspection;"
    "  2. run acoustic_fembem_vol_mesh_summary for MCP/headless preflight;"
    "  3. use plotVolMesh only as a MATLAB quick-look figure."
    ""
    "Analytic references are not gone, but they are no longer enough by"
    "themselves. Keep the exact sphere/ball, Laplace capacity, point-source"
    "reproduction, reciprocity, symmetry, and convergence identities as"
    "anchors; for the remaining coverage, compare MATLAB/Gypsilab against"
    "radia-ngsolve/NGSolve.BEM and store result manifests with versions,"
    "run dates, timing breakdowns, schema ids, and convention ids."
    ""
    "Promotion rule: a 50-case result teaches the MCP only after mesh summary,"
    "operator convention gate, reference comparison, and manifest metadata all"
    "pass. Do not claim an analytic solution when the reference is a measured"
    "cross-code artifact."
    ], newline);
end


function s = CATALOG_100()
s = strjoin([
    "# 100-case acoustic FEM/BEM catalog"
    ""
    "The repository already has a runnable 100-case teaching catalog:"
    "  acoustic_fembem.fembem_crossval_gate(""catalog_100"")"
    ""
    "Coverage:"
    "  GYP-001..010  .vol mesh/topology, including negative quad/hex rejection"
    "  GYP-011..020  H1 P1 tetra FEM"
    "  GYP-021..030  HCurl/Nedelec edge FEM"
    "  GYP-031..040  Laplace P1 BEM dense operators"
    "  GYP-041..050  readable Laplace H-matrix blocks"
    "  GYP-051..060  acoustic low-frequency kernels"
    "  GYP-061..070  Helmholtz/acoustic BEM"
    "  GYP-071..080  scalar FEM/BEM coupling"
    "  GYP-081..090  RWG/HCurl trace maps"
    "  GYP-091..100  NGSolve/NGSolve.BEM reference smoke gates"
    ""
    "This catalog is a regression/teaching ladder, not the end of the story."
    "The next high-value examples should be more physical: baffled vibrating"
    "disk/drum radiation, vibro-acoustic FEM/BEM coupling, and radiation"
    "impedance checks."
    ], newline);
end


function s = VIBROACOUSTIC_DRUM()
s = strjoin([
    "# Vibro-acoustic drum / struck membrane example"
    ""
    "Yes, this is exactly the kind of FEM/BEM example the teaching lane should"
    "grow next.  The educational model is a baffled circular membrane or thin"
    "elastic plate: the structural FEM supplies the normal velocity on the"
    "radiating face, and the acoustic P1 BEM supplies the exterior radiation"
    "condition and pressure field."
    "Standing split for drums: the drum structure is FEM; the air radiation"
    "is acoustic BEM.  Do not make the drum body an acoustic volume-FEM"
    "problem unless the lesson is specifically an interior-air/cavity"
    "transmission problem."
    ""
    "Minimal readable rung:"
    "  1. Membrane/plate FEM eigenmode on a circular disk (start with mode 0,1)."
    "  2. Use the mode shape as prescribed normal velocity on the disk."
    "  3. Solve exterior Helmholtz P1 BEM for the radiated pressure."
    "  4. Check acoustic power/radiation impedance and far-field directivity."
    "  5. Compare with NGSolve.BEM on the same .vol/.surface labels."
    ""
    "Time-domain MATLAB visualization rung:"
    "  field = drumStepTimeField();"
    "  plotDrumStepTimeField(field, 30);"
    "  writeDrumStepTimeGif(field, fullfile(tempdir, ""drum_step_time_field.gif""));"
    "  scene = drumHighOrderImpedanceScene(field);"
    "  writeDrumHighOrderImpedanceGif(scene, fullfile(tempdir, ""drum_high_order_impedance.gif""));"
    "  realScene = drumFemBemCoupledDemo();"
    "  writeDrumHighOrderImpedanceGif(realScene, fullfile(tempdir, ""drum_fembem_coupled.gif""));"
    "  field = drumScatterField();  % drum + sphere scatterer, MATLAB-only, no gmsh"
    "  writeSoftSphereScatterGif(field, fullfile(tempdir, ""drum_scatter.gif""));"
    "This uses a step-force structural modal response plus the causal Rayleigh"
    "retarded-potential integral, then draws the r-z pressure snapshot inside"
    "MATLAB.  The GIF writer maps the pressure array directly to an indexed"
    "image, so it works in batch/headless MATLAB and does not require Gmsh."
    "The higher-context scene uses axis-equal x-z pixels, shows the full"
    "spherical truncation boundary (not a hemisphere), and treats the drum"
    "top membrane at z=0 as the struck surface.  The spherical truncation is"
    "labelled as the high-order impedance absorbing-boundary lane for"
    "Radia-style open-boundary experiments.  For acoustic waves, this"
    "high-order impedance boundary is mandatory; do not use or name a Kelvin"
    "boundary in the acoustic lane."
    "In this first rung the lower half-space is intentionally quiet: the"
    "cylindrical body is a rigid baffle, not a second radiating drum head."
    "The real-drum rung, drumFemBemCoupledDemo, removes that simplification:"
    "top membrane FEM mode, bottom membrane FEM mode, and shell side-leakage"
    "mode are damped oscillators with explicit damping ratios, then the"
    "exterior top/bottom/side radiation is evaluated by BEM-style retarded"
    "boundary integrals.  Do not model this teaching drum as an internal"
    "cavity-pressure oscillator unless a true acoustic cavity mesh is present."
    "The"
    "BEM layer must NOT split the observation field by source direction:"
    "top, bottom, and side boundary sources are all evaluated at every"
    "exterior air observation point with the same retarded Green kernel,"
    "then superposed.  A top source may wrap to the side/lower field and a"
    "side source may radiate upward/downward; direction-only painting is a"
    "visualization bug, not FEM/BEM coupling."
    "transient solve is a reduced FEM ODE integrated by ode45 with an impact"
    "pulse and damping-ratio terms; the exterior field is reconstructed from"
    "causal retarded boundary potentials.  It is a readable FEM/BEM coupling"
    "demo, not yet a full time-domain P1 volume-FEM/P1 surface-BEM solver."
    "The visualized color field is only the propagating air pressure outside"
    "the drum; the interior air volume is geometric only and is not a cavity"
    "pressure DOF.  The result intentionally shows lower-half radiation,"
    "side leakage, and decaying membrane/shell vibration."
    "The same modeling split can be implemented in NGSolve: membrane/shell"
    "FEM or modal damped oscillators provide Neumann data, ngsolve.bem"
    "provides the exterior radiation operator, and time response comes from"
    "frequency sweeps with inverse FFT or from an externally assembled CQ"
    "loop over Laplace-domain BEM operators."
    "Be precise about dimensionality: the drum GIF is a 3D axisymmetric"
    "physical picture rendered on an r-z slice.  Disk and side sources are"
    "integrated over azimuth, so it is not a 2D Cartesian wave cartoon, but"
    "it is not yet a full 3D structural-FEM/acoustic-BEM drum mesh."
    "The parallel acoustic-volume teaching lane is volFemBemIfftResponse:"
    "read a Netgen .vol tri/tet mesh, assemble H1/P1 acoustic volume FEM"
    "plus boundary P1 BEM, solve the frequency-domain Helmholtz FEM/BEM"
    "system on many frequency bins, multiply by a real pulse spectrum,"
    "enforce Hermitian symmetry, and inverse FFT to obtain a real pressure"
    "time history.  This is not a periodic sine-wave animation and not yet"
    "convolution-quadrature TD-BEM; it is the frequency-sweep/iFFT route for"
    "acoustic transmission and interior-fluid examples, not the preferred"
    "drum structural model."
    "The first convolution-quadrature TD-BEM rung is"
    "volTdBemConvolutionQuadrature.  It uses the .vol boundary triangles as"
    "a P1 Galerkin BEM surface, samples the BDF generating function"
    "delta(zeta) at CQ points, evaluates Laplace-domain single-layer matrices"
    "V(s) with Re(s)>0, solves V(s) qhat = ghat, and FFT-recovers q(t) plus"
    "the exterior pressure.  This is real Lubich CQ TD-BEM for the exterior"
    "Dirichlet single-layer problem; the production volume-FEM/interior"
    "coupling CQ version is the next step."
    "The production-form acoustic-volume/interior coupled CQ rung is"
    "volFemBemCoupledConvolutionQuadrature.  It uses the same .vol volume"
    "tetrahedra as H1/P1 interior wave FEM and the same boundary triangles as"
    "P1 BEM, then solves [A+(s/c1)^2 M, -T'Mb; (1/2 Mb-K(s))*T, V(s)]"
    "at every CQ Laplace point.  A volume source pulse drives the interior,"
    "and the BEM flux density radiates through the exterior representation"
    "-S(s)q + D(s)Tu.  This is now the Calderon/Johnson-Nedelec coupled CQ"
    "system with the retarded double-layer K(s).  CouplingForm="
    "SingleLayerTeaching keeps the old [Mb*T, -V(s)] rung for regression,"
    "but production examples should use the JohnsonNedelec default."
    "For MATLAB-side movies use the native GIF path above (drumStepTimeField /"
    "drumScatterField -> writeSoftSphereScatterGif); for a drum, replace the interior acoustic volume FEM with structural membrane/shell FEM and use"
    "the BEM unknowns for exterior sound.  Gmsh field movies of the acoustic-"
    "volume CQ result are produced on the radia-acoustic (Python) side via"
    "GmshPostExport (.msh v4.1 NodeData), not in this MATLAB lane.  Draw the"
    "sound field as a plane plot and let the drum boundary itself be the 3D"
    "object that moves; avoid the misleading sphere-in-plane picture."
    ""
    "This is better than another sphere-only case because it feels like a"
    "student experiment: hit a drum, watch the membrane mode, then hear the"
    "radiated sound.  It should become a post-100 catalog extension before"
    "performance tuning."
    ], newline);
end


function s = CURVED_VOL_GEOMETRY()
s = strjoin([
    "# Curve-only high-order .vol geometry policy"
    ""
    "Good idea, but keep it optional and isolated.  The readable teaching"
    "solver should keep first-order unknowns: H1 P1 for FEM and P1 boundary"
    "BEM.  A high-order .vol should be used only as a geometry/quadrature view"
    "for curved surfaces, not as hidden high-order solution DOFs."
    ""
    "Readable design:"
    "  - VolMesh remains the first-order tri/tet connectivity and node-id"
    "    contract."
    "  - CurvedBoundaryView, or an equivalent small adapter, owns high-order"
    "    boundary nodes, curved triangle mapping, normals, Jacobians, and"
    "    quadrature points."
    "  - FEM/BEM operators accept either SurfaceMesh or CurvedBoundaryView."
    "  - The default parser still fails loud on high-order records unless"
    "    EnableCurvedGeometry=true is explicitly requested."
    ""
    "This does lower readability if mixed into readVolTriTet or every operator."
    "It stays readable if all curved geometry logic lives behind one adapter and"
    "the examples show both: linear geometry first, curved geometry as the"
    "accuracy extension."
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
