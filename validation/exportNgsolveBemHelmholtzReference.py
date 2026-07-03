"""Export ngsolve.bem HELMHOLTZ reference data for one .vol fixture.

    python validation/exportNgsolveBemHelmholtzReference.py \
        fixtures/mesh_topology/unit_sphere_coarse.vol \
        validation/data/ngbem_helmholtz_reference_unit_sphere_coarse.mat

Per wavenumber k (default 0.5 and 2.0, both away from the unit sphere's
first interior Dirichlet eigenvalue kR = pi where the first-kind V_k
equation goes singular), this stores a case struct with:

  V            complex Galerkin single layer on continuous P1
               (HelmholtzSingleLayerPotentialOperator, Sauter-Schwab,
               intorder 16; e^{+ikr} kernel - measured to MATCH the MATLAB
               teaching convention, lambda_0 = sin(k) e^{+ik}/k)
  intorderConvergenceV   12-vs-16 self-convergence, stored so a stale or
               hand-edited artifact fails loudly in the MATLAB gate
  gPointSource / qPointSource / probePointSource
               exterior Dirichlet solve with boundary data from an interior
               point source (the EXACT analytic gate) and NGSolve's own
               GetPotential values at the probe points - a solution-level
               reference fully independent of the MATLAB evaluator
  gPlaneWave / qPlaneWave / probePlaneWave
               sound-soft plane-wave scattering data (g = -exp(ikz)),
               compared against the analytic partial-wave series

Boundary data uses NODAL interpolation (H1 order-1 dof k = .vol point k+1),
matching the MATLAB side exactly, so density and probe values compare
without projection ambiguity. Requires NGSolve >= 6.2.2604 and scipy.
"""
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.io import savemat

import ngsolve
import ngsolve.bem as bem
from netgen.occ import Box, OCCGeometry, Pnt
from ngsolve import BND, H1, BilinearForm, GridFunction, Mesh, TaskManager, ds

SOURCE_POINT = np.array([0.3, 0.2, -0.25])   # inside the unit sphere
PROBE_POINTS = np.array([[2.0, 0.0, 0.0],
                         [0.0, 0.0, 3.0],
                         [-1.2, 1.6, 0.0]])


def point_source(k, x0, pts):
    r = np.linalg.norm(pts - x0, axis=1)
    return np.exp(1j * k * r) / (4 * np.pi * r)


def probe_potential(pot, points):
    values = np.zeros(len(points), dtype=complex)
    for i, pt in enumerate(points):
        geo = OCCGeometry(Box(Pnt(*(pt - 0.2)), Pnt(*(pt + 0.2))))
        pmesh = Mesh(geo.GenerateMesh(maxh=0.4))
        values[i] = pot(pmesh(*pt))
    return values


def export(vol_file, out_mat, wavenumbers=(0.5, 2.0), intorder=16, check_intorder=12):
    data = {
        "ngsolveVersion": ngsolve.__version__,
        "intorder": float(intorder),
        "sourcePoint": SOURCE_POINT.reshape(1, 3),
        "probePoints": PROBE_POINTS,
    }
    with TaskManager():
        mesh = Mesh(vol_file)
        fes = H1(mesh, order=1, complex=True)
        u, v = fes.TnT()
        mass = BilinearForm(fes, check_unused=False)
        mass += u * v * ds
        mass.Assemble()
        rows, cols, vals = mass.mat.COO()
        M = np.asarray(sp.csr_matrix(
            (vals, (rows, cols)), shape=(fes.ndof, fes.ndof)).todense())
        data["M"] = M.real

        vtx = np.array([p.point for p in mesh.vertices])
        data["vtx"] = vtx
        bndmask = np.zeros(fes.ndof, dtype=bool)
        for el in fes.Elements(BND):
            for d in el.dofs:
                bndmask[d] = True

        print(f"{vol_file}: ndof {fes.ndof}, boundary {bndmask.sum()}")
        for k in wavenumbers:
            case = {"k": float(k)}
            t0 = time.perf_counter()
            Vc = np.array(bem.HelmholtzSingleLayerPotentialOperator(
                fes, fes, kappa=k, intorder=check_intorder).mat.ToDense(), copy=True)
            Vop = bem.HelmholtzSingleLayerPotentialOperator(
                fes, fes, kappa=k, intorder=intorder)
            Vd = np.array(Vop.mat.ToDense(), copy=True)
            case["V"] = Vd
            case["intorderConvergenceV"] = (
                np.linalg.norm(Vc - Vd) / np.linalg.norm(Vd))
            print(f"  k={k}: assembled in {time.perf_counter() - t0:.1f}s, "
                  f"12-vs-16 {case['intorderConvergenceV']:.3e}")

            Vb = Vd[np.ix_(bndmask, bndmask)]
            Mb = M[np.ix_(bndmask, bndmask)]
            for tag, g_full in (
                    ("PointSource", point_source(k, SOURCE_POINT, vtx)),
                    ("PlaneWave", -np.exp(1j * k * vtx[:, 2]))):
                g = np.where(bndmask, g_full, 0.0)
                q = np.zeros(fes.ndof, dtype=complex)
                q[bndmask] = np.linalg.solve(Vb, Mb @ g[bndmask])
                gf = GridFunction(fes)
                gf.vec.FV().NumPy()[:] = q
                pot = Vop.GetPotential(gf)
                case[f"g{tag}"] = g
                case[f"q{tag}"] = q
                case[f"probe{tag}"] = probe_potential(pot, PROBE_POINTS)
            data["case_k" + str(k).replace(".", "p")] = case

    savemat(out_mat, data, do_compression=True)
    print(f"  saved {out_mat}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python exportNgsolveBemHelmholtzReference.py "
                 "<volFile> <outMat> [k1,k2,...]")
    ks = tuple(float(s) for s in sys.argv[3].split(",")) if len(sys.argv) > 3 \
        else (0.5, 2.0)
    export(sys.argv[1], sys.argv[2], ks)
