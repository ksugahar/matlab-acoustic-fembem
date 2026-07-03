"""Export NGSolve BEM + FEM references for the sonic-crystal sphere chain.

    python validation/exportNgsolveChainReference.py

Two independent NGSolve legs on the soft_sphere_chain_5 fixture (five
sound-soft spheres, R = 0.3, d = 1.5), at k = 1.0 (sub-wavelength
metamaterial regime) and k = 2.0944 (the chain's Bragg wavenumber pi/d -
where the free-space chain measurably does NOT develop a stop band):

  BEM   ngsolve.bem HelmholtzSingleLayerPotentialOperator on the SAME .vol
        and continuous-P1 space as the MATLAB side (intorder 16, 12-vs-16
        self-convergence stored), point-source + plane-wave solves, probe
        values via GetPotential on auxiliary probe meshes.
  FEM   volume Helmholtz FEM, an entirely different discretization: air
        ball R_out = 6 with the five spheres subtracted (OCC), scattered
        field with p_s = -p_inc on the sphere surfaces (Dirichlet),
        first-order Sommerfeld ABC (du/dn - i k u = 0) on the outer sphere,
        order 2, curved; total field evaluated at the same probe points.

The FEM leg deviates from BEM by ABC reflection + volume discretization;
the measured gap is recorded in the .mat consumers, not assumed.
"""
import time

import numpy as np
import scipy.sparse as sp
from scipy.io import savemat

import ngsolve
import ngsolve.bem as bem
from netgen.occ import Box, OCCGeometry, Pnt, Sphere
from ngsolve import (BND, H1, BilinearForm, GridFunction, LinearForm, Mesh,
                     TaskManager, ds, dx, exp, grad, x, y, z)

VOL = "fixtures/mesh_topology/soft_sphere_chain_5.vol"
OUT = "validation/data/ngsolve_chain_reference_soft_sphere_chain_5.mat"
RADIUS = 0.3
CENTERS = [(0.0, 0.0, zc) for zc in (-3.0, -1.5, 0.0, 1.5, 3.0)]
SOURCE_POINT = np.array([0.0, 0.0, 0.05])       # inside the middle sphere
PROBE_POINTS = np.array([[0.0, 0.0, 4.2],
                         [0.0, 0.0, 5.0],
                         [0.6, 0.0, 4.5]])
WAVENUMBERS = (1.0, 2.0944)
FEM_OUTER_RADIUS = 6.0
FEM_ORDER = 2
FEM_MAXH = 0.4
FEM_SCATTERER_MAXH = 0.12


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


def bem_leg(data):
    mesh = Mesh(VOL)
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
    print(f"BEM: ndof {fes.ndof}, boundary {bndmask.sum()}")

    for k in WAVENUMBERS:
        case = {"k": float(k)}
        t0 = time.perf_counter()
        Vc = np.array(bem.HelmholtzSingleLayerPotentialOperator(
            fes, fes, kappa=k, intorder=12).mat.ToDense(), copy=True)
        Vop = bem.HelmholtzSingleLayerPotentialOperator(
            fes, fes, kappa=k, intorder=16)
        Vd = np.array(Vop.mat.ToDense(), copy=True)
        case["V"] = Vd
        case["intorderConvergenceV"] = np.linalg.norm(Vc - Vd) / np.linalg.norm(Vd)
        print(f"  k={k}: BEM assembled {time.perf_counter() - t0:.0f}s, "
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
            case[f"g{tag}"] = g
            case[f"q{tag}"] = q
            case[f"probe{tag}"] = probe_potential(Vop.GetPotential(gf),
                                                  PROBE_POINTS)
        data["case_k" + str(k).replace(".", "p")] = case


def fem_leg(data):
    solid = Sphere(Pnt(0, 0, 0), FEM_OUTER_RADIUS)
    for c in CENTERS:
        solid = solid - Sphere(Pnt(*c), RADIUS)
    for face in solid.faces:
        if face.mass > 10.0:                      # outer: area 4*pi*36
            face.name = "outer"
        else:                                     # scatterers: area ~1.13
            face.name = "scatterer"
            face.maxh = FEM_SCATTERER_MAXH
    geo = OCCGeometry(solid)
    mesh = Mesh(geo.GenerateMesh(maxh=FEM_MAXH))
    mesh.Curve(FEM_ORDER)
    fes = H1(mesh, order=FEM_ORDER, complex=True, dirichlet="scatterer")
    print(f"FEM: ne {mesh.ne}, ndof {fes.ndof}")
    data["femNdof"] = float(fes.ndof)
    data["femNe"] = float(mesh.ne)
    data["femOrder"] = float(FEM_ORDER)
    data["femMaxh"] = FEM_MAXH
    data["femOuterRadius"] = FEM_OUTER_RADIUS
    data["femBc"] = "first_order_sommerfeld_abc"

    u, v = fes.TnT()
    for k in WAVENUMBERS:
        t0 = time.perf_counter()
        pinc = exp(1j * k * z)
        a = BilinearForm(fes)
        a += (grad(u) * grad(v) - k * k * u * v) * dx
        a += (-1j * k) * u * v * ds("outer")
        a.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(-pinc, BND, definedon=mesh.Boundaries("scatterer"))
        res = (-a.mat * gfu.vec).Evaluate()
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * res
        vals = np.array([gfu(mesh(*pt)) + np.exp(1j * k * pt[2])
                         for pt in PROBE_POINTS], dtype=complex)
        key = "case_k" + str(k).replace(".", "p")
        data[key]["femProbePlaneWave"] = vals
        print(f"  k={k}: FEM solved {time.perf_counter() - t0:.0f}s, "
              f"probes {np.abs(vals)}")


def main():
    data = {
        "ngsolveVersion": ngsolve.__version__,
        "intorder": 16.0,
        "sourcePoint": SOURCE_POINT.reshape(1, 3),
        "probePoints": PROBE_POINTS,
        "radius": RADIUS,
        "centers": np.array(CENTERS),
    }
    with TaskManager():
        bem_leg(data)
        fem_leg(data)
    savemat(OUT, data, do_compression=True)
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
