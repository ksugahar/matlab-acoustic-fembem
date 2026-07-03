"""Stage 9: sonic-crystal band gap in a rigid duct (NGSolve FEM reference).

    python validation/exportNgsolveDuctBandReference.py

The confinement experiment the free-space chain (stage 8) proved
necessary: the SAME scatterer family (sound-soft sphere R = 0.3, lattice
constant d = 1.5) placed in a rigid-walled duct of cross-section a x a
(a = 1.0, first transverse cutoff pi/a ~ 3.14, so the working window
k < 3 is single-mode).

Part A - Bloch band structure of the duct unit cell (a x a x d minus the
sphere, Dirichlet on the sphere, rigid walls natural, z-faces
OCC-identified periodic with phase exp(-i q d); order-2 complex H1,
ArnoldiSolver). Analytic gate: the EMPTY cell's lowest bands are exactly
|q + 2 pi n / d| (plane duct modes below the transverse cutoff).

Part B - transmission through the finite crystal: duct of length L with
N = 5 spheres, port Robin conditions (inlet dp/dn - i k p = -2 i k,
outlet dp/dn - i k p = 0, exact for the single propagating mode),
plane-mode amplitude at the outlet vs k. Gate: the EMPTY duct must be
transparent (|T| ~ 1), and the transmission dip of the loaded duct must
sit inside the Bloch gap from part A.

Results land in validation/data/ngsolve_duct_band_reference.mat and are
locked by tests/testDuctBandGap.m on the MATLAB side.
"""
import time

import numpy as np
from scipy.io import savemat

import ngsolve
from netgen.occ import Box, IdentificationType, OCCGeometry, Pnt, Sphere, Z
from ngsolve import (BND, H1, ArnoldiSolver, BilinearForm, GridFunction,
                     Integrate, LinearForm, Mesh, Periodic, TaskManager, ds,
                     dx, exp, grad, x, y, z)

A_SIDE = 1.0
D_CELL = 1.5
RADIUS = 0.3
N_SPHERES = 5
PAD = 3.75
MAXH = 0.18
SPHERE_MAXH = 0.09
ORDER = 2
NEV = 6
QD_GRID = np.linspace(0.0, np.pi, 9)
K_GRID = np.arange(1.2, 2.91, 0.1)


def on_box_wall(c, length):
    eps = 1e-9
    return (abs(c.x) < eps or abs(c.x - A_SIDE) < eps
            or abs(c.y) < eps or abs(c.y - A_SIDE) < eps
            or abs(c.z) < eps or abs(c.z - length) < eps)


def cell_mesh(with_sphere):
    box = Box(Pnt(0, 0, 0), Pnt(A_SIDE, A_SIDE, D_CELL))
    solid = box
    if with_sphere:
        sph = Sphere(Pnt(A_SIDE / 2, A_SIDE / 2, D_CELL / 2), RADIUS)
        solid = box - sph
        for face in solid.faces:
            # planar wall faces have their center ON a wall plane; the
            # sphere face's center is the sphere centroid, strictly interior
            if not on_box_wall(face.center, D_CELL):
                face.name = "scatterer"
                face.maxh = SPHERE_MAXH
    solid.faces.Min(Z).Identify(solid.faces.Max(Z), "zper",
                                IdentificationType.PERIODIC)
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=MAXH))


def bloch_bands(mesh, dirichlet):
    bands = np.zeros((len(QD_GRID), NEV))
    for iq, qd in enumerate(QD_GRID):
        fes = Periodic(H1(mesh, order=ORDER, complex=True,
                          dirichlet=dirichlet),
                       phase=[np.exp(-1j * qd)])
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        m = BilinearForm(fes)
        m += u * v * dx
        m.Assemble()
        gf = GridFunction(fes, multidim=NEV * 3)
        lam = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
                            list(gf.vecs), shift=4.0)
        ks = np.sort(np.sqrt(np.abs(np.array([complex(l).real
                                              for l in lam]))))
        bands[iq, :] = ks[:NEV]
    return bands


def duct_transmission(with_spheres):
    box = Box(Pnt(0, 0, 0), Pnt(A_SIDE, A_SIDE,
                                2 * PAD + N_SPHERES * D_CELL))
    length = 2 * PAD + N_SPHERES * D_CELL
    solid = box
    if with_spheres:
        for j in range(N_SPHERES):
            zc = PAD + (j + 0.5) * D_CELL
            solid = solid - Sphere(Pnt(A_SIDE / 2, A_SIDE / 2, zc), RADIUS)
    for face in solid.faces:
        c = face.center
        if abs(c.z) < 1e-9:
            face.name = "inlet"
        elif abs(c.z - length) < 1e-9:
            face.name = "outlet"
        elif not on_box_wall(c, length):
            face.name = "scatterer"
            face.maxh = SPHERE_MAXH
    dirichlet = "scatterer" if with_spheres else ""
    mesh = Mesh(OCCGeometry(solid).GenerateMesh(maxh=MAXH))
    fes = H1(mesh, order=ORDER, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    print(f"duct with_spheres={with_spheres}: ne {mesh.ne}, ndof {fes.ndof}")

    T = np.zeros(len(K_GRID))
    for ik, k in enumerate(K_GRID):
        a = BilinearForm(fes)
        a += (grad(u) * grad(v) - k * k * u * v) * dx
        a += (-1j * k) * u * v * ds("inlet|outlet")
        a.Assemble()
        f = LinearForm(fes)
        f += (-2j * k) * v * ds("inlet")
        f.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                     inverse="pardiso") * f.vec
        amp = Integrate(gfu, mesh, BND,
                        definedon=mesh.Boundaries("outlet")) / A_SIDE**2
        T[ik] = abs(amp)
    return T


def main():
    data = {
        "ngsolveVersion": ngsolve.__version__,
        "aSide": A_SIDE, "dCell": D_CELL, "radius": RADIUS,
        "nSpheres": float(N_SPHERES), "order": float(ORDER),
        "maxh": MAXH, "qd": QD_GRID, "k": K_GRID,
    }
    with TaskManager():
        t0 = time.perf_counter()
        empty = bloch_bands(cell_mesh(False), "")
        data["bandsEmpty"] = empty
        # analytic empty duct-cell bands:
        # k = sqrt((m1^2+m2^2) pi^2/a^2 + ((qd + 2 pi n)/d)^2),
        # transverse families (0,0),(1,0),(0,1),(1,1), axial n = 0,-1,+1
        families = []
        for m2sum in (0.0, 1.0, 1.0, 2.0):
            for n in (0, -1, 1):
                families.append(np.sqrt(
                    m2sum * (np.pi / A_SIDE)**2
                    + ((QD_GRID + 2 * np.pi * n) / D_CELL)**2))
        ana = np.sort(np.stack(families, axis=1), axis=1)
        data["bandsEmptyAnalytic"] = ana
        err = np.abs(empty[:, :3] - ana[:, :3]) / np.abs(ana[:, :3] + 1e-30)
        # skip q = 0 first band (k = 0 constant mode)
        data["emptyLatticeMaxRelErr"] = float(np.max(err[1:, :]))
        print(f"empty-lattice gate: max rel err {data['emptyLatticeMaxRelErr']:.3e} "
              f"({time.perf_counter() - t0:.0f}s)")

        t0 = time.perf_counter()
        loaded = bloch_bands(cell_mesh(True), "scatterer")
        data["bands"] = loaded
        gap_low = float(np.max(loaded[:, 0]))
        gap_high = float(np.min(loaded[:, 1]))
        data["gapLow"] = gap_low
        data["gapHigh"] = gap_high
        print(f"Bloch bands: gap [{gap_low:.4f}, {gap_high:.4f}] "
              f"({time.perf_counter() - t0:.0f}s)")

        t0 = time.perf_counter()
        data["transmissionEmpty"] = duct_transmission(False)
        data["transmission"] = duct_transmission(True)
        print(f"duct sweeps done ({time.perf_counter() - t0:.0f}s)")

    savemat("validation/data/ngsolve_duct_band_reference.mat", data,
            do_compression=True)
    print("saved validation/data/ngsolve_duct_band_reference.mat")


if __name__ == "__main__":
    main()
