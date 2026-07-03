"""Export ngsolve.bem Laplace reference operators for one .vol fixture.

    python validation/exportNgsolveBemReference.py \
        fixtures/mesh_topology/unit_sphere_coarse.vol \
        validation/data/ngbem_reference_unit_sphere_coarse.mat

Assembles the Galerkin single layer V and double layer K on continuous P1
(NGSolve ``H1(mesh, order=1)``, whose dof numbering equals the .vol point
order) plus the boundary P1 surface mass M, densifies them, and saves a
MATLAB .mat consumed by ``verifyGalerkinAgainstNgsolve.m`` and
``testNgsolveBemCrossCheck.m``.

Convention check (measured on the sphere fixtures before this exporter was
committed): ngsolve.bem shares this repo's operator conventions exactly -
G = 1/(4*pi*r), K is the OUTWARD-normal double layer with the principal
value on the diagonal (K[1] = -1/2 to 1e-9, K[Y_1] -> -1/6), and interior
volume-vertex rows/columns are exactly zero, so restricting to the boundary
nodes via SurfaceMesh.volNodeIds is the whole reindexing.

ngsolve.bem integrates with Sauter-Schwab-type numerical singular
quadrature.  The export assembles at TWO integration orders (12 and 16) and
stores their relative difference in the .mat, so the artifact itself
records that the reference is quadrature-converged (measured ~2e-9 for V,
~2e-8 for K on the sphere fixtures).  The stored V/K are the intorder-16
matrices.

Requires NGSolve >= 6.2.2604 (ngsolve.bem) and scipy.  Python is needed
only to (re)generate the .mat; the MATLAB verifier reads the committed
artifact.
"""
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.io import savemat

import ngsolve
import ngsolve.bem as bem
from ngsolve import H1, BilinearForm, Mesh, TaskManager, ds


def export(vol_file, out_mat, intorder=16, check_intorder=12):
    with TaskManager():
        mesh = Mesh(vol_file)
        fes = H1(mesh, order=1)  # dof k (0-based) = .vol point k+1 (1-based)
        u, v = fes.TnT()
        mass = BilinearForm(fes, check_unused=False)
        mass += u * v * ds
        mass.Assemble()
        rows, cols, vals = mass.mat.COO()
        M = np.asarray(sp.csr_matrix(
            (vals, (rows, cols)), shape=(fes.ndof, fes.ndof)).todense())

        def assemble(io):
            t0 = time.perf_counter()
            V = bem.SingleLayerPotentialOperator(fes, fes, intorder=io)
            K = bem.DoubleLayerPotentialOperator(fes, fes, intorder=io)
            Vd = np.array(V.mat.ToDense(), copy=True)
            Kd = np.array(K.mat.ToDense(), copy=True)
            print(f"  intorder {io}: {time.perf_counter() - t0:.1f}s")
            return Vd, Kd

        print(f"{vol_file}: ndof {fes.ndof}")
        Vc, Kc = assemble(check_intorder)
        Vd, Kd = assemble(intorder)
        conv_v = np.linalg.norm(Vc - Vd) / np.linalg.norm(Vd)
        conv_k = np.linalg.norm(Kc - Kd) / np.linalg.norm(Kd)
        print(f"  intorder {check_intorder}-vs-{intorder}: "
              f"V {conv_v:.3e}  K {conv_k:.3e}")

        vtx = np.array([p.point for p in mesh.vertices])

    savemat(out_mat, {
        "V": Vd,
        "K": Kd,
        "M": M,
        "vtx": vtx,
        "intorder": float(intorder),
        "intorderConvergenceV": conv_v,
        "intorderConvergenceK": conv_k,
        "ngsolveVersion": ngsolve.__version__,
        "volFile": vol_file,
    }, do_compression=True)
    print(f"  saved {out_mat}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python exportNgsolveBemReference.py "
                 "<volFile> <outMat> [intorder]")
    export(sys.argv[1], sys.argv[2],
           int(sys.argv[3]) if len(sys.argv) > 3 else 16)
