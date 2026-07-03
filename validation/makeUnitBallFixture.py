"""Generate the finer unit-ball fixture for the coupled acoustic FEM/BEM rung.

    python validation/makeUnitBallFixture.py

Writes fixtures/mesh_topology/unit_ball_maxh018.vol: the unit ball meshed
at maxh = 0.18. The existing unit_sphere_fine (maxh 0.3, 18 interior
nodes) resolves the interior Helmholtz FEM at k1 ~ 2.9 only to the
(k1 h)^2 / 8 ~ 9% P1 interpolation class; this finer ball halves h so the
Anderson fluid-sphere gate can LOCK the mesh-convergence of the coupled
solve instead of just a loose band (the pollution/resolution teaching
moment made measurable).
"""
from netgen.csg import CSGeometry, Pnt, Sphere

MAXH = 0.18
OUT = "fixtures/mesh_topology/unit_ball_maxh018.vol"


def main():
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0).maxh(MAXH))
    mesh = geo.GenerateMesh(maxh=MAXH)
    mesh.Save(OUT)
    print(f"saved {OUT}: {len(mesh.Points())} points, "
          f"{len(mesh.Elements2D())} boundary triangles, "
          f"{len(mesh.Elements3D())} tetrahedra")


if __name__ == "__main__":
    main()
