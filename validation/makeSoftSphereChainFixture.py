"""Generate the sonic-crystal chain fixture: 5 sound-soft spheres in a line.

    python validation/makeSoftSphereChainFixture.py

Writes fixtures/mesh_topology/soft_sphere_chain_5.vol: five spheres of
radius R = 0.3 with lattice constant d = 1.5 centered on the z axis at
z = -3, -1.5, 0, 1.5, 3 (netgen CSG, maxh = 0.15, ~120 boundary triangles
per sphere - the same faceting class as unit_sphere_coarse). The Bragg
wavenumber of the chain is k* = pi/d ~ 2.094; the finite-array insertion
loss shows the stop-band dip there (the sonic-crystal teaching rung,
COMSOL "Sonic Crystal" model class reproduced in the geometry this
first-order lane supports: sound-soft spheres, lossless).

Each sphere is its own CSG top-level object, so the .vol carries five
disjoint closed boundary surfaces plus the interior tetrahedra the
first-order parser expects.
"""
from netgen.csg import CSGeometry, Pnt, Sphere

RADIUS = 0.3
SPACING = 1.5
COUNT = 5
MAXH = 0.15
OUT = "fixtures/mesh_topology/soft_sphere_chain_5.vol"


def main():
    geo = CSGeometry()
    z0 = -SPACING * (COUNT - 1) / 2.0
    for i in range(COUNT):
        geo.Add(Sphere(Pnt(0, 0, z0 + i * SPACING), RADIUS).maxh(MAXH))
    mesh = geo.GenerateMesh(maxh=MAXH)
    mesh.Save(OUT)
    print(f"saved {OUT}: {len(mesh.Points())} points, "
          f"{len(mesh.Elements2D())} boundary triangles, "
          f"{len(mesh.Elements3D())} tetrahedra")


if __name__ == "__main__":
    main()
