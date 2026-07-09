"""Export a convention-free curved-boundary-node companion for a curved mesh.

This is the netgen/NGSolve side of the MATLAB curved-panel BEM's optional
netgen path.  A Netgen curved .vol stores its curving as a `curvedelements`
COEFFICIENT basis (integrated-Legendre / Jacobi) that is fragile to reparse in
MATLAB.  Instead, NGSolve -- which owns the curved-element engine -- evaluates
each boundary triangle's geometry map at the Lagrange reference nodes and writes
the physical curved node COORDINATES (convention free) to a JSON companion.  The
saved .vol stays P1 (straight); the curving lives entirely in the companion, so
MATLAB never touches the coefficient basis.  MATLAB reads both with
curvedQuadratureFromNetgen and builds a CurvedPanelQuadrature.

Usage:
    python export_curved_boundary_nodes.py --make-sphere --radius 1.0 --maxh 0.8 \
        --curve-order 2 --vol-out ng_sphere.vol --json-out ng_sphere_curved_p2.json
    python export_curved_boundary_nodes.py --vol existing_with_geometry.vol ...

Only curve order 2 (6-node quadratic panels) is emitted; that is the primary
case the MATLAB reader consumes.  Node contract per triangle:
    corners[3]    : the 3 P1 corner coordinates (on the surface)
    edgemids[3]   : curved midpoint of edge e, between corners[e] and
                    corners[(e+1) % 3]  (positional pairing, convention free)
"""
import argparse
import json

import numpy as np
from ngsolve import Mesh, BND, IntegrationRule, TaskManager


def _eval_ref(trafo, u, v):
    for ip in IntegrationRule([(u, v)], [1.0]):
        mip = trafo(ip)
        return [float(mip.point[0]), float(mip.point[1]), float(mip.point[2])]


def export_curved_nodes(mesh, curve_order):
    if curve_order != 2:
        raise ValueError("only curve order 2 is exported (got %d)" % curve_order)
    corner_ref = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    edge_ref = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]   # edge e joins corner e and (e+1)%3
    triangles = []
    for el in mesh.Elements(BND):
        trafo = mesh.GetTrafo(el)
        corners = [_eval_ref(trafo, u, v) for u, v in corner_ref]
        edgemids = [_eval_ref(trafo, u, v) for u, v in edge_ref]
        triangles.append({"corners": corners, "edgemids": edgemids})
    return triangles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vol", help="input .vol carrying geometry (Curve projects onto it)")
    ap.add_argument("--make-sphere", action="store_true", help="generate an OCC sphere instead")
    ap.add_argument("--radius", type=float, default=1.0)
    ap.add_argument("--maxh", type=float, default=0.8)
    ap.add_argument("--curve-order", type=int, default=2)
    ap.add_argument("--vol-out", required=True, help="P1 .vol to write (read by MATLAB VolMesh)")
    ap.add_argument("--json-out", required=True, help="curved-node companion JSON to write")
    args = ap.parse_args()

    with TaskManager():
        if args.make_sphere:
            from netgen.occ import Sphere, Pnt, OCCGeometry
            geo = OCCGeometry(Sphere(Pnt(0, 0, 0), args.radius))
            ngm = geo.GenerateMesh(maxh=args.maxh)
        elif args.vol:
            ngm = Mesh(args.vol).ngmesh
        else:
            raise SystemExit("pass --make-sphere or --vol")
        ngm.Save(args.vol_out)               # P1 .vol (straight); curving is in the JSON
        mesh = Mesh(ngm)
        mesh.Curve(args.curve_order)
        triangles = export_curved_nodes(mesh, args.curve_order)

    # sanity: curved edge mids should sit ~on the sphere if this is a sphere
    if args.make_sphere:
        mids = np.array([m for t in triangles for m in t["edgemids"]])
        dev = float(np.max(np.abs(np.linalg.norm(mids, axis=1) - args.radius)))
        print("max curved edge-mid radius deviation: %.2e" % dev)

    out = {
        "kind": "curved_boundary_nodes",
        "source": "make-sphere" if args.make_sphere else args.vol,
        "curve_order": args.curve_order,
        "radius": args.radius if args.make_sphere else None,
        "n_triangles": len(triangles),
        "triangles": triangles,
    }
    with open(args.json_out, "w") as f:
        json.dump(out, f)
    print("wrote %s (%d triangles) + %s" % (args.json_out, len(triangles), args.vol_out))


if __name__ == "__main__":
    main()
