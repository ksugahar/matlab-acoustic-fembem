from __future__ import annotations

import argparse
import json
from pathlib import Path


def dense_to_list(mat):
    dense = mat.ToDense().NumPy()
    return [[float(dense[i, j]) for j in range(dense.shape[1])] for i in range(dense.shape[0])]


def summarize(vol_file: Path) -> dict[str, object]:
    out: dict[str, object] = {"ok": False}
    try:
        import ngsolve
        import ngsolve.bem as bem
        from ngsolve import BilinearForm, HCurl, H1, Mesh, curl, dx, grad

        mesh = Mesh(str(vol_file))

        h1 = H1(mesh, order=1)
        u, v = h1.TnT()
        k = BilinearForm(h1)
        k += grad(u) * grad(v) * dx
        k.Assemble()
        m = BilinearForm(h1)
        m += u * v * dx
        m.Assemble()

        hcurl = HCurl(mesh, order=0)
        a, b = hcurl.TnT()
        cm = BilinearForm(hcurl)
        cm += a * b * dx
        cm.Assemble()
        cc = BilinearForm(hcurl)
        cc += curl(a) * curl(b) * dx
        cc.Assemble()

        out.update(
            {
                "ok": True,
                "version": getattr(ngsolve, "__version__", ""),
                "meshVertices": int(mesh.nv),
                "meshElements": int(mesh.ne),
                "meshEdges": int(mesh.nedge),
                "h1Dofs": int(h1.ndof),
                "h1Stiffness": dense_to_list(k.mat),
                "h1Mass": dense_to_list(m.mat),
                "hcurlDofs": int(hcurl.ndof),
                "hcurlMass": dense_to_list(cm.mat),
                "hcurlCurlCurl": dense_to_list(cc.mat),
                "hasBem": True,
                "hasLaplaceSL": hasattr(bem, "LaplaceSL"),
                "hasHelmholtzSL": hasattr(bem, "HelmholtzSL"),
                "hasMaxwellSL": hasattr(bem, "MaxwellSL"),
            }
        )
    except Exception as exc:  # noqa: BLE001
        out.update({"ok": False, "errorType": type(exc).__name__, "errorMessage": str(exc)})
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol-file", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(summarize(Path(args.vol_file)), indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
