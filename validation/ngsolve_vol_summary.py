from __future__ import annotations

import argparse
import json
from pathlib import Path


def summarize(path: Path) -> dict[str, object]:
    try:
        from ngsolve import BND, Mesh

        mesh = Mesh(str(path))
        boundary_elements = list(mesh.Elements(BND))
        return {
            "source": str(path),
            "ok": True,
            "errorType": "",
            "errorMessage": "",
            "points": int(mesh.nv),
            "triangles": int(len(boundary_elements)),
            "tets": int(mesh.ne),
            "materials": int(len(mesh.GetMaterials())),
            "boundaries": int(len(mesh.GetBoundaries())),
            "hcurlEdges": int(mesh.nedge),
            "faces": int(mesh.nface),
            "materialNames": list(mesh.GetMaterials()),
            "boundaryNames": list(mesh.GetBoundaries()),
        }
    except Exception as exc:  # noqa: BLE001 - validation must report any loader failure.
        return {
            "source": str(path),
            "ok": False,
            "errorType": type(exc).__name__,
            "errorMessage": str(exc),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("vol_files", nargs="+")
    args = parser.parse_args()

    rows = [summarize(Path(vol_file)) for vol_file in args.vol_files]
    Path(args.output).write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
