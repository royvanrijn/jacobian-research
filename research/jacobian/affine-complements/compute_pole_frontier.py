"""Exact necessary linear differential-pole constraints at the new (3,4) frontier."""
from pathlib import Path
import hashlib,json,time
import compute_pole_layers as core

HERE=Path(__file__).resolve().parent
start=time.monotonic()
layer=core.layer(3,4)
kernel=layer.pop("kernel_basis")
layer["sparse_kernel_basis"]=[
    [[j,c] for j,c in enumerate(row) if c!="0"] for row in kernel]
report=dict(as_of="2026-11-14",outcome="PASS",layer=layer,
            scope="Necessary linear constraints only; rank-one coefficient tensors and finite maps remain unsolved",
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__),HERE/"compute_pole_layers.py"]},
            elapsed_seconds=round(time.monotonic()-start,6))
(HERE/"pole_frontier.json").write_text(json.dumps(report,indent=2)+"\n")
print("pole frontier",layer["pole_orders"],"sections",layer["section_dimensions"],
      "coefficient columns",len(layer["columns"]),"constraint rank",layer["constraint_rank"],
      "kernel dimension",layer["kernel_dimension"],"seconds",report["elapsed_seconds"])

