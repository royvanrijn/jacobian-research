#!/usr/bin/env python3
"""Cluster recurring exact obstruction signatures into a compact report."""
import json
import sympy as sp
from collections import defaultdict
from pathlib import Path

root=Path(__file__).resolve().parents[1];results=root/"results"
pole=json.loads((results/"pole_kernels_deg4_box.json").read_text())["records"]
clusters=defaultdict(list)
v=sp.symbols("v")
for r in pole:
  degree=0 if r["p"]=="0" else sp.Poly(sp.sympify(r["p"]),v).degree()
  clusters[(degree,r["rank_Q"],tuple(r["earliest_forbidden"]))].append(r["p"])
stage=json.loads((results/"stage_a_support_families.json").read_text())
toric=json.loads((results/"stage_c_toric.json").read_text())
translated=json.loads((results/"stage_c_translated.json").read_text())
newton=json.loads((results/"stage_d_newton_bands.json").read_text())
core=json.loads((results/"triangle_d3_minimal_core.json").read_text())
regression=json.loads((results/"newton_9_27_regression.json").read_text())
reduced=json.loads((results/"reduced_72_108_systems_broad.json").read_text())
edgepowers=json.loads((results/"reduced_72_108_systems_edgepowers.json").read_text())
lines=["# Automatically clustered obstruction report","","## One-divisor pole matrices",""]
for (degree,rank,earliest),ps in sorted(clusters.items()):
  lines.append(f"- degree `{degree}`: rank `{rank}`, earliest pole `{earliest}`, representatives `{ps}`")
lines += ["","No tested prime changed these rational ranks.","","## Collision-normalized F4 families",""]
for r in stage:
  lines.append(f"- `{r['family']}` ({r['terms_per_coordinate']} terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.")
lines += ["","## Two-divisor toric families",""]
for r in toric:
  lines.append(f"- chart `{r['chart']}`, {r['terms_per_coordinate']} terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.")
lines += ["","## Translated two-divisor families",""]
for r in translated:
  lines.append(f"- `{r['family']}`: support {r['support_size']}, cancellation kernel {r['kernel_dimension']}, pole rank {r['pole_rank']}; exact `Q` basis `[1]` and three modular runs empty.")
lines += ["","## Newton-directed unequal-coordinate families",""]
for r in newton:
  dimensions=f"{r['P_kernel_dimension']}/{r['Q_kernel_dimension']}"
  if "timeout" in r["failure_class"]:
    outcome="recorded modular F4 timeout; no exact emptiness claim"
  else:
    outcome="exact collision-normalized `Q` basis `[1]`; three modular runs empty"
  lines.append(f"- `{r['family']}`: kernel dimensions {dimensions}; {outcome}.")
lines += ["","## Recurring reason","",
  "Pole cancellation alone has a large kernel and is not the obstruction. In every exactly completed sparse family, collision normalization turns the ideal into the unit ideal; all Stage D kernels still contain the inverse-chart Keller control. The dense total-degree Stage D control is explicitly unclassified after a scaling timeout."]
lines += ["","## Reduced degree-three contradiction core","",
  f"The original 23 equations reduce greedily to {core['core_size']}. The constant Jacobian coefficient and `Py(0)=0` are redundant because the retained derivative normalization fixes the constant term. Retained labels:","",
  "`"+"`, `".join(core["labels"])+"`.","",
  "The retained system is exactly empty over `Q`. It is greedily reduced, not claimed inclusion-minimal, because single-deletion ideals become costly positive-dimensional systems."]
lines += ["","## Published `(9,27)` regression","",
  f"The nine-equation Section 5 elimination reproduces equation (5.9) over `Q`: terminal relation found = `{regression['exact_terminal_relation_found']}`. The subsequent valuation split gives the published order contradiction `70/66/72` and degree contradiction `76/78/75`.",
  "This assumes the paper's prior reductions; it is not an independent re-proof of every cited proposition."]
lines += ["","## Reduced `(72,108)` target","",]
for r in reduced:
  run=r.get("run")
  outcome=(f"timeout `{run['timeout']}s` (complexity only)" if run and "timeout" in run else "not launched")
  lines.append(f"- broad `{r['case']}`: {r['variables']} variables, {r['jacobian_coefficient_equations']} equations; {outcome}.")
for r in edgepowers:
  run=r.get("run")
  outcome=(f"timeout `{run['timeout']}s` (complexity only)" if run and "timeout" in run else "not launched")
  lines.append(f"- experimental edge-power `{r['case']}`: {r['variables']} variables, {r['jacobian_coefficient_equations']} equations; {outcome}.")
lines += ["", "No `(72,108)` survivor or exclusion has been obtained. Original affine collision constraints are not valid in these reduced coordinates unless explicitly transported."]
path=results/"OBSTRUCTION_CLUSTERS.md";path.write_text("\n".join(lines)+"\n");print("wrote",path)
