"""Directed slab using compressed residual and a weighted annular error."""

import argparse
import hashlib
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import sys
import time

import numpy as np


LOWER = Q(89, 200)
UPPER = Q(447, 1000)
ORIGINAL_SHA256 = "c974ff60f8c1bf4723cd22d1a3074f932437b2303eb8219bcd29741a950aba06"
LOWER_FLOAT = math.nextafter(float(LOWER), -math.inf)
UPPER_FLOAT = math.nextafter(float(UPPER), math.inf)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_original():
    location = Path(__file__).resolve()
    original = location.with_name("check_three_cell_finite_amplitude.py")
    if not original.is_file():
        original = location.parents[1] / "research/plurality-stability" / original.name
    assert digest(original) == ORIGINAL_SHA256
    content = original.read_text()
    # Retarget angular integration and widen its trial-norm denominator.
    replacements = (
        ("rho = I(math.nextafter(.44, -math.inf), math.nextafter(.442, math.inf))",
         f"rho = I({LOWER_FLOAT!r}, {UPPER_FLOAT!r})"),
        ("norm = I.arb(arb(271)/1000)", "norm = I.arb(arb(274)/1000)"),
        ("ans = gamma_box(box)", "ans = gamma_box(box, Q(2, 125))"),
    )
    for old, new in replacements:
        assert content.count(old) == 1
        content = content.replace(old, new)
    namespace = {"__file__": str(original), "__name__": "_finite_amplitude_original"}
    exec(compile(content, str(original), "exec"), namespace)
    return namespace


def projected_assembly(original):
    """FA8-FA9, full-vector clipping AP3, and final quadratic determinant."""
    arb = original["arb"]
    pi = original["pi"]
    j1 = original["j1"]
    cs, ls = original["cs"], original["ls"]
    mb = original["M"]
    base = original["s"].base
    rat = lambda n, d=1: arb(n) / d
    rho = base.ball(LOWER_FLOAT, UPPER_FLOAT)
    d = 1-rho*rho
    alpha = rho/d.sqrt()
    kappa = rho*rho/d
    x0 = rat(1, 10)
    q = alpha*x0
    special = lambda z: pi.sqrt()*(1+2*z*z)*(z*z).exp()*z.erfc()-2*z
    weighted_j = j1*alpha*(-(1+kappa)*x0*x0).exp()*(special(q/2)+special(q)/2)
    slope = sum((c*l for c, l in zip(cs, ls) if c > 0), arb(0))
    pole = (sum(cs, arb(0))-slope*x0*x0)/rat(274, 1000).sqrt()
    constant = abs(mb)/rat(268, 1000).sqrt()
    A0, c, C, j = rat(1, 40), rat(373, 200), rat(14, 25), rat(107, 1000)
    gap, e, R2, P = rat(34, 25), rat(71, 4000), rat(1, 60000), rat(97, 20)
    H2 = ((R2.sqrt()+(R2+gap*e).sqrt())/gap)**2
    radial = rat(2, 3).sqrt()+A0*C
    X = A0*c/radial
    Phi = 2*j*c**3*A0/(3*radial)
    zeta = rat(1, 6)-rat(2, 3).sqrt()*A0*P
    annular_coefficient = rat(7, 625)
    cube_coefficient = rat(2, 125)
    theta = annular_coefficient/cube_coefficient
    rest = rat(16, 25)**2-(rat(89, 20)-rat(4023, 2000)**2)
    final = (annular_coefficient*rat(163, 100)-e)*(gap-rat(11, 10))-R2
    assert 2*cube_coefficient*theta/(1-theta) < rat(3, 40)
    assert weighted_j > j and pole > c and constant < C
    assert X < x0 and Phi > H2 and zeta > 0
    assert j1 > rat(9, 100) and rest > 0 and final > 0
    values = dict(weighted_J_local=weighted_j, pole_lower=pole,
                  pole_constant_upper=constant, integration_endpoint=X,
                  feasibility_Phi=Phi, feasibility_H_squared=H2,
                  feasibility_margin=Phi-H2, outside_cube_distance=zeta,
                  orthogonal_block_margin=rest, final_determinant=final)
    return {name: list(base.bounds(value)) for name, value in values.items()}


def weighted_angular_slack(original, table, exp_table, panels=512):
    """Uniform bound for E(x) times the annular row divided by J_+(x^2)."""
    I = original["I"]
    base = original["s"].base
    arb = original["arb"]
    edges = np.linspace(5/16, 2.0, panels+1)
    assert Q.from_float(float(edges[0])) == Q(5, 16)
    assert Q.from_float(float(edges[-1])) == 2
    xx = I(edges[:-1], edges[1:])
    widths = I(base.sub_lo(edges[1:], edges[:-1]), base.sub_hi(edges[1:], edges[:-1]))
    values = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = base.ball(float(lo), float(hi))
        value = sum((c*(-l*x*x).exp() for c, l in zip(original["cs"], original["ls"])), arb(0))
        value = (value/x-original["M"])/(arb(268)/1000).sqrt()
        values.append(base.bounds(abs(value))[1])
    norm = I(values)
    assert Q.from_float(max(values)) < Q(97, 20)
    one_sixth = I.arb(arb(1)/6)
    amplitude = I.arb(arb(1)/40)
    diameter = I.arb(arb(97)/20)+norm
    zeta = one_sixth-I.arb((arb(2)/3).sqrt())*amplitude*norm
    assert float(np.min(zeta.lo)) > 0
    error = I.arb(arb(7)/625)*(amplitude*diameter/zeta).square()
    error = I(np.maximum(error.lo, .075), np.maximum(error.hi, .075))

    rho = I(LOWER_FLOAT, UPPER_FLOAT)
    d = 1-rho.square()
    alpha = rho/I(base.sqrt_pos_lo(d.lo), base.sqrt_pos_hi(d.hi))
    q = alpha*xx
    jb = alpha*(table.query("HPLUS", q)+I.arb(original["j1"]/2)*table.query("S", q))
    assert float(np.min(jb.lo)) > 0
    # V(x^2,y^2)dnu(y^2)/(x exp(-kappa*x^2)) =
    # (4 rho/d^2)y^2 exp(-y^2/d-(2rho/d)xy)dy.
    xx_col = I(xx.lo[:,None], xx.hi[:,None])
    yy_row = I(xx.lo[None,:], xx.hi[None,:])
    exponent = yy_row.square()/d+(2*rho/d)*xx_col*yy_row
    expval = I(*exp_table.query("E", exponent.lo, exponent.hi))
    product = (4*rho/d.square())*yy_row.square()*expval*I(widths.lo[None,:], widths.hi[None,:])
    rows = np.array([base.sum_enclosure(values)[1] for values in product.hi])
    row_upper = (I(rows)/jb).hi
    weighted_upper = (I(row_upper)*error).hi
    assert Q.from_float(float(np.max(weighted_upper))) < Q(11, 20), float(np.max(weighted_upper))
    return {"panels": panels, "maximum_row_upper": float(np.max(row_upper)),
            "maximum_weighted_row_upper": float(np.max(weighted_upper))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert Q.from_float(LOWER_FLOAT) < LOWER < UPPER < Q.from_float(UPPER_FLOAT)
    start = time.monotonic()
    original = load_original()
    core = original["s"]
    base = core.base
    positive = core.positive
    constants = projected_assembly(original)
    print("Analytic assembly", constants, flush=True)

    print("Rebuilding validated one-variable tables", flush=True)
    table = core.InterpolatedTable(4., 8192)
    exp_table = base.build_exp_table(64., 131072, 2.)
    controls = original["numerical_controls"](table)
    weighted = weighted_angular_slack(original, table, exp_table)
    print("Weighted angular error", weighted, flush=True)
    gamma = original["gamma_cover"](limit=200000)
    assert gamma["status"] == "pass", gamma
    assert gamma["attempts"] == 2*(gamma["accepted"]+gamma["skipped"])-1
    print("Angular cube", gamma, flush=True)
    angular_band = original["angular_band"](exp_table)
    assert Q.from_float(angular_band["normalized_dirichlet_lower"]) > Q(163, 100)
    print("Angular band", angular_band, flush=True)

    mesh = original["Mesh"](16384)
    parents = []
    all_children = []
    edges = np.linspace(LOWER_FLOAT, UPPER_FLOAT, 6)
    for lo, hi in zip(edges[:-1], edges[1:]):
        lo, hi = float(lo), float(hi)
        parent = original["cover"](lo, hi, mesh, table, exp_table, 32)
        parent["rho"] = [lo, hi]
        assert original["interval_covers"](parent["slabs"], (Q.from_float(lo), Q.from_float(hi)))
        for child in parent["slabs"]:
            assert lo <= child["rho"][0] < child["rho"][1] <= hi
            lo_d, hi_d = map(Q.from_float, child["D"])
            assert Q(268, 1000) <= lo_d <= hi_d <= Q(274, 1000), child
            lo_a, hi_a = map(Q.from_float, child["a"])
            assert Q(4023, 2000) <= lo_a <= hi_a <= Q(8071, 4000), child
            compressed_bound = Q.from_float(child["r2"])-(lo_a-2)**2
            assert compressed_bound <= Q(1, 60000), (child, float(compressed_bound))
        parents.append(parent)
        all_children.extend(parent["slabs"])
        print("Spectral parent", {k: v for k, v in parent.items() if k != "slabs"}, flush=True)
    assert original["interval_covers"](all_children, (LOWER, UPPER))
    assert not original["interval_covers"](all_children[1:], (LOWER, UPPER))

    qt = positive.build_positive_q_table(4., 1024, 64, 2.)
    at = positive.build_positive_a_table(64., 4096, 64, 2.)
    upstream_constants = positive.build_constants()
    hs = []
    edges = np.linspace(LOWER_FLOAT, UPPER_FLOAT, 9)
    for lo, hi in zip(edges[:-1], edges[1:]):
        lo, hi = float(lo), float(hi)
        value = positive.certify_slab(lo, hi, 7., 6400, 48, at, qt, exp_table, upstream_constants)
        assert Q.from_float(value.total) <= Q(89, 20), value
        hs.append(dict(rho=[lo, hi], hs_squared=value.total, finite=value.finite,
                       tail=value.tail, kernel_tail=value.kernel_tail,
                       additive_tail=value.additive_tail))
        print("Hilbert--Schmidt parent", hs[-1], flush=True)
    assert original["interval_covers"](hs, (LOWER, UPPER))
    assert not original["interval_covers"](hs[1:], (LOWER, UPPER))
    controls.update(exact_rational_endpoint_cover=True, child_parent_containment=True,
                    incomplete_spectral_cover_rejected=True,
                    incomplete_hilbert_schmidt_cover_rejected=True,
                    all_spectral_caps=True, all_hilbert_schmidt_caps=True,
                    weighted_annular_error=True)
    sources = {src: digest(original["SOURCE_DIR"] / src)
               for src in ("plur_cert.py", "plur_cert_positive.py",
                           "check_positive_spectral_extension.py")}
    summary = dict(minimum_trial_norm=min(v["minD"] for v in parents),
                   maximum_trial_norm=max(v["maxD"] for v in parents),
                   minimum_rayleigh=min(v["mina"] for v in parents),
                   maximum_rayleigh=max(v["maxa"] for v in parents),
                   maximum_raw_residual_squared=max(v["maxr2"] for v in parents),
                   maximum_compressed_residual_squared=float(max(
                       Q.from_float(child["r2"])-(Q.from_float(child["a"][0])-2)**2
                       for child in all_children)),
                   maximum_hilbert_schmidt_squared=max(v["hs_squared"] for v in hs))
    result = dict(status="pass", exact_correlation_range=[str(LOWER), str(UPPER)],
                  original_checker_sha256=ORIGINAL_SHA256,
                  extension_checker_sha256=digest(Path(__file__)),
                  spectral_cover=parents, hilbert_schmidt_cover=hs,
                  angular_band=angular_band, angular_cube=gamma, weighted_angular_error=weighted,
                  analytic_assembly=constants, controls=controls,
                  input_sources=sources, summary=summary,
                  environment=dict(python=sys.version.split()[0],
                                   python_flint=__import__("flint").__version__,
                                   numpy=np.__version__),
                  elapsed_seconds=time.monotonic()-start)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status="pass", elapsed_seconds=result["elapsed_seconds"],
                          summary=summary)), flush=True)


if __name__ == "__main__":
    main()
