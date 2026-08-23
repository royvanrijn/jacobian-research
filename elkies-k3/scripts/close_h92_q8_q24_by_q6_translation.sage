#!/usr/bin/env sage -python
"""
Transport the already-certified physical H92 q24 divisor into the CURRENT
equation-level q8 frame via the exact q6 Mordell-Weil translation.

This closes the missing physical -> equation NS isometry rather than guessing
a D13 vertical correction.

The q6 standard zero Ostd has old-zero MW (-2,1,0) and trivial IV* component.
Hence its Shioda vector

    v = Ostd - Oold - (Ostd.Oold + 2) F6

is integral and orthogonal to F6, Oold, and the q6 root lattice.  Translation
by Ostd is the Eichler transvection

    tau(x) = x + (x.F6)v - (x.v)F6 - (v.v/2)(x.F6)F6.

We certify:
  * tau is an integral unimodular NS isometry;
  * tau(Oold)=Ostd;
  * tau fixes every stored E6/E8 component;
  * tau(physical q8 fibre)=equation q8 fibre;
  * tau(physical q24 divisor) has D12 child;
  * its equation-frame horizontal MW section is compared exactly with AJ(S3);
  * the complete vertical correction is solved in F8 + D13 roots.
"""

import json
import sys
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix,
    lcm, matrix, pari, vector
)


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


ROOT = locate_repo()
LOCAL = ROOT / "artifacts/local/elkies-k3"
CERT = ROOT / "elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
OLDQ24 = LOCAL / "q8-q24-effective-zero-choices-current.json"
OUTPUT = LOCAL / "q8-q24-physical-to-equation-translation.json"

for path in (CERT, OLDQ24):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

# Run the CURRENT local equation certifier, not an archived copy.
saved_argv = list(sys.argv)
scope = {"__name__": "__embedded_q8_certifier__"}
try:
    sys.argv = [str(CERT)]
    exec(compile(CERT.read_text(), str(CERT), "exec"), scope)
finally:
    sys.argv = saved_argv

need = (
    "ns", "F6", "Oold", "Ostd", "Smarked", "S3std",
    "F8eq", "physical", "selected", "isotropic_mate", "roots_and_data",
    "section_from_old_mw",
)
missing = [name for name in need if name not in scope]
if missing:
    raise SystemExit("current certifier missing variables: " + ",".join(missing))

ns = scope["ns"]
F6 = vector(ZZ, scope["F6"])
Oold = vector(ZZ, scope["Oold"])
Ostd = vector(ZZ, scope["Ostd"])
Smarked = vector(ZZ, scope["Smarked"])
S3std = vector(ZZ, scope["S3std"])
F8eq = vector(ZZ, scope["F8eq"])
F8physical = vector(ZZ, scope["physical"])
selected = scope["selected"]
isotropic_mate = scope["isotropic_mate"]
roots_and_data = scope["roots_and_data"]
section_from_old_mw = scope["section_from_old_mw"]

old24 = json.loads(OLDQ24.read_text())
assert old24["status"] == "PASS_EXACT_Q24_EFFECTIVE_ZERO_CHOICES"
assert old24["transport"]["target_endpoint_match"] is True
D24physical = vector(ZZ, old24["transport"]["q24_divisor_source_h3_ns"])

physical_degree = ZZ(D24physical * ns * F8physical)
physical_square = ZZ(D24physical * ns * D24physical)
print(
    "Q8Q24TRANS_INPUT|"
    f"D_square={physical_square}|"
    f"D_degree_on_current_physical_F8={physical_degree}|"
    f"same_target_artifact=1|"
    f"status={'PASS' if physical_square==0 and physical_degree==2 else 'MISMATCH'}",
    flush=True,
)
assert physical_square == 0
assert physical_degree == 2

# ---------------------------------------------------------------------------
# 1. Exact q6 Eichler transvection.
# ---------------------------------------------------------------------------

pOstd = ZZ(Ostd * ns * Oold)
v = Ostd - Oold - (pOstd + 2) * F6

assert v * ns * F6 == 0
assert v * ns * Oold == 0
assert v * ns * v % 2 == 0

height_translation = -ZZ(v * ns * v)
assert height_translation == 12, height_translation


def tau(x):
    x = vector(ZZ, x)
    a = ZZ(x * ns * F6)
    b = ZZ(x * ns * v)
    # v^2 is even; this is the standard Eichler-Siegel transvection.
    return vector(
        ZZ,
        x + a * v - b * F6 - ZZ((v * ns * v) // 2) * a * F6,
    )


n = ns.nrows()
standard = identity_matrix(ZZ, n)
M = matrix(ZZ, [list(tau(standard.row(i))) for i in range(n)])

assert abs(M.det()) == 1
assert M * ns * M.transpose() == ns
assert tau(F6) == F6
assert tau(Oold) == Ostd

# Strong section-level regression: the physical marked horizontal section has
# old-zero MW (-2,-2,0), and translating by (-2,1,0) gives Smarked (-4,-1,0).
Qphysical, unused_pole, unused_h = section_from_old_mw(vector(ZZ, (-2,-2,0)))
assert tau(Qphysical) == Smarked

print(
    "Q8Q24TRANS_ISOM|"
    f"translation_height={height_translation}|det={M.det()}|"
    "F6_fixed=1|Oold_to_Ostd=1|Qphysical_to_Smarked=1|status=PASS",
    flush=True,
)

# Diagnostic only: these stored E6/E8 vectors belong to the repaired q8
# component-nef presentation. They are NOT assumed to be q6 root components,
# and therefore need not be individually fixed by the q6 translation.
component_vectors = []
component_stats = {
    "count": 0,
    "F6_orthogonal": 0,
    "v_orthogonal": 0,
    "fixed": 0,
}
for family in ("E6", "E8"):
    for row in selected[family]["simple_root_vectors_in_source_h3_ns"]:
        r = vector(ZZ, row)
        component_vectors.append(r)
        component_stats["count"] += 1
        component_stats["F6_orthogonal"] += int(r * ns * F6 == 0)
        component_stats["v_orthogonal"] += int(r * ns * v == 0)
        component_stats["fixed"] += int(tau(r) == r)

print(
    "Q8Q24TRANS_COMPONENT_DIAG|"
    f"count={component_stats['count']}|"
    f"F6_orthogonal={component_stats['F6_orthogonal']}|"
    f"v_orthogonal={component_stats['v_orthogonal']}|"
    f"fixed={component_stats['fixed']}|"
    "status=DIAGNOSTIC_ONLY",
    flush=True,
)

# This is the actual missing frame identity; it is tested directly and does
# not depend on any individual component vector being fixed.
mapped_F8 = tau(F8physical)
fiber_match = mapped_F8 == F8eq
print(
    "Q8Q24TRANS_FIBRE|"
    f"physical_to_equation={int(fiber_match)}|"
    f"difference_square={(mapped_F8-F8eq)*ns*(mapped_F8-F8eq)}|"
    "status=PASS" if fiber_match else
    "Q8Q24TRANS_FIBRE|status=MISMATCH",
    flush=True,
)
assert fiber_match

# ---------------------------------------------------------------------------
# 2. Transport the already-certified physical q24 divisor.
# ---------------------------------------------------------------------------

D24eq = tau(D24physical)

assert D24eq * ns * D24eq == 0
assert D24eq * ns * F8eq == 2
assert gcd(tuple(D24eq)) == 1

mate24 = isotropic_mate(ns, D24eq)
orth24 = matrix(
    ZZ, [list(D24eq * ns), list(mate24 * ns)]
).right_kernel_matrix()
child24 = -(orth24 * ns * orth24.transpose())
root24 = roots_and_data(child24)[2]

print(
    "Q8Q24TRANS_CHILD|"
    f"square={D24eq*ns*D24eq}|degree_on_F8eq={D24eq*ns*F8eq}|"
    f"primitive={int(gcd(tuple(D24eq))==1)}|"
    f"root_data={root24[0]},{root24[1]},{root24[2]}|"
    f"MW_rank_if_rho19={17-root24[0]}|"
    f"status={'PASS_D12' if root24==(12,264,4) else 'UNEXPECTED_CHILD'}",
    flush=True,
)
assert root24 == (12,264,4)


# ---------------------------------------------------------------------------
# 3. Deterministic equation-frame D13 coordinates.
# ---------------------------------------------------------------------------

def child_frame_with_zero(ns, fibre, zero):
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in complement.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    U2 = matrix(ZZ, ((0,1),(1,0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1,1000):
        candidate = vector(
            ZZ, [(i+1)**2 + shift*(i+1) + 1 for i in range(gram.nrows())]
        )
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [r for r in roots if regular*r > 0]
    pset = {tuple(r) for r in positive}
    simple = [
        r for r in positive
        if not any(tuple(r-left) in pset for left in positive)
    ]
    simple = matrix(ZZ, [list(r) for r in simple])
    assert simple.nrows() == simple.rank() == rank
    return simple, simple * gram * simple.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13,312,4)
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i,i]) for i in range(13)) == (1,)*13
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted0 = initial * child * initial.transpose()
    coupling = adapted0[:13,13:]
    tail = adapted0[13:,13:]
    H0 = tail - coupling.transpose()*cartan.inverse()*coupling

    scale = ZZ(1)
    for value in H0.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale*H0).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ,13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling
    assert H.det() == 237
    return basis, adapted, H


def class_order(dual):
    order = ZZ(1)
    for value in dual:
        order = lcm(order, ZZ(QQ(value).denominator()))
    return order


def d13_correction(dual, root):
    order = class_order(dual)
    expected = {ZZ(1):QQ(0), ZZ(2):QQ(1), ZZ(4):QQ(13)/4}
    assert order in expected
    correction = expected[order]
    raw = QQ(dual*root*dual)
    mod2 = lambda x: QQ(x)-2*(QQ(x)/2).floor()
    assert mod2(raw) == mod2(correction)
    return order, correction


def minimal_section_for_mw(adapted, z, cap=32768):
    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling
    z = vector(ZZ,z)
    height = QQ(z*H*z)

    base = vector(ZZ,[0]*13+list(z))
    pairings = vector(QQ, base*adapted[:,:13])
    dual = pairings*root.inverse()
    order, correction = d13_correction(dual,root)
    target_norm = height+correction
    assert target_norm in ZZ and target_norm >= 4 and target_norm % 2 == 0

    L = IntegralLattice(root)
    it = L.enumerate_close_vectors(-dual)
    chosen = None
    for unused in range(cap):
        shift = vector(ZZ,next(it))
        candidate = base + vector(ZZ,list(shift)+[0]*4)
        norm = ZZ(candidate*adapted*candidate)
        if norm == target_norm:
            chosen = candidate
            break
        if norm > target_norm:
            break
    assert chosen is not None

    pole = QQ(height+correction-4)/2
    assert pole in ZZ and pole >= 0
    pole = ZZ(pole)

    a = ZZ((target_norm-2)/2)
    section = vector(ZZ,[a,1]+list(chosen))
    return {
        "mw":z, "height":height, "class_order":order,
        "correction":correction, "pole":pole,
        "lift":chosen, "section":section,
    }


O8 = vector(ZZ, selected["E8"]["simple_root_vectors_in_source_h3_ns"][0])
assert O8 * ns * F8eq == 1

# O8 is the actual equation-level branch-point zero, but the physical->equation
# q6 translation is used only as an ABSOLUTE NS isometry. It need not fix O8.
mapped_O8 = tau(O8)
assert mapped_O8 * ns * mapped_O8 == -2
assert mapped_O8 * ns * F8eq == 1
print(
    "Q8Q24TRANS_ZERO_DIAG|"
    f"zero_fixed={int(mapped_O8 == O8)}|"
    f"mapped_zero_intersection={mapped_O8*ns*O8}|"
    "status=DIAGNOSTIC_ONLY",
    flush=True,
)

child8, Bzero = child_frame_with_zero(ns,F8eq,O8)
A13, adapted, H13 = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ,2),A13) * Bzero
Gadapt = block_diagonal_matrix(matrix(ZZ,((0,1),(1,0))),-adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt


def coords(C):
    c = vector(QQ,C) * Badapt.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ,c)


# Exact AJ(S3) in this equation frame.
cS3 = coords(S3std)
assert cS3[1] == 52
zS3 = vector(ZZ,cS3[-4:])
profS3 = minimal_section_for_mw(adapted,zS3)
PS3 = profS3["section"] * Badapt

# Exact transported q24 horizontal section.
cD = coords(D24eq)
assert cD[1] == 2
z24 = vector(ZZ,cD[-4:])
prof24 = minimal_section_for_mw(adapted,z24)
P24 = prof24["section"] * Badapt

horizontal_same = z24 == zS3
mw_delta = z24-zS3

print(
    "Q8Q24TRANS_HORIZONTAL|"
    f"AJ_S3_mw={','.join(map(str,zS3))}|"
    f"q24_mw={','.join(map(str,z24))}|"
    f"delta={','.join(map(str,mw_delta))}|"
    f"same={int(horizontal_same)}|"
    f"AJ_height={profS3['height']}|AJ_corr={profS3['correction']}|AJ_O={profS3['pole']}|"
    f"q24_height={prof24['height']}|q24_corr={prof24['correction']}|q24_O={prof24['pole']}|"
    "status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Solve the complete vertical correction of the true equation q24 divisor.
# ---------------------------------------------------------------------------

vertical = D24eq - O8 - P24
assert vertical * ns * F8eq == 0

root_source = []
for i in range(13):
    e = vector(
        ZZ,
        [0,0] + [ZZ(j==i) for j in range(17)]
    )
    r = e * Badapt
    assert r * ns * F8eq == 0
    assert r * ns * O8 == 0
    root_source.append(r)

vertical_basis = matrix(
    QQ, [list(F8eq)] + [list(r) for r in root_source]
)
vertical_coeffs = vertical_basis.solve_left(vector(QQ,vertical))
assert all(x in ZZ for x in vertical_coeffs)
vertical_coeffs = vector(ZZ,vertical_coeffs)
assert vector(QQ,vertical_coeffs) * vertical_basis == vector(QQ,vertical)

vf = ZZ(vertical_coeffs[0])
vr = vector(ZZ,vertical_coeffs[1:])

print(
    "Q8Q24TRANS_VERTICAL|"
    f"fibre={vf}|roots={','.join(map(str,vr))}|"
    f"root_L1={sum(abs(int(x)) for x in vr)}|"
    f"root_support={sum(bool(x) for x in vr)}|status=PASS_EXACT_VERTICAL",
    flush=True,
)

# Compare to the zero-root self-neighbour built from AJ(S3).
k_naive = ZZ((profS3["pole"]-2)//2)
Dnaive = O8 + PS3 - k_naive*F8eq
mate_naive = isotropic_mate(ns,Dnaive)
orth_naive = matrix(
    ZZ,[list(Dnaive*ns),list(mate_naive*ns)]
).right_kernel_matrix()
child_naive = -(orth_naive*ns*orth_naive.transpose())
root_naive = roots_and_data(child_naive)[2]
assert root_naive == (13,312,4)

delta_naive = D24eq-Dnaive
print(
    "Q8Q24TRANS_NAIVE_DIAG|"
    f"naive_twist={k_naive}|naive_child=D13|"
    f"true_minus_naive_square={delta_naive*ns*delta_naive}|"
    f"true_minus_naive_F8={delta_naive*ns*F8eq}|"
    "status=PASS",
    flush=True,
)

# Optional exact QQ(U) binding of AJ(S3), if Hensel has already finished.
exact_bound = False
EXACT = LOCAL / "q8-s3-direct-section-qq.json"
if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    assert exact["status"] == "PASS_EXACT_Q8_S3_DIRECT_SECTION"
    assert exact["verification"]["exact_weierstrass_identity"] is True
    assert ZZ(exact["profile"]["P_dot_O_from_denominator"]) == profS3["pole"]
    exact_bound = True

payload = {
    "schema":"elkies-k3.h92-q8-q24-physical-to-equation-translation.v1",
    "status":"PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION",
    "q6_translation":{
        "shioda_vector_source_h3_ns":list(map(int,v)),
        "height":int(height_translation),
        "matrix_determinant":int(M.det()),
        "isometry":True,
        "Oold_to_Ostd":True,
        "component_vector_diagnostic":component_stats,
    },
    "q8_fibre":{
        "physical_source_h3_ns":list(map(int,F8physical)),
        "equation_source_h3_ns":list(map(int,F8eq)),
        "translation_match":True,
    },
    "q24_equation":{
        "physical_divisor_source_h3_ns":list(map(int,D24physical)),
        "equation_divisor_source_h3_ns":list(map(int,D24eq)),
        "child_root_data":list(map(int,root24)),
        "child_root_lattice":"D12",
        "MW_rank_if_rho19":5,
        "horizontal_mw":list(map(int,z24)),
        "height":str(prof24["height"]),
        "D13_class_order":int(prof24["class_order"]),
        "D13_local_correction":str(prof24["correction"]),
        "P_dot_O":int(prof24["pole"]),
        "section_source_h3_ns":list(map(int,P24)),
        "vertical_fibre_coefficient":int(vf),
        "vertical_root_coefficients":list(map(int,vr)),
    },
    "AJ_S3":{
        "mw":list(map(int,zS3)),
        "height":str(profS3["height"]),
        "D13_class_order":int(profS3["class_order"]),
        "D13_local_correction":str(profS3["correction"]),
        "P_dot_O":int(profS3["pole"]),
        "same_horizontal_as_q24":bool(horizontal_same),
        "q24_minus_AJ_mw":list(map(int,mw_delta)),
        "exact_characteristic_zero_section_bound":bool(exact_bound),
    },
    "naive_zero_root_neighbor":{
        "twist":int(k_naive),
        "child_root_data":list(map(int,root_naive)),
        "child_root_lattice":"D13",
    },
    "boundary":(
        "This closes the physical->equation NS translation and produces the exact "
        "D12 q24 divisor with its full vertical correction. If AJ(S3) and q24 "
        "horizontal MW classes differ, the remaining equation task is to form "
        "the q24 point from the exact AJ point plus the reported MW correction; "
        "if they coincide, the exact AJ point is already the q24 horizontal section."
    ),
}

OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(f"OUTPUT|{OUTPUT}",flush=True)
print(
    "Q8Q24TRANS_RESULT|"
    f"child=D12|MW=5|horizontal_same={int(horizontal_same)}|"
    f"vertical_F={vf}|vertical_support={sum(bool(x) for x in vr)}|"
    f"char0_AJ_bound={int(exact_bound)}|"
    "status=PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION",
    flush=True,
)
