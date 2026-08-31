#!/usr/bin/env sage -python
"""Enumerate and mark the q12/o5867 P.O=0 compiler shell modulo p.

This is a bounded polynomial-square enumeration on the exact P1229-pointed
q8/o376 child.  A tiny generated C helper performs only the finite evaluation
loop; Sage verifies every returned section, applies the embedding-15 I2
profile and inverse pointed-quartic parent-degree gates, computes the ordinary
13-by-12 coefficient Jacobian, and tests exact lattice-derived group-word
fingerprints.  No Groebner basis or elimination is used.

The default prime is 67.  At that prime the q4/o164 equation-basis scaling is
zero, so direct evaluation of the parent B-basis words is deliberately not
claimed; the exact QQ Q3 section on the q8 child remains a valid word anchor.
"""

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
FRONTIER = GEN / "elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
INPUTS = (Q8, Q3, FRONTIER, MODEL, BASIS, C8, HORIZONTAL)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=79)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--direct-words-only", action="store_true",
    help="evaluate the supplied q4 MW words and stop before any shell",
)
parser.add_argument(
    "--preflight-only", action="store_true",
    help="record the parent/q8 good-prime gates and stop",
)
parser.add_argument(
    "--include-all-records", action="store_true",
    help="retain the complete signed polynomial shell for downstream diagnostics",
)
args = parser.parse_args()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
OUTPUT = (
    args.output.resolve() if args.output else
    LOCAL / f"q12o5867-p0-shell-word-fingerprints-mod{prime}.json"
)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
q3_exact = json.loads(Q3.read_text())
frontier = json.loads(FRONTIER.read_text())
model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
c8 = json.loads(C8.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert q3_exact["status"] == "PASS_EXACT_QQ_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD"
assert frontier["status"] == "PASS_EXACT_ROOTLESS_P0_SECTION_WORD_FRONTIER"

F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()
RX = PolynomialRing(F, "x")
xvar = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % prime == 0:
        raise ZeroDivisionError(f"bad reduction denominator at p={prime}")
    return F(value.numerator()) / F(denominator)


def record_function(record):
    numerator = R([reduce_qq(v) for v in record["numerator_coefficients_low_to_high"]])
    denominator = R([reduce_qq(v) for v in record["denominator_coefficients_low_to_high"]])
    return K(numerator) / K(denominator)


child = q8["child"]
A = R([reduce_qq(v) for v in child["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(v) for v in child["minimal_B_coefficients_low_to_high"]])
assert (A.degree(), B.degree()) == (8, 12)

# Equation order is the stored finite-I2 order followed by infinity.  In the
# selected embedding 15, marked 0100 is equation profile 1000.
supports = []
nodes = []
for fibre in child["finite_reducible_fibres"]:
    factor_text = fibre["factor"]
    support = F.zero() if factor_text == "u" else reduce_qq(factor_text.split(" - ", 1)[1])
    cubic = xvar**3 + A(support)*xvar + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    supports.append(support)
    nodes.append(-repeated[0]/repeated[1])
assert len(supports) == 3 and len(set(supports)) == 3
infinity_cubic = xvar**3 + A[8]*xvar + B[12]
infinity_repeated = infinity_cubic.gcd(infinity_cubic.derivative())
assert infinity_repeated.degree() == 1
infinity_node = -infinity_repeated[0]/infinity_repeated[1]


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def coefficient_jacobian(X, Y):
    dx = -3*X**2-A
    dy = 2*Y
    columns = [padded(u**power*dx, 13) for power in range(5)]
    columns += [padded(u**power*dy, 13) for power in range(7)]
    return matrix(F, columns).transpose()


def component_profile(X, Y):
    result = []
    for support, node in zip(supports, nodes):
        result.append(int(X(support) == node and Y(support) == 0))
    infinity_x = X[4] if X.degree() == 4 else F.zero()
    infinity_y = Y[6] if Y.degree() == 6 else F.zero()
    result.append(int(infinity_x == infinity_node and infinity_y == 0))
    return result


pointing = q8["preferred_pointed_zero"]
a1, a2, a3, unused_a4, unused_a6 = [
    record_function(record) for record in pointing["generalized_weierstrass_a_invariants"]
]
ordinate = record_function(pointing["quartic_ordinate"])
b2 = a1**2+4*a2
old_base_support = reduce_qq(pointing["old_base_coordinate"])


def inverse_parent_base(X, Y):
    x_general = K(X)/9-b2/12
    y_general = K(Y)/27-(a1*x_general+a3)/2
    if not y_general:
        return None
    return K(old_base_support+2*ordinate*(x_general+a2)/y_general)


def rational_degree(value):
    return max(value.numerator().degree(), value.denominator().degree())


q3_section = q3_exact["section"]
Q3_point = (
    R([reduce_qq(v) for v in q3_section["x_coefficients_low_to_high"]]),
    R([reduce_qq(v) for v in q3_section["y_coefficients_low_to_high"]]),
)
assert Q3_point[1]**2 == Q3_point[0]**3+A*Q3_point[0]+B
q3_parent = inverse_parent_base(*Q3_point)
q3_gate = {
    "equation_component_profile": component_profile(*Q3_point),
    "inverse_parent_degree": int(rational_degree(q3_parent)),
    "ordinary_coefficient_jacobian_rank": int(coefficient_jacobian(*Q3_point).rank()),
}


def rational_gate(label, value):
    value = QQ(value)
    return {
        "label": label,
        "numerator_mod_p": int(ZZ(value.numerator()) % prime),
        "denominator_mod_p": int(ZZ(value.denominator()) % prime),
        "p_adic_valuation": int(value.valuation(prime)),
        "is_unit": bool(value.valuation(prime) == 0),
    }


coordinate_change = model["exact_coordinate_change"]
coordinate_gates = [
    rational_gate("base_scale_c", coordinate_change["c"]),
    rational_gate("xy_scale_s", coordinate_change["s"]),
]


def bad_denominators(label, values):
    answer = []
    for index, value in enumerate(values):
        value = QQ(value)
        if value.denominator() % prime == 0:
            answer.append({
                "label": label,
                "coefficient_index": index,
                "denominator_p_adic_valuation": int(ZZ(value.denominator()).valuation(prime)),
            })
    return answer


parent_denominator_failures = []
for index, section in enumerate(basis["resolved_hensel"]["sections"]):
    parent_denominator_failures += bad_denominators(
        f"B{index}.x", section["x_coefficients_low_to_high"]
    )
    parent_denominator_failures += bad_denominators(
        f"B{index}.y", section["y_coefficients_low_to_high"]
    )
for coordinate in ("x", "y"):
    parent_denominator_failures += bad_denominators(
        f"H.{coordinate}.numerator",
        horizontal["section"][coordinate]["numerator_coefficients_low_to_high"],
    )
    parent_denominator_failures += bad_denominators(
        f"H.{coordinate}.denominator",
        horizontal["section"][coordinate]["denominator_coefficients_low_to_high"],
    )
    parent_denominator_failures += bad_denominators(
        f"C8opposite.{coordinate}.numerator",
        c8["opposite_constant_support_section"][coordinate]["numerator_coefficients_low_to_high"],
    )
    parent_denominator_failures += bad_denominators(
        f"C8opposite.{coordinate}.denominator",
        c8["opposite_constant_support_section"][coordinate]["denominator_coefficients_low_to_high"],
    )

preflight_pass = (
    all(item["is_unit"] for item in coordinate_gates)
    and not parent_denominator_failures
    and q3_gate["equation_component_profile"] == [0, 0, 0, 0]
    and q3_gate["inverse_parent_degree"] == 1
    and q3_gate["ordinary_coefficient_jacobian_rank"] == 12
)
if not preflight_pass:
    preflight_payload = {
        "schema": "elkies-k3.h92-q12o5867-direct-word-preflight-modp.v1",
        "status": "REJECTED_MODP_Q12O5867_DIRECT_WORD_PREFLIGHT",
        "prime": int(prime),
        "preflight": {
            "pass": False,
            "q4_compact_coordinate_change": coordinate_gates,
            "parent_equation_coordinate_bad_denominators": parent_denominator_failures,
            "q8_Q3_gate": q3_gate,
        },
        "direct_words": {
            "Q1": "H-B3",
            "Q2": "H+B0+B2-B3+B4",
            "Q4": "H-C8opp-B1-B2+B3+B5-B6-B7",
            "evaluated": False,
        },
        "shell": {"run": False, "reason": "preflight rejection; bounded cross-check not needed"},
        "method": {
            "large_Groebner_required": False,
            "elimination_required": False,
            "QQ_lift_attempted": False,
            "runtime_seconds": time.monotonic()-started,
        },
        "proof_boundary": (
            "This artifact is a good-prime rejection only. The direct marked words were not "
            "evaluated and no shell or characteristic-zero lift was run."
        ),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
        },
    }
    OUTPUT.write_text(json.dumps(preflight_payload, indent=2, sort_keys=True)+"\n")
    print(
        "Q12O5867DIRECTPREFLIGHT|prime={}|coordinate_units={}|bad_denominators={}|"
        "Q3_profile={}|Q3_degree={}|Q3_rank={}|status={}|output={}".format(
            prime, [item["is_unit"] for item in coordinate_gates],
            len(parent_denominator_failures), q3_gate["equation_component_profile"],
            q3_gate["inverse_parent_degree"], q3_gate["ordinary_coefficient_jacobian_rank"],
            preflight_payload["status"], OUTPUT,
        ), flush=True,
    )
    sys.exit(0)


if args.preflight_only:
    preflight_payload = {
        "schema": "elkies-k3.h92-q12o5867-direct-word-preflight-modp.v1",
        "status": "PASS_MODP_Q12O5867_DIRECT_WORD_PREFLIGHT",
        "prime": int(prime),
        "preflight": {
            "pass": True,
            "q4_compact_coordinate_change": coordinate_gates,
            "parent_equation_coordinate_bad_denominators": [],
            "q8_Q3_gate": q3_gate,
        },
        "direct_words": {"evaluated": False},
        "shell": {"run": False, "reason": "preflight-only was requested"},
        "method": {
            "large_Groebner_required": False,
            "elimination_required": False,
            "QQ_lift_attempted": False,
            "runtime_seconds": time.monotonic()-started,
        },
        "proof_boundary": (
            "This is a good-prime preflight only. No direct word, shell, or "
            "characteristic-zero lift was run."
        ),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
        },
    }
    OUTPUT.write_text(json.dumps(preflight_payload, indent=2, sort_keys=True)+"\n")
    print(
        "Q12O5867DIRECTPREFLIGHT|prime={}|coordinate_units=[True, True]|"
        "bad_denominators=0|Q3_profile={}|Q3_degree={}|Q3_rank={}|status={}|output={}".format(
            prime, q3_gate["equation_component_profile"], q3_gate["inverse_parent_degree"],
            q3_gate["ordinary_coefficient_jacobian_rank"], preflight_payload["status"], OUTPUT,
        ), flush=True,
    )
    sys.exit(0)


if args.direct_words_only:
    # The supplied words live in the q4 Mordell--Weil group.  Test whether
    # their actual q4 section curves have the required degrees over the q8
    # pencil before attempting the degree-one pointed-quartic inversion.
    PT = PolynomialRing(F, "T")
    parent_T = PT.gen()
    PK = PT.fraction_field()

    def parent_polynomial(values):
        return PT([reduce_qq(value) for value in values])

    parent_A = parent_polynomial(model["compact_model"]["A_coefficients_low_to_high"])
    parent_B = parent_polynomial(model["compact_model"]["B_coefficients_low_to_high"])

    def parent_record_function(record):
        return PK(parent_polynomial(record["numerator_coefficients_low_to_high"])) / PK(
            parent_polynomial(record["denominator_coefficients_low_to_high"])
        )

    def checked_parent_point(x_coordinate, y_coordinate):
        point = PK(x_coordinate), PK(y_coordinate)
        assert point[1]**2 == point[0]**3+PK(parent_A)*point[0]+PK(parent_B)
        return point

    def parent_neg(point):
        return None if point is None else (point[0], -point[1])

    def parent_add(left, right):
        if left is None:
            return right
        if right is None:
            return left
        x1, y1 = left
        x2, y2 = right
        if x1 == x2:
            if y1 == -y2:
                return None
            slope = (3*x1**2+PK(parent_A))/(2*y1)
        else:
            slope = (y2-y1)/(x2-x1)
        x3 = slope**2-x1-x2
        return checked_parent_point(x3, slope*(x1-x3)-y1)

    parent_basis = [
        checked_parent_point(
            parent_polynomial(section["x_coefficients_low_to_high"]),
            parent_polynomial(section["y_coefficients_low_to_high"]),
        )
        for section in basis["resolved_hensel"]["sections"]
    ]
    parent_horizontal = checked_parent_point(
        parent_record_function(horizontal["section"]["x"]),
        parent_record_function(horizontal["section"]["y"]),
    )
    c8_record = c8["opposite_constant_support_section"]
    c8_old = (
        parent_record_function(c8_record["x"]),
        parent_record_function(c8_record["y"]),
    )
    base_scale = reduce_qq(coordinate_change["c"])
    xy_scale = reduce_qq(coordinate_change["s"])

    def substitute_parent_scaled(function):
        numerator = PT(function.numerator()(base_scale*parent_T))
        denominator = PT(function.denominator()(base_scale*parent_T))
        return PK(numerator)/PK(denominator)

    c8_opposite = checked_parent_point(
        substitute_parent_scaled(c8_old[0])/xy_scale**2,
        substitute_parent_scaled(c8_old[1])/xy_scale**3,
    )
    named_words = {
        "Q1": [(1, parent_horizontal), (-1, parent_basis[3])],
        "Q2": [
            (1, parent_horizontal), (1, parent_basis[0]), (1, parent_basis[2]),
            (-1, parent_basis[3]), (1, parent_basis[4]),
        ],
        "Q4": [
            (1, parent_horizontal), (-1, c8_opposite), (-1, parent_basis[1]),
            (-1, parent_basis[2]), (1, parent_basis[3]), (1, parent_basis[5]),
            (-1, parent_basis[6]), (-1, parent_basis[7]),
        ],
    }
    expected_degrees = {"Q1": 3, "Q2": 2, "Q4": 2}
    parent_points = {}
    for name, word in named_words.items():
        point = None
        for coefficient, addend in word:
            point = parent_add(point, addend if coefficient == 1 else parent_neg(addend))
        parent_points[name] = point

    horizontal_x, horizontal_y = parent_horizontal
    denominator_x = parent_polynomial(
        horizontal["section"]["x"]["denominator_coefficients_low_to_high"]
    )
    denominator_y = parent_polynomial(
        horizontal["section"]["y"]["denominator_coefficients_low_to_high"]
    )
    Z = PT.one()
    for factor, exponent in denominator_x.factor():
        assert int(exponent) % 2 == 0
        Z *= factor.monic()**(int(exponent)//2)
    assert Z**2 == denominator_x and Z**3 == denominator_y
    resolved_pairs = [
        (
            parent_polynomial(item["AA_coefficients_low_to_high"]),
            parent_polynomial(item["BB_coefficients_low_to_high"]),
        )
        for item in q8["resolved_RR"]["resolved_basis_pairs"]
    ]
    direct_records = {}
    for name, (parent_x, parent_y) in parent_points.items():
        slope = (parent_y+horizontal_y)/(parent_x-horizontal_x)
        restrictions = [
            PK(AA)+PK(BB)*PK(Z)*slope for AA, BB in resolved_pairs
        ]
        new_base = -restrictions[0]/restrictions[1]
        restriction_degree = int(max(
            new_base.numerator().degree(), new_base.denominator().degree()
        ))
        direct_records[name] = {
            "word": {
                "Q1": "H-B3",
                "Q2": "H+B0+B2-B3+B4",
                "Q4": "H-C8opp-B1-B2+B3+B5-B6-B7",
            }[name],
            "parent_MW_section_x_degrees_numerator_denominator": [
                int(parent_x.numerator().degree()), int(parent_x.denominator().degree())
            ],
            "parent_MW_section_y_degrees_numerator_denominator": [
                int(parent_y.numerator().degree()), int(parent_y.denominator().degree())
            ],
            "q8_pencil_restriction_degree": restriction_degree,
            "expected_physical_multisection_degree": expected_degrees[name],
            "degree_matches_physical_class": restriction_degree == expected_degrees[name],
            "exact_parent_weierstrass_identity_mod_p": True,
        }
    degrees_match = all(
        record["degree_matches_physical_class"] for record in direct_records.values()
    )
    direct_payload = {
        "schema": "elkies-k3.h92-q12o5867-direct-q4-mw-word-audit-modp.v1",
        "status": (
            "PASS_MODP_Q12O5867_DIRECT_WORDS_READY_FOR_POINTED_TRANSPORT"
            if degrees_match else
            "REJECTED_MODP_Q12O5867_MW_WORDS_ARE_NOT_PHYSICAL_MULTISECTIONS"
        ),
        "prime": int(prime),
        "preflight": {
            "pass": True,
            "q4_compact_coordinate_change": coordinate_gates,
            "parent_equation_coordinate_bad_denominators": parent_denominator_failures,
            "q8_Q3_gate": q3_gate,
        },
        "direct_words": direct_records,
        "pointed_q8_transport_run": False,
        "pointed_q8_transport_reason": (
            None if degrees_match else
            "the q4 MW section representatives have different q8-pencil degrees from the physical compiler multisection classes; degree-one Mobius inversion would be invalid"
        ),
        "shell": {"run": False, "reason": "direct-words-only was requested"},
        "method": {
            "large_Groebner_required": False,
            "elimination_required": False,
            "QQ_lift_attempted": False,
            "runtime_seconds": time.monotonic()-started,
        },
        "proof_boundary": (
            "This tests the supplied exact q4 Mordell--Weil words as actual section curves. "
            "A marked MW tail is the Abel/Jacobian representative of a physical multisection "
            "class and does not retain its full U/root data. A failed degree comparison blocks "
            "the degree-one pointed-quartic transport and constructs no q8-child seed."
        ),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
        },
    }
    OUTPUT.write_text(json.dumps(direct_payload, indent=2, sort_keys=True)+"\n")
    print(
        "Q12O5867DIRECTWORDS|prime={}|degrees={}|expected={}|status={}|output={}".format(
            prime,
            {name: row["q8_pencil_restriction_degree"] for name, row in direct_records.items()},
            expected_degrees, direct_payload["status"], OUTPUT,
        ), flush=True,
    )
    sys.exit(0)

# Add one ordinary evaluation point; the leading coefficient is the fifth
# interpolation functional.
generic = next(F(i) for i in range(int(prime)) if F(i) not in supports)
evaluation_points = supports + [generic]
interpolation = matrix(F, [
    [point**degree for degree in range(5)]
    for point in evaluation_points
] + [[0, 0, 0, 0, 1]]).inverse()


def c_array(values):
    return ",".join(str(int(v)) for v in values)


C_SOURCE = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define P @P@
static int A[9]={@A@}, B[13]={@B@};
static int pts[4]={@PTS@}, nodes[4]={@NODES@};
static int invm[5][5]={@INVM@};
static int md(long long x){ x%=P; if(x<0)x+=P; return (int)x; }
static int pw(int a,int n){int r=1;while(n){if(n&1)r=md((long long)r*a);a=md((long long)a*a);n>>=1;}return r;}
static int val(const int *f,int d,int z){int r=0;for(int i=d;i>=0;i--)r=md((long long)r*z+f[i]);return r;}
static int square(int a){return a==0 || pw(a,(P-1)/2)==1;}
static int sqrt0(int a){for(int y=0;y<P;y++)if(md((long long)y*y)==a)return y;return -1;}
static int rootsquare(const int *r,int *y){
 int d=12;while(d>=0 && r[d]==0)d--; if(d<0){for(int i=0;i<7;i++)y[i]=0;return 1;}
 if(d&1)return 0; int m=d/2,s=sqrt0(r[d]);if(s<0)return 0;
 for(int i=0;i<7;i++)y[i]=0;y[m]=s;int den=md(2*s),iden=pw(den,P-2);
 for(int k=d-1;k>=m;k--){int j=k-m;long long known=0;for(int a=j+1;a<=m;a++){int b=k-a;if(b>=0&&b<=m)known+=(long long)y[a]*y[b];}y[j]=md((long long)(r[k]-md(known))*iden);}
 for(int k=0;k<=12;k++){long long q=0;for(int i=0;i<=6;i++){int j=k-i;if(j>=0&&j<=6)q+=(long long)y[i]*y[j];}if(md(q)!=r[k])return 0;}return 1;
}
static void emit(int mode,int *x,int *y){printf("M%d|X=",mode);for(int i=0;i<5;i++)printf("%s%d",i?",":"",x[i]);printf("|Y=");for(int i=0;i<7;i++)printf("%s%d",i?",":"",y[i]);puts("");}
int main(void){
 int allowed[2][5][P],n[2][5]={0};
 for(int mode=0;mode<2;mode++)for(int j=0;j<5;j++)for(int x=0;x<P;x++){
  int rhs;if(j<4)rhs=md((long long)x*x%P*x+(long long)val(A,8,pts[j])*x+val(B,12,pts[j]));else rhs=md((long long)x*x%P*x+(long long)A[8]*x+B[12]);
  if(!square(rhs))continue;
  if(j<3 || j==4){int profile=(mode==1 && j==0);int node=(j==4?nodes[3]:nodes[j]);if((x==node)!=profile)continue;}
  allowed[mode][j][n[mode][j]++]=x;
 }
 for(int mode=0;mode<2;mode++){
  fprintf(stderr,"MODE%d|lists=%d,%d,%d,%d,%d\n",mode,n[mode][0],n[mode][1],n[mode][2],n[mode][3],n[mode][4]);
  for(int i0=0;i0<n[mode][0];i0++)for(int i1=0;i1<n[mode][1];i1++)for(int i2=0;i2<n[mode][2];i2++)for(int i3=0;i3<n[mode][3];i3++)for(int i4=0;i4<n[mode][4];i4++){
   int vv[5]={allowed[mode][0][i0],allowed[mode][1][i1],allowed[mode][2][i2],allowed[mode][3][i3],allowed[mode][4][i4]},x[5],x2[9]={0},x3[13]={0},rhs[13],y[7];
   for(int a=0;a<5;a++){long long s=0;for(int b=0;b<5;b++)s+=(long long)invm[a][b]*vv[b];x[a]=md(s);}
   for(int i=0;i<5;i++)for(int j=0;j<5;j++)x2[i+j]=md(x2[i+j]+(long long)x[i]*x[j]);
   for(int i=0;i<9;i++)for(int j=0;j<5;j++)x3[i+j]=md(x3[i+j]+(long long)x2[i]*x[j]);
   for(int k=0;k<13;k++){long long s=x3[k]+B[k];for(int i=0;i<9;i++){int j=k-i;if(j>=0&&j<5)s+=(long long)A[i]*x[j];}rhs[k]=md(s);}
   if(!rootsquare(rhs,y))continue;emit(mode,x,y);int nonzero=0;for(int i=0;i<7;i++)nonzero|=y[i];if(nonzero){for(int i=0;i<7;i++)y[i]=md(-y[i]);emit(mode,x,y);}
  }
 }
 return 0;
}
'''
C_SOURCE = C_SOURCE.replace("@P@", str(int(prime)))
C_SOURCE = C_SOURCE.replace("@A@", c_array([A[i] for i in range(9)]))
C_SOURCE = C_SOURCE.replace("@B@", c_array([B[i] for i in range(13)]))
C_SOURCE = C_SOURCE.replace("@PTS@", c_array(evaluation_points))
C_SOURCE = C_SOURCE.replace("@NODES@", c_array(nodes + [infinity_node]))
C_SOURCE = C_SOURCE.replace("@INVM@", ",".join(
    "{" + c_array(row) + "}" for row in interpolation.rows()
))

with tempfile.TemporaryDirectory(prefix=f"q12o5867-p{prime}-") as directory:
    directory = Path(directory)
    source = directory / "shell.c"
    executable = directory / "shell"
    source.write_text(C_SOURCE)
    compile_result = subprocess.run(
        ["gcc", "-O3", "-std=c99", str(source), "-o", str(executable)],
        check=True, capture_output=True, text=True,
    )
    shell_started = time.monotonic()
    shell_result = subprocess.run(
        [str(executable)], check=True, capture_output=True, text=True,
    )
    shell_runtime = time.monotonic()-shell_started
    compiler_version = subprocess.run(
        ["gcc", "--version"], check=True, capture_output=True, text=True,
    ).stdout.splitlines()[0]


pointing = q8["preferred_pointed_zero"]
a1, a2, a3, unused_a4, unused_a6 = [
    record_function(record) for record in pointing["generalized_weierstrass_a_invariants"]
]
ordinate = record_function(pointing["quartic_ordinate"])
b2 = a1**2+4*a2
old_base_support = reduce_qq(pointing["old_base_coordinate"])


def inverse_parent_base(X, Y):
    x_general = K(X)/9-b2/12
    y_general = K(Y)/27-(a1*x_general+a3)/2
    if not y_general:
        return None
    return K(old_base_support+2*ordinate*(x_general+a2)/y_general)


def rational_degree(value):
    return max(value.numerator().degree(), value.denominator().degree())


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def coefficient_jacobian(X, Y):
    dx = -3*X**2-A
    dy = 2*Y
    columns = [padded(u**power*dx, 13) for power in range(5)]
    columns += [padded(u**power*dy, 13) for power in range(7)]
    return matrix(F, columns).transpose()


def component_profile(X, Y):
    result = []
    for support, node in zip(supports, nodes):
        result.append(int(X(support) == node and Y(support) == 0))
    infinity_x = X[4] if X.degree() == 4 else F.zero()
    infinity_y = Y[6] if Y.degree() == 6 else F.zero()
    result.append(int(infinity_x == infinity_node and infinity_y == 0))
    return result


records = []
seen = set()
for line in shell_result.stdout.splitlines():
    mode_text, x_text, y_text = line.split("|")
    X = R([F(int(v)) for v in x_text[2:].split(",")])
    Y = R([F(int(v)) for v in y_text[2:].split(",")])
    key = (tuple(X.list()), tuple(Y.list()))
    if key in seen:
        continue
    seen.add(key)
    assert Y**2 == X**3+A*X+B
    profile = component_profile(X, Y)
    assert profile == ([0, 0, 0, 0] if mode_text == "M0" else [1, 0, 0, 0])
    parent = inverse_parent_base(X, Y)
    records.append({
        "x_coefficients_low_to_high": list(map(int, X.list())),
        "y_coefficients_low_to_high": list(map(int, Y.list())),
        "equation_component_profile": profile,
        "inverse_parent_degree": None if parent is None else int(rational_degree(parent)),
        "ordinary_coefficient_jacobian_rank": int(coefficient_jacobian(X, Y).rank()),
    })

branch_specs = {
    "Q1": {"profile": [1, 0, 0, 0], "degree": 3, "marked_profile": [0, 1, 0, 0]},
    "Q2": {"profile": [0, 0, 0, 0], "degree": 2, "marked_profile": [0, 0, 0, 0]},
    "Q4": {"profile": [1, 0, 0, 0], "degree": 2, "marked_profile": [0, 1, 0, 0]},
}
branches = {}
for name, spec in branch_specs.items():
    candidates = [
        dict(record, shell_index=index)
        for index, record in enumerate(records)
        if record["equation_component_profile"] == spec["profile"]
        and record["inverse_parent_degree"] == spec["degree"]
    ]
    ranks = sorted({record["ordinary_coefficient_jacobian_rank"] for record in candidates})
    branches[name] = {
        "expected_equation_component_profile": spec["profile"],
        "expected_marked_component_profile": spec["marked_profile"],
        "expected_inverse_parent_degree": spec["degree"],
        "candidate_count": len(candidates),
        "rank_histogram": {
            str(rank): sum(record["ordinary_coefficient_jacobian_rank"] == rank for record in candidates)
            for rank in ranks
        },
        "candidates": candidates,
    }


def point_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = map(K, left)
    x2, y2 = map(K, right)
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+K(A))/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    y3 = slope*(x1-x3)-y1
    assert y3**2 == x3**3+K(A)*x3+K(B)
    return x3, y3


def point_neg(point):
    return point[0], -point[1]


def section_point(record):
    return R(record["x_coefficients_low_to_high"]), R(record["y_coefficients_low_to_high"])


def pole_shape(point, P_dot_O):
    if point is None:
        return False
    x_value, y_value = point
    return (
        x_value.numerator().degree() <= 4+2*P_dot_O
        and x_value.denominator().degree() <= 2*P_dot_O
        and y_value.numerator().degree() <= 6+3*P_dot_O
        and y_value.denominator().degree() <= 3*P_dot_O
    )


def directional_parent_degrees(point):
    if point is None:
        return None
    forward = inverse_parent_base(point[0], point[1])
    reverse = inverse_parent_base(point[0], -point[1])
    if forward is None or reverse is None:
        return None
    return [int(rational_degree(forward)), int(rational_degree(reverse))]


# Q3 is exact over QQ and regular at p=67.  These directional fingerprints
# are computed from the exact marked lattice classes; at a good prime they
# identify named reductions without assuming shell isolation.
q3_section = q3_exact["section"]
Q3_point = (
    R([reduce_qq(v) for v in q3_section["x_coefficients_low_to_high"]]),
    R([reduce_qq(v) for v in q3_section["y_coefficients_low_to_high"]]),
)
assert Q3_point[1]**2 == Q3_point[0]**3+A*Q3_point[0]+B
q3_rank = int(coefficient_jacobian(*Q3_point).rank())
q3_fingerprints = {
    "Q1": {"P_dot_O": 1, "directional_parent_degrees": [8, 4]},
    "Q2": {"P_dot_O": 1, "directional_parent_degrees": [7, 5]},
    "Q4": {"P_dot_O": 1, "directional_parent_degrees": [7, 5]},
}
for name, fingerprint in q3_fingerprints.items():
    survivors = []
    for index, record in enumerate(branches[name]["candidates"]):
        difference = point_add(section_point(record), point_neg(Q3_point))
        if not pole_shape(difference, fingerprint["P_dot_O"]):
            continue
        if directional_parent_degrees(difference) != fingerprint["directional_parent_degrees"]:
            continue
        survivors.append(index)
    branches[name]["exact_Q3_word_fingerprint"] = fingerprint
    branches[name]["Q3_fingerprint_survivor_indices"] = survivors
    branches[name]["Q3_fingerprint_survivor_count"] = len(survivors)

# Pairwise P.O=0 fingerprints further cross-check signs and named ordering.
pair_specs = {
    "Q1_minus_Q4": ("Q1", "Q4", [5, 3]),
    "Q1_minus_Q2": ("Q1", "Q2", [5, 3]),
    "Q2_minus_Q4": ("Q2", "Q4", [4, 4]),
}
pair_results = {}
for label, (left_name, right_name, degrees) in pair_specs.items():
    matches = []
    left_indices = branches[left_name]["Q3_fingerprint_survivor_indices"]
    right_indices = branches[right_name]["Q3_fingerprint_survivor_indices"]
    for left_index in left_indices:
        left = branches[left_name]["candidates"][left_index]
        for right_index in right_indices:
            right = branches[right_name]["candidates"][right_index]
            difference = point_add(section_point(left), point_neg(section_point(right)))
            if not pole_shape(difference, 0):
                continue
            if directional_parent_degrees(difference) != degrees:
                continue
            matches.append({f"{left_name}_candidate_index": left_index, f"{right_name}_candidate_index": right_index})
    pair_results[label] = {
        "P_dot_O": 0,
        "directional_parent_degrees": degrees,
        "match_count_after_Q3_filter": len(matches),
        "matches": matches,
    }

# p=67 is good for the child and Q3, but bad for direct use of the compact
# q4 equation basis: its exact xy scaling vanishes.  Keep this as an explicit
# construction diagnostic rather than silently treating the word route as run.
xy_scale_value = QQ(model["exact_coordinate_change"]["s"])
parent_word_direct_reduction_good = (
    xy_scale_value.numerator() % prime != 0
    and xy_scale_value.denominator() % prime != 0
)

fully_isolated = all(
    branches[name]["Q3_fingerprint_survivor_count"] == 1 for name in branch_specs
)
try:
    output_display = str(OUTPUT.relative_to(ROOT))
except ValueError:
    output_display = str(OUTPUT)
payload = {
    "schema": "elkies-k3.h92-q12o5867-p0-shell-word-fingerprints-modp.v1",
    "status": (
        "PASS_MODP_Q12O5867_NAMED_P0_BRANCH_SEEDS"
        if fully_isolated else
        "PASS_MODP_Q12O5867_P0_SHELL_WITH_RESIDUAL_NAMING_AMBIGUITY"
    ),
    "prime": int(prime),
    "model": {
        "A_coefficients_low_to_high": list(map(int, A.list())),
        "B_coefficients_low_to_high": list(map(int, B.list())),
        "equation_I2_supports_then_infinity": [int(v) for v in supports] + ["infinity"],
        "equation_I2_nodes": [int(v) for v in nodes] + [int(infinity_node)],
        "embedding_15_marked_0100_to_equation_1000": True,
    },
    "direct_shell": {
        "unique_signed_section_count": len(records),
        "profile_0000_count": sum(record["equation_component_profile"] == [0, 0, 0, 0] for record in records),
        "profile_1000_count": sum(record["equation_component_profile"] == [1, 0, 0, 0] for record in records),
        "generated_C_helper_stderr": shell_result.stderr.splitlines(),
        "enumeration_runtime_seconds": shell_runtime,
        "method": "five-functional interpolation and exact univariate polynomial-square test",
    },
    "branches": branches,
    "pairwise_exact_word_fingerprints": pair_results,
    "exact_Q3_anchor": {
        "x_coefficients_low_to_high": list(map(int, Q3_point[0].list())),
        "y_coefficients_low_to_high": list(map(int, Q3_point[1].list())),
        "ordinary_coefficient_jacobian_rank": q3_rank,
    },
    "parent_equation_word_reduction_gate": {
        "good": bool(parent_word_direct_reduction_good),
        "reason_if_bad": (
            None if parent_word_direct_reduction_good else
            "the exact q4/o164 compact-model xy scaling is zero modulo p; direct B-basis word evaluation is not valid at this prime"
        ),
    },
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "QQ_lift_attempted": False,
        "sage_version": SAGE_VERSION,
        "C_compiler": compiler_version,
        "reproducing_command": (
            "sage -python elkies-k3/scripts/construct_h92_q12o5867_p0_shell_modp.sage "
            f"--prime {prime} "
            f"{'--include-all-records ' if args.include_all_records else ''}"
            f"--output {output_display}"
        ),
        "runtime_seconds": time.monotonic()-started,
    },
    "residual_ambiguity": {
        name: branches[name]["Q3_fingerprint_survivor_count"]
        for name in branch_specs
        if branches[name]["Q3_fingerprint_survivor_count"] != 1
    },
    "proof_boundary": (
        "This is an exact bounded computation over the displayed finite field. "
        "Every shell point is literally checked on the exact reduced P1229-pointed q8 child, "
        "and all degree/profile/Jacobian and group-law fingerprints are exact modulo p. "
        "Residual modular naming ambiguity is retained explicitly. No characteristic-zero "
        "section lift or q12 resolved Riemann--Roch plane is claimed."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
if args.include_all_records:
    payload["all_records"] = records
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867P0SHELL|prime={}|shell={}|Q1={}/{}|Q2={}/{}|Q4={}/{}|"
    "runtime={:.3f}|status={}|output={}".format(
        prime, len(records),
        branches["Q1"]["candidate_count"], branches["Q1"]["Q3_fingerprint_survivor_count"],
        branches["Q2"]["candidate_count"], branches["Q2"]["Q3_fingerprint_survivor_count"],
        branches["Q4"]["candidate_count"], branches["Q4"]["Q3_fingerprint_survivor_count"],
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
