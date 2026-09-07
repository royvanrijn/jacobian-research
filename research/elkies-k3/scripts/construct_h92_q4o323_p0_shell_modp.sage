#!/usr/bin/env sage -python
"""Enumerate the polynomial P.O=0 shell on the q4/o323 child modulo p.

status: ACTIVE_SEARCH
claim: complete bounded polynomial P.O=0 shell at a fibre-safe prime
inputs: exact component-2-pointed q4/o323 A3+2A2 child equation
outputs: artifacts/local/elkies-k3/q4o323-p0-shell-mod<prime>.json

This is the small finite-field preflight for the reflected physical q12 edge.
It uses five interpolation functionals and a univariate polynomial-square
test.  It deliberately does not identify resolved I3/I4 component signs yet.
No Groebner basis or elimination is used.
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
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o208-q4o323-a3-2a2-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=61)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
OUTPUT = (
    args.output.resolve() if args.output else
    LOCAL / f"q4o323-p0-shell-mod{prime}.json"
)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce_qq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % prime == 0:
        raise ZeroDivisionError(f"bad reduction denominator at p={prime}")
    return F(value.numerator()) / F(value.denominator())


data = json.loads(MODEL.read_text())
assert data["status"] == "PASS_EXACT_QQ_Q4O323_A3_2A2_RR_AND_JACOBIAN"
child = data["child"]
F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
A = R([reduce_qq(value) for value in child["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in child["minimal_B_coefficients_low_to_high"]])
assert (A.degree(), B.degree()) == (8, 12)

# The two displayed linear finite factors are the I3 supports.  Add two
# ordinary evaluation points; the fifth functional is the leading coefficient.
supports = []
for fibre in child["finite_fibres"]:
    if fibre["kodaira"] != "I3":
        continue
    c0, c1 = map(reduce_qq, fibre["factor_coefficients_low_to_high"])
    supports.append(-c0/c1)
assert len(supports) == 2 and len(set(supports)) == 2

# Preserve the marked fibre types under specialization.  Several superficially
# usable small primes (including 31) make a residual discriminant factor collide
# with an I3 support and must not be used for component matching.
discriminant = -16*(4*A**3+27*B**2)
support_valuations = [int(discriminant.valuation(u-support)) for support in supports]
good_fibre_specialization = discriminant.degree() == 20 and support_valuations == [3, 3]
if not good_fibre_specialization:
    payload = {
        "schema": "elkies-k3.h92-q4o323-p0-shell-modp.v1",
        "status": "REJECTED_MODP_Q4O323_BAD_FIBRE_SPECIALIZATION",
        "prime": int(prime),
        "fibre_gate": {
            "discriminant_degree": int(discriminant.degree()),
            "finite_I3_support_valuations": support_valuations,
            "expected": {"discriminant_degree": 20, "support_valuations": [3, 3]},
        },
        "shell": {"run": False},
        "method": {
            "large_Groebner_required": False,
            "elimination_required": False,
            "runtime_seconds": time.monotonic()-started,
        },
        "proof_boundary": "Bad-prime rejection only; no section shell was enumerated.",
        "inputs": {
            "paths": [str(MODEL.relative_to(ROOT))],
            "sha256": {str(MODEL.relative_to(ROOT)): sha256(MODEL)},
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(
        "Q4O323P0SHELL|prime={}|valuations={}|delta_degree={}|status={}|output={}".format(
            prime, support_valuations, discriminant.degree(), payload["status"], OUTPUT,
        ), flush=True,
    )
    sys.exit(0)
generic = [F(i) for i in range(int(prime)) if F(i) not in supports][:2]
evaluation_points = supports + generic
interpolation = matrix(F, [
    [point**degree for degree in range(5)]
    for point in evaluation_points
] + [[0, 0, 0, 0, 1]]).inverse()


def c_array(values):
    return ",".join(str(int(value)) for value in values)


C_SOURCE = r'''
#include <stdio.h>
#include <stdint.h>
#define P @P@
static int A[9]={@A@}, B[13]={@B@};
static int pts[4]={@PTS@};
static int invm[5][5]={@INVM@};
static int md(long long x){x%=P;if(x<0)x+=P;return (int)x;}
static int pw(int a,int n){int r=1;while(n){if(n&1)r=md((long long)r*a);a=md((long long)a*a);n>>=1;}return r;}
static int val(const int *f,int d,int z){int r=0;for(int i=d;i>=0;i--)r=md((long long)r*z+f[i]);return r;}
static int square(int a){return a==0 || pw(a,(P-1)/2)==1;}
static int sqrt0(int a){for(int y=0;y<P;y++)if(md((long long)y*y)==a)return y;return -1;}
static int rootsquare(const int *r,int *y){
 int d=12;while(d>=0 && r[d]==0)d--;if(d<0){for(int i=0;i<7;i++)y[i]=0;return 1;}
 if(d&1)return 0;int m=d/2,s=sqrt0(r[d]);if(s<0)return 0;
 for(int i=0;i<7;i++)y[i]=0;y[m]=s;int iden=pw(md(2*s),P-2);
 for(int k=d-1;k>=m;k--){int j=k-m;long long known=0;for(int a=j+1;a<=m;a++){int b=k-a;if(b>=0&&b<=m)known+=(long long)y[a]*y[b];}y[j]=md((long long)(r[k]-md(known))*iden);}
 for(int k=0;k<=12;k++){long long q=0;for(int i=0;i<=6;i++){int j=k-i;if(j>=0&&j<=6)q+=(long long)y[i]*y[j];}if(md(q)!=r[k])return 0;}return 1;
}
static void emit(int *x,int *y){printf("X=");for(int i=0;i<5;i++)printf("%s%d",i?",":"",x[i]);printf("|Y=");for(int i=0;i<7;i++)printf("%s%d",i?",":"",y[i]);puts("");}
int main(void){
 int allowed[5][P],n[5]={0};
 for(int j=0;j<5;j++)for(int x=0;x<P;x++){
  int rhs=j<4?md((long long)x*x%P*x+(long long)val(A,8,pts[j])*x+val(B,12,pts[j])):md((long long)x*x%P*x+(long long)A[8]*x+B[12]);
  if(square(rhs))allowed[j][n[j]++]=x;
 }
 fprintf(stderr,"lists=%d,%d,%d,%d,%d\n",n[0],n[1],n[2],n[3],n[4]);
 for(int i0=0;i0<n[0];i0++)for(int i1=0;i1<n[1];i1++)for(int i2=0;i2<n[2];i2++)for(int i3=0;i3<n[3];i3++)for(int i4=0;i4<n[4];i4++){
  int vv[5]={allowed[0][i0],allowed[1][i1],allowed[2][i2],allowed[3][i3],allowed[4][i4]},x[5],x2[9]={0},x3[13]={0},rhs[13],y[7];
  for(int a=0;a<5;a++){long long s=0;for(int b=0;b<5;b++)s+=(long long)invm[a][b]*vv[b];x[a]=md(s);}
  for(int i=0;i<5;i++)for(int j=0;j<5;j++)x2[i+j]=md(x2[i+j]+(long long)x[i]*x[j]);
  for(int i=0;i<9;i++)for(int j=0;j<5;j++)x3[i+j]=md(x3[i+j]+(long long)x2[i]*x[j]);
  for(int k=0;k<13;k++){long long s=x3[k]+B[k];for(int i=0;i<9;i++){int j=k-i;if(j>=0&&j<5)s+=(long long)A[i]*x[j];}rhs[k]=md(s);}
  if(!rootsquare(rhs,y))continue;emit(x,y);int nz=0;for(int i=0;i<7;i++)nz|=y[i];if(nz){for(int i=0;i<7;i++)y[i]=md(-y[i]);emit(x,y);}
 }
 return 0;
}
'''
C_SOURCE = C_SOURCE.replace("@P@", str(int(prime)))
C_SOURCE = C_SOURCE.replace("@A@", c_array(A[i] for i in range(9)))
C_SOURCE = C_SOURCE.replace("@B@", c_array(B[i] for i in range(13)))
C_SOURCE = C_SOURCE.replace("@PTS@", c_array(evaluation_points))
C_SOURCE = C_SOURCE.replace("@INVM@", ",".join(
    "{" + c_array(row) + "}" for row in interpolation.rows()
))

with tempfile.TemporaryDirectory(prefix=f"q4o323-p{prime}-") as directory:
    directory = Path(directory)
    source = directory / "shell.c"
    executable = directory / "shell"
    source.write_text(C_SOURCE)
    subprocess.run(
        ["gcc", "-O3", "-std=c99", str(source), "-o", str(executable)],
        check=True, capture_output=True, text=True,
    )
    shell_started = time.monotonic()
    result = subprocess.run([str(executable)], check=True, capture_output=True, text=True)
    shell_runtime = time.monotonic()-shell_started
    compiler = subprocess.run(
        ["gcc", "--version"], check=True, capture_output=True, text=True,
    ).stdout.splitlines()[0]


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def jacobian_rank(X, Y):
    columns = [padded(u**power*(-3*X**2-A), 13) for power in range(5)]
    columns += [padded(u**power*(2*Y), 13) for power in range(7)]
    return matrix(F, columns).transpose().rank()


# Node incidence is only a coarse component invariant for I3/I4.  The signed
# resolved-component labels will be attached in the next, target-specific pass.
RX = PolynomialRing(F, "x")
xvar = RX.gen()
nodes = []
for support in supports:
    cubic = xvar**3+A(support)*xvar+B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    nodes.append(-repeated[0]/repeated[1])
infinity_cubic = xvar**3+A[8]*xvar+B[12]
repeated = infinity_cubic.gcd(infinity_cubic.derivative())
assert repeated.degree() == 1
infinity_node = -repeated[0]/repeated[1]

records = []
seen = set()
for line in result.stdout.splitlines():
    x_text, y_text = line.split("|")
    X = R([F(int(value)) for value in x_text[2:].split(",")])
    Y = R([F(int(value)) for value in y_text[2:].split(",")])
    key = (tuple(X.list()), tuple(Y.list()))
    if key in seen:
        continue
    seen.add(key)
    assert Y**2 == X**3+A*X+B
    profile = [int(X(s)==node and Y(s)==0) for s, node in zip(supports, nodes)]
    profile.append(int(X.degree() == 4 and Y.degree() == 6 and X[4] == infinity_node and Y[6] == 0))
    records.append({
        "x_coefficients_low_to_high": list(map(int, X.list())),
        "y_coefficients_low_to_high": list(map(int, Y.list())),
        "node_incidence_finite_I3_then_infinity_I4": profile,
        "ordinary_coefficient_jacobian_rank": int(jacobian_rank(X, Y)),
    })

rank_histogram = {
    str(rank): sum(record["ordinary_coefficient_jacobian_rank"] == rank for record in records)
    for rank in sorted({record["ordinary_coefficient_jacobian_rank"] for record in records})
}
payload = {
    "schema": "elkies-k3.h92-q4o323-p0-shell-modp.v1",
    "status": "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL",
    "prime": int(prime),
    "model": {
        "A_coefficients_low_to_high": list(map(int, A.list())),
        "B_coefficients_low_to_high": list(map(int, B.list())),
        "finite_I3_supports": list(map(int, supports)),
        "finite_I3_nodes": list(map(int, nodes)),
        "infinity_I4_node": int(infinity_node),
    },
    "shell": {
        "unique_signed_section_count": len(records),
        "ordinary_coefficient_jacobian_rank_histogram": rank_histogram,
        "records": records,
        "C_helper_stderr": result.stderr.splitlines(),
        "enumeration_runtime_seconds": shell_runtime,
    },
    "method": {
        "description": "five-functional interpolation and exact univariate polynomial-square test",
        "large_Groebner_required": False,
        "elimination_required": False,
        "sage_version": SAGE_VERSION,
        "C_compiler": compiler,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "This is a complete bounded polynomial P.O=0 shell over the displayed finite field. "
        "Every record is checked on the reduced exact q4/o323 child and has its ordinary "
        "coefficient-Jacobian rank computed. Resolved I3/I4 component signs, marked NS "
        "classes, characteristic-zero lifts, and the physical q12 horizontal remain open."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): sha256(MODEL)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323P0SHELL|prime={}|signed={}|ranks={}|runtime={:.3f}|status={}|output={}".format(
        prime, len(records), rank_histogram, payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
