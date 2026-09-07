#!/usr/bin/env sage -python
"""Enumerate the terminal q12/o5867 polynomial P.O=0 shell modulo 131.

status: ACTIVE_PROOF
claim: complete good-reduction polynomial shell on the exact rootless child
inputs: artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json
outputs: artifacts/local/elkies-k3/q12o5867-rootless-p0-shell-mod131.json

Four finite evaluations and the leading coefficient determine X of degree at
most four.  A generated OpenMP C loop tests whether X^3+A*X+B is the square of
a polynomial of degree at most six.  Sage then verifies every returned section
and its ordinary 13-by-12 coefficient Jacobian.  This is exhaustive finite
evaluation, not elimination; no Groebner basis is used.
"""

import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
OUTPUT = LOCAL / "q12o5867-rootless-p0-shell-mod131.json"
PRIME = ZZ(131)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


data = json.loads(MODEL.read_text())
assert data["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
child = data["child"]
F = GF(PRIME)
R = PolynomialRing(F, "u")
u = R.gen()


def reduce_qq(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % PRIME == 0:
        raise ZeroDivisionError(f"bad reduction denominator at p={PRIME}")
    return F(value.numerator()) / F(denominator)


A = R([reduce_qq(value) for value in child["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in child["minimal_B_coefficients_low_to_high"]])
Delta = -16 * (4*A**3 + 27*B**2)
assert (A.degree(), B.degree(), Delta.degree()) == (8, 12, 24)
assert Delta.gcd(Delta.derivative()).degree() == 0


def is_square(value):
    return value == 0 or value**((PRIME-1)//2) == 1


def admissible_x_count(support):
    coefficient_A = A(support)
    coefficient_B = B(support)
    return sum(
        1 for xvalue in F
        if is_square(xvalue**3 + coefficient_A*xvalue + coefficient_B)
    )


# Choose the four cheapest distinct finite interpolation functionals.  The
# fifth functional is the leading coefficient and sees the infinity cubic.
finite_counts = sorted(
    (admissible_x_count(F(integer)), integer)
    for integer in range(int(PRIME))
)
evaluation_points = [F(integer) for unused_count, integer in finite_counts[:4]]
interpolation = matrix(F, [
    [point**degree for degree in range(5)] for point in evaluation_points
] + [[0, 0, 0, 0, 1]]).inverse()


def c_array(values):
    return ",".join(str(int(value)) for value in values)


C_SOURCE = r'''
#include <stdio.h>
#include <stdint.h>
#include <omp.h>
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
 int d=12;while(d>=0&&r[d]==0)d--;if(d<0){for(int i=0;i<7;i++)y[i]=0;return 1;}
 if(d&1)return 0;int m=d/2,s=sqrt0(r[d]);if(s<0)return 0;
 for(int i=0;i<7;i++)y[i]=0;y[m]=s;int iden=pw(md(2*s),P-2);
 for(int k=d-1;k>=m;k--){int j=k-m;long long known=0;for(int a=j+1;a<=m;a++){int b=k-a;if(b>=0&&b<=m)known+=(long long)y[a]*y[b];}y[j]=md((long long)(r[k]-md(known))*iden);}
 for(int k=0;k<=12;k++){long long q=0;for(int i=0;i<=6;i++){int j=k-i;if(j>=0&&j<=6)q+=(long long)y[i]*y[j];}if(md(q)!=r[k])return 0;}return 1;
}
static void emit(int *x,int *y){
 #pragma omp critical
 {
  printf("X=");for(int i=0;i<5;i++)printf("%s%d",i?",":"",x[i]);
  printf("|Y=");for(int i=0;i<7;i++)printf("%s%d",i?",":"",y[i]);puts("");
 }
}
int main(void){
 int allowed[5][P],n[5]={0};
 for(int j=0;j<5;j++)for(int x=0;x<P;x++){
  int rhs=(j<4)?md((long long)x*x%P*x+(long long)val(A,8,pts[j])*x+val(B,12,pts[j])):md((long long)x*x%P*x+(long long)A[8]*x+B[12]);
  if(square(rhs))allowed[j][n[j]++]=x;
 }
 fprintf(stderr,"LISTS=%d,%d,%d,%d,%d\n",n[0],n[1],n[2],n[3],n[4]);
 #pragma omp parallel for schedule(dynamic)
 for(int i0=0;i0<n[0];i0++)for(int i1=0;i1<n[1];i1++)for(int i2=0;i2<n[2];i2++)for(int i3=0;i3<n[3];i3++)for(int i4=0;i4<n[4];i4++){
  int vv[5]={allowed[0][i0],allowed[1][i1],allowed[2][i2],allowed[3][i3],allowed[4][i4]},x[5],x2[9]={0},x3[13]={0},rhs[13],y[7];
  for(int a=0;a<5;a++){long long s=0;for(int b=0;b<5;b++)s+=(long long)invm[a][b]*vv[b];x[a]=md(s);}
  for(int i=0;i<5;i++)for(int j=0;j<5;j++)x2[i+j]=md(x2[i+j]+(long long)x[i]*x[j]);
  for(int i=0;i<9;i++)for(int j=0;j<5;j++)x3[i+j]=md(x3[i+j]+(long long)x2[i]*x[j]);
  for(int k=0;k<13;k++){long long s=x3[k]+B[k];for(int i=0;i<9;i++){int j=k-i;if(j>=0&&j<5)s+=(long long)A[i]*x[j];}rhs[k]=md(s);}
  if(!rootsquare(rhs,y))continue;emit(x,y);int nonzero=0;for(int i=0;i<7;i++)nonzero|=y[i];if(nonzero){for(int i=0;i<7;i++)y[i]=md(-y[i]);emit(x,y);}
 }
 return 0;
}
'''
C_SOURCE = C_SOURCE.replace("@P@", str(int(PRIME)))
C_SOURCE = C_SOURCE.replace("@A@", c_array([A[index] for index in range(9)]))
C_SOURCE = C_SOURCE.replace("@B@", c_array([B[index] for index in range(13)]))
C_SOURCE = C_SOURCE.replace("@PTS@", c_array(evaluation_points))
C_SOURCE = C_SOURCE.replace("@INVM@", ",".join(
    "{"+c_array(row)+"}" for row in interpolation.rows()
))

with tempfile.TemporaryDirectory(prefix="q12o5867-rootless-p131-") as directory:
    directory = Path(directory)
    source = directory / "shell.c"
    executable = directory / "shell"
    source.write_text(C_SOURCE)
    subprocess.run(
        ["gcc", "-O3", "-march=native", "-fopenmp", "-std=c99", str(source), "-o", str(executable)],
        check=True, capture_output=True, text=True,
    )
    shell_started = time.monotonic()
    result = subprocess.run(
        [str(executable)], check=True, capture_output=True, text=True,
        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", str(os.cpu_count() or 1))},
    )
    shell_runtime = time.monotonic()-shell_started
    compiler = subprocess.run(
        ["gcc", "--version"], check=True, capture_output=True, text=True,
    ).stdout.splitlines()[0]


def padded(poly, length):
    return [poly[index] if index <= poly.degree() else F.zero() for index in range(length)]


def jacobian_rank(X, Y):
    dx = -3*X**2-A
    dy = 2*Y
    columns = [padded(u**power*dx, 13) for power in range(5)]
    columns += [padded(u**power*dy, 13) for power in range(7)]
    return int(matrix(F, columns).transpose().rank())


records = []
seen = set()
for line in result.stdout.splitlines():
    xtext, ytext = line.split("|", 1)
    xcoeffs = tuple(int(value) for value in xtext.removeprefix("X=").split(","))
    ycoeffs = tuple(int(value) for value in ytext.removeprefix("Y=").split(","))
    if (xcoeffs, ycoeffs) in seen:
        continue
    seen.add((xcoeffs, ycoeffs))
    X = R([F(value) for value in xcoeffs])
    Y = R([F(value) for value in ycoeffs])
    assert Y**2 == X**3+A*X+B
    records.append({
        "x_coefficients_low_to_high": list(xcoeffs),
        "y_coefficients_low_to_high": list(ycoeffs),
        "coefficient_jacobian_rank": jacobian_rank(X, Y),
    })
records.sort(key=lambda row: (row["x_coefficients_low_to_high"], row["y_coefficients_low_to_high"]))
rank_counts = {}
for record in records:
    rank = str(record["coefficient_jacobian_rank"])
    rank_counts[rank] = rank_counts.get(rank, 0)+1

payload = {
    "schema": "elkies-k3.h92-q12o5867-rootless-p0-shell-modp.v1",
    "status": "PASS_COMPLETE_MOD131_Q12O5867_ROOTLESS_P0_SHELL",
    "prime": int(PRIME),
    "good_reduction": {
        "degrees_A_B_Delta": [int(A.degree()), int(B.degree()), int(Delta.degree())],
        "squarefree_discriminant": True,
        "infinity_smooth": True,
    },
    "interpolation": {
        "finite_points": [int(value) for value in evaluation_points],
        "finite_admissible_x_counts": [admissible_x_count(value) for value in evaluation_points],
        "C_stderr": result.stderr.strip(),
    },
    "signed_section_count": len(records),
    "unoriented_section_pair_count": len(records)//2,
    "coefficient_jacobian_rank_counts": rank_counts,
    "records": records,
    "method": {
        "C_compiler": compiler,
        "openmp_threads": int(os.environ.get("OMP_NUM_THREADS", os.cpu_count() or 1)),
        "shell_runtime_seconds": shell_runtime,
        "runtime_seconds": time.monotonic()-started,
        "large_Groebner_required": False,
        "elimination_required": False,
        "sage_version": SAGE_VERSION,
    },
    "proof_boundary": (
        "Complete polynomial P.O=0 shell at one exact good prime. This does not "
        "lift sections to QQ or prove a saturated characteristic-zero MW basis."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): sha256(MODEL)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ROOTLESSSHELL|prime={}|signed={}|pairs={}|ranks={}|seconds={:.3f}|"
    "status={}|output={}".format(
        PRIME, len(records), len(records)//2, rank_counts,
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
