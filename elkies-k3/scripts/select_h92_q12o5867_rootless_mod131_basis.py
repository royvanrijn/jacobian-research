#!/usr/bin/env python3
"""Select a rank-17 short-vector basis from the q12/o5867 mod-131 shell.

status: ACTIVE_PROOF
claim: a regular modular section basis with the pinned short-basis height Gram
inputs: q12o5867-rootless-p0-shell-mod131.json, pinned lattice files
outputs: q12o5867-rootless-mod131-selected-basis.json

For polynomial sections on the rootless model, <P,Q>=2-P.Q.  The intersection
P.Q is the finite gcd degree of (x_P-x_Q,y_P-y_Q), plus the analogous minimum
order at infinity.  A tiny OpenMP C helper computes this exact finite-field
pairing matrix; bitset-style MRV backtracking matches the pinned short basis.
No Groebner basis or elimination is used.
"""

import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SHELL = LOCAL / "q12o5867-rootless-p0-shell-mod131.json"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SIGNED_SHORT = ROOT / "elkies-k3/data/relations/all_2622_signed_short_basis.npy"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = LOCAL / "q12o5867-rootless-mod131-selected-basis.json"
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_integer_matrix(path):
    return np.loadtxt(path, dtype=np.int64)


shell = json.loads(SHELL.read_text())
assert shell["status"] == "PASS_COMPLETE_MOD131_Q12O5867_ROOTLESS_P0_SHELL"
assert shell["prime"] == 131
records = shell["records"]
n = len(records)
assert n == 2622
assert all(record["coefficient_jacobian_rank"] == 12 for record in records)
short_gram = load_integer_matrix(SHORT_GRAM)
short_coords = load_integer_matrix(SHORT_COORDS)
signed_short = np.load(SIGNED_SHORT).astype(np.int64, copy=False)
# A deterministic unimodular minimal-vector basis found by a rare-signature
# greedy rank-16 selection followed by primitive determinant-one completion.
target_signed_indices = [
    0, 2, 74, 328, 78, 12, 170, 376, 290,
    164, 1120, 6, 300, 316, 326, 4, 14,
]
target_vectors = signed_short[target_signed_indices]
target = target_vectors @ short_gram @ target_vectors.T
coords = target_vectors @ short_coords
pinned = load_integer_matrix(PINNED)
assert target.shape == short_gram.shape == coords.shape == pinned.shape == (17, 17)
assert round(np.linalg.det(target_vectors)) == 1
assert round(np.linalg.det(coords)) == -1
assert np.array_equal(coords @ pinned @ coords.T, target)
assert round(np.linalg.det(target)) == 948


def c_rows(name, width):
    return ",\n".join(
        "{"+",".join(map(str, record[name]))+"}" for record in records
    )


C_SOURCE = r'''
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#define P 131
#define N @N@
static uint8_t X[N][5]={@X@};
static uint8_t Y[N][7]={@Y@};
static int md(int x){x%=P;if(x<0)x+=P;return x;}
static int inv(int a){int r=1,n=P-2;while(n){if(n&1)r=md(r*a);a=md(a*a);n>>=1;}return r;}
static int degree(int *a,int cap){while(cap>=0&&a[cap]==0)cap--;return cap;}
static int gcd_degree(int *aa,int da,int *bb,int db){
 int a[7]={0},b[7]={0},r[7]={0};for(int i=0;i<=da;i++)a[i]=aa[i];for(int i=0;i<=db;i++)b[i]=bb[i];
 while(db>=0){for(int i=0;i<7;i++)r[i]=a[i];int dr=da,ib=inv(b[db]);while(dr>=db){int c=md(r[dr]*ib);int shift=dr-db;for(int j=0;j<=db;j++)r[j+shift]=md(r[j+shift]-c*b[j]);dr=degree(r,dr-1);}for(int i=0;i<7;i++){a[i]=b[i];b[i]=r[i];}da=db;db=dr;}
 return da;
}
int main(int argc,char **argv){
 if(argc!=2)return 2;FILE *out=fopen(argv[1],"wb");if(!out)return 3;int8_t *M=(int8_t*)malloc((size_t)N*N);if(!M)return 4;
 #pragma omp parallel for schedule(dynamic)
 for(int i=0;i<N;i++){
  M[(size_t)i*N+i]=4;
  for(int j=i+1;j<N;j++){
   int dx[7]={0},dy[7]={0};for(int k=0;k<5;k++)dx[k]=md((int)X[i][k]-(int)X[j][k]);for(int k=0;k<7;k++)dy[k]=md((int)Y[i][k]-(int)Y[j][k]);
   int ddx=degree(dx,4),ddy=degree(dy,6);int finite=gcd_degree(dx,ddx,dy,ddy);int ox=(ddx<0)?99:4-ddx,oy=(ddy<0)?99:6-ddy;int infinity=(ox<oy)?ox:oy;int pairing=2-finite-infinity;
   M[(size_t)i*N+j]=(int8_t)pairing;M[(size_t)j*N+i]=(int8_t)pairing;
  }
 }
 fwrite(M,1,(size_t)N*N,out);fclose(out);free(M);return 0;
}
'''
C_SOURCE = C_SOURCE.replace("@N@", str(n))
C_SOURCE = C_SOURCE.replace("@X@", c_rows("x_coefficients_low_to_high", 5))
C_SOURCE = C_SOURCE.replace("@Y@", c_rows("y_coefficients_low_to_high", 7))

with tempfile.TemporaryDirectory(prefix="q12o5867-pairing-p131-") as directory:
    directory = Path(directory)
    source = directory / "pairing.c"
    executable = directory / "pairing"
    binary = directory / "pairing.bin"
    source.write_text(C_SOURCE)
    subprocess.run(
        ["gcc", "-O3", "-march=native", "-fopenmp", "-std=c99", str(source), "-o", str(executable)],
        check=True, capture_output=True, text=True,
    )
    pairing_started = time.monotonic()
    subprocess.run(
        [str(executable), str(binary)], check=True,
        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", str(os.cpu_count() or 1))},
    )
    pairing_seconds = time.monotonic()-pairing_started
    pairing = np.fromfile(binary, dtype=np.int8).reshape((n, n)).copy()

assert np.array_equal(pairing, pairing.T)
assert np.all(np.diag(pairing) == 4)
assert set(np.unique(pairing)).issubset({-4, -3, -2, -1, 0, 1, 2, 4})


# The most distinctive target row is used first. Candidate anchors are tried
# deterministically; the exact pair constraints prune the remaining domains.
def profile(row, skip):
    values = np.delete(row, skip)
    return tuple((value, int(np.sum(values == value))) for value in sorted(set(values.tolist())))


target_profiles = [profile(target[index], index) for index in range(17)]
anchor_target = 0
nodes = 0


def search(assignment, used):
    global nodes
    nodes += 1
    if len(assignment) == 17:
        return assignment
    best_vertex = min(vertex for vertex in range(17) if vertex not in assignment)
    domain = np.ones(n, dtype=bool)
    for old_vertex, candidate in assignment.items():
        domain &= pairing[:, candidate] == target[best_vertex, old_vertex]
    if used:
        domain[list(used)] = False
    best_domain = np.flatnonzero(domain)
    if best_domain is None or not len(best_domain):
        return None
    for candidate in best_domain:
        assignment[best_vertex] = int(candidate)
        answer = search(assignment, used | {int(candidate)})
        if answer is not None:
            return answer
        del assignment[best_vertex]
    return None


solution = None
anchor_trials = 0
target_anchor_shell_pairings = signed_short @ short_gram @ target_vectors[anchor_target]
target_anchor_histogram = {
    int(value): int(np.sum(target_anchor_shell_pairings == value))
    for value in np.unique(target_anchor_shell_pairings)
}
for anchor_candidate in range(n):
    candidate_histogram = {
        int(value): int(np.sum(pairing[anchor_candidate] == value))
        for value in np.unique(pairing[anchor_candidate])
    }
    if candidate_histogram != target_anchor_histogram:
        continue
    anchor_trials += 1
    solution = search({anchor_target: anchor_candidate}, {anchor_candidate})
    if solution is not None:
        break
assert solution is not None
selected_indices = [solution[index] for index in range(17)]
selected_gram = pairing[np.ix_(selected_indices, selected_indices)].astype(np.int64)
assert np.array_equal(selected_gram, target)

selected = []
for basis_index, record_index in enumerate(selected_indices):
    selected.append({
        "basis_index": basis_index,
        "shell_record_index": record_index,
        "x_coefficients_low_to_high": records[record_index]["x_coefficients_low_to_high"],
        "y_coefficients_low_to_high": records[record_index]["y_coefficients_low_to_high"],
        "coefficient_jacobian_rank": records[record_index]["coefficient_jacobian_rank"],
    })

payload = {
    "schema": "elkies-k3.h92-q12o5867-rootless-modp-selected-basis.v1",
    "status": "PASS_MOD131_Q12O5867_ROOTLESS_REGULAR_SHORT_BASIS",
    "prime": 131,
    "signed_shell_count": n,
    "selected_sections": selected,
    "height_gram": selected_gram.tolist(),
    "height_gram_determinant": 948,
    "target_signed_short_vector_indices": target_signed_indices,
    "short_basis_to_pinned_basis": coords.tolist(),
    "short_basis_to_pinned_determinant": -1,
    "search": {
        "anchor_target_index": anchor_target,
        "anchor_trials": anchor_trials,
        "dfs_nodes": nodes,
        "pairing_matrix_seconds": pairing_seconds,
        "runtime_seconds": time.monotonic()-started,
    },
    "method": {
        "pairing_formula": "<P,Q>=2-deg(gcd(xP-xQ,yP-yQ))-min(ord_infinity(xP-xQ),ord_infinity(yP-yQ))",
        "large_Groebner_required": False,
        "elimination_required": False,
    },
    "proof_boundary": (
        "This selects a regular rank-17 modular basis with the exact pinned short-basis "
        "Gram. Characteristic-zero lifting, literal QQ identities, and saturation are separate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SHELL, SHORT_GRAM, SHORT_COORDS, SIGNED_SHORT, PINNED)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SHELL, SHORT_GRAM, SHORT_COORDS, SIGNED_SHORT, PINNED)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ROOTLESSBASIS|prime=131|sections=17|det=948|anchor_trials={}|nodes={}|"
    "seconds={:.3f}|status={}|output={}".format(
        anchor_trials, nodes, payload["search"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
