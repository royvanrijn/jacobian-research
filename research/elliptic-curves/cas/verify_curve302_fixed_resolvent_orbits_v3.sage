#!/usr/bin/env sage -python
"""Independent exact replay of the retained v3 fixed-resolvent pairs."""

from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from math import gcd
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
TARGET = ART / "curve302_fixed_resolvent_orbits_v3.json"
STRICT = ART / "rank_jump_curve302_strict_constructor_arithmetic_v1.json"
OUTPUT = ART / "curve302_fixed_resolvent_orbits_v3_verification_v1.json"
SOURCE = Path(__file__).with_name("curve302_fixed_resolvent_orbits_v3.sage")
PROTOCOL = Path(__file__).with_name("CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL_V3.json")


def read(path): return json.loads(path.read_text())
def digest(path): return sha256(path.read_bytes()).hexdigest()
def relative(path): return str(path.relative_to(ROOT))
def put_new(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True); handle.write("\n")


def discriminant(a, b, c, d):
    return b*b*c*c - 4*a*c**3 - 4*b**3*d - 27*a*a*d*d + 18*a*b*c*d


def recovered_maximal_form():
    data = read(STRICT); R = PolynomialRing(QQ, "z"); f = R(data["cubic_ascending"])
    nf = pari.nfinit([pari(f), data["S_finite"]]); zk = list(nf.nf_get_zk()); assert zk[0] == 1
    pair = list(map(ZZ, pari.nfalgtobasis(nf, zk[1]*zk[2])))
    w, t = zk[1]-pair[2], zk[2]-pair[1]
    def coords(value):
        row = list(map(ZZ, pari.nfalgtobasis(nf, value)))
        return [row[0]+row[1]*pair[2]+row[2]*pair[1], row[1], row[2]]
    table = [[coords(u*v) for v in (1,w,t)] for u in (1,w,t)]
    form = [-table[1][1][2], table[1][1][1], -table[2][2][2], table[2][2][1]]
    assert discriminant(*form) == ZZ(nf.disc()) and ZZ(nf[3]) > 1
    return f, nf, form


def form_matrices(values):
    assert len(values) == 12
    values = list(map(QQ, values))
    def build(offset):
        a11,a22,a33,a12,a23,a13 = values[offset:offset+6]
        return matrix(QQ, [[a11,a12/2,a13/2],[a12/2,a22,a23/2],[a13/2,a23/2,a33]])
    return build(0), build(6)


def sign(values):
    assert sorted(values) == [1,2,3]
    return -1 if sum(values[i] > values[j] for i in range(3) for j in range(i+1,3)) % 2 else 1


def direct_quartic_table(values):
    """Separate transcription of the lambda/c multiplication rules."""
    a,b = [QQ(x) for x in values[:6]], [QQ(x) for x in values[6:]]
    pos={(1,1):0,(2,2):1,(3,3):2,(1,2):3,(2,1):3,(2,3):4,(3,2):4,(1,3):5,(3,1):5}
    lam=lambda i,j,k,l:a[pos[i,j]]*b[pos[k,l]]-b[pos[i,j]]*a[pos[k,l]]
    C=lambda i:{1:lam(2,3,1,1),2:-lam(1,3,2,2),3:lam(1,2,3,3)}[i]
    cache={}
    def c(i,j,k):
        key=(i,j,k)
        if key in cache:return cache[key]
        if k == 0:
            if i>j:i,j=j,i
            rest=[1,2,3];rest.remove(j)
            if i!=j:rest.remove(i)
            kk=rest[0];ans=sum(c(j,kk,r)*c(r,i,kk)-c(i,j,r)*c(r,kk,kk) for r in range(1,4))
        else:
            if j==k and i!=j:j,i=i,k
            if i==j and j==k:
                rest=[1,2,3];rest.remove(i);j,k=rest;ans=sign([i,j,k])*lam(i,k,i,j)+C(i)
            elif i==j:
                j=k;rest=[1,2,3];rest.remove(i);rest.remove(j);k=rest[0];ans=sign([i,j,k])*lam(i,i,i,k)
            elif i==k:
                rest=[1,2,3];rest.remove(i);rest.remove(j);k=rest[0];ans=sign([i,j,k])*lam(i,k,j,j)/2+C(j)/2
            else:ans=sign([i,j,k])*lam(j,j,i,i)
        cache[key]=ans;return ans
    table=[]
    for i in range(4):
        row=[]
        for j in range(4):
            if i==0: row.append(vector(QQ,[int(k==j) for k in range(4)]))
            elif j==0: row.append(vector(QQ,[int(k==i) for k in range(4)]))
            else: row.append(vector(QQ,[c(i,j,k) for k in range(4)]))
        table.append(row)
    assert all(entry in ZZ for row in table for product in row for entry in product)
    return [[vector(ZZ, product) for product in row] for row in table]


def compute():
    stored = read(TARGET); assert stored["status"] == "BOUNDED_SLICE_NO_ACCEPTED_CHARACTER"
    raw, nf, form = recovered_maximal_form()
    assert [str(x) for x in form] == stored["maximal_cubic_preflight"]["binary_cubic_descending"]
    checked = 0
    for row in stored["bounded_search"]["raw_pair_candidates"]:
        values = [ZZ(x) for x in row["pair_coefficients"]]; A,B=form_matrices(values)
        R=PolynomialRing(QQ,names=("x","y"));x,y=R.gens();det=4*(x*A-y*B).det()
        actual=[str(det.monomial_coefficient(x**(3-i)*y**i)) for i in range(4)]
        assert actual == stored["maximal_cubic_preflight"]["binary_cubic_descending"]
        content=0
        for value in values:content=gcd(content,abs(int(value)))
        assert content == 1
        table=direct_quartic_table(values)
        for i,j,k in itertools.product(range(4),repeat=3):
            left=sum((table[r][k]*table[i][j][r] for r in range(4)),vector(ZZ,4))
            right=sum((table[i][r]*table[j][k][r] for r in range(4)),vector(ZZ,4))
            assert left == right
        # Every retained witness is e2.  Its left-multiplication matrix has
        # zero determinant, so e2 is an explicit nonzero zero divisor.
        left_e2=matrix(ZZ,4,4,lambda out,col:table[2][col][out])
        assert left_e2.det() == 0 and row["algebra_probe"]["coordinates"] == [0,1,0]
        checked += 1
    return {"schema":"elliptic-curves.curve302-fixed-resolvent-orbits-v3-verification.v1","status":"PASS",
            "maximal_order_index_positive":str(nf[3]),"maximal_binary_cubic_discriminant":str(nf.disc()),
            "pairs_with_exact_fixed_resolvent":checked,"primitive_pairs":checked,"associative_quartic_orders":checked,
            "explicit_nonzero_zero_divisor_certificates":checked,"accepted_strict_characters":0,"sealed_cover_jobs":0,
            "bindings":{relative(path):digest(path) for path in (TARGET,STRICT,PROTOCOL,SOURCE,Path(__file__))},
            "boundary":"The verification proves only that this finite retained slice lies in the nonfield locus. It supplies no class-group or rank conclusion."}


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("build","check"));args=parser.parse_args();result=compute()
    if args.mode=="build":put_new(OUTPUT,result)
    else:assert result==read(OUTPUT)
    print("PASS independent maximal-order and zero-divisor replay",flush=True)
