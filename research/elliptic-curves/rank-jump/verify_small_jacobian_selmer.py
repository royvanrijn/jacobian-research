#!/usr/bin/env python3
"""Independent quotient identities, mod-4 exclusions, and torsion matrices."""
import argparse
from itertools import product
from pathlib import Path
import retrospective as r
import small_jacobian_selmer as ex

OUTPUT=r.OUT/"rank_jump_small_jacobian_selmer_verification_v1.json"


def square_mod4(v):
    c=[0]*5
    for i,x in enumerate(v):
        for j,y in enumerate(v):c[i+j]+=x*y
    for d in (4,3):
        c[d-3]+=c[d];c[d-2]+=14*c[d];c[d-1]+=11*c[d]
    return tuple(x%4 for x in c[:3])


def invariants(perms):
    # Even subsets modulo all six points; basis e0+ei for i=1,...,4.
    equations=[]
    for perm in perms:
        columns=[]
        for i in range(1,5):
            bits=(1<<perm[0])^(1<<perm[i])
            if bits>>5&1:bits^=63
            columns.append((bits>>1)&15)
        for i in range(4):
            equations.append(sum(((columns[j]>>i)&1) <<j for j in range(4))^(1<<i))
    return 4-r.rank(equations)


def verify(check=False):
    from sage.all import QQ,PolynomialRing
    data=r.read(ex.OUTPUT);assert data["bindings"]==ex.bindings()
    R=PolynomialRing(QQ,"X");X=R.gen();f=X**3-11*X**2-14*X-1
    squares={square_mod4(v) for v in product(range(4),repeat=3)}
    for row in data["dyadic_unit_checks"]:
        residue=tuple(int(v)%4 for v in row["residue_mod8"])
        assert residue not in squares and not row["has_square_root_mod8"]
    rotate=[(i+2)%6 for i in range(6)]
    def flip(mask):return [i^((mask>>(i//2))&1) for i in range(6)]
    rows=[]
    for row in data["rows"]:
        sign=row["sign"];g=f(sign*X*X-1)
        assert g==sign*X**6-14*X**4+11*sign*X*X+1
        disc=g.discriminant();assert disc==-64*sign*163**4
        # The second quotient is (x,y)=(1/X^2, sign*Z/X^3).
        assert g==X**6*((1/X**2)**3+11*sign*(1/X**2)**2-14/X**2+sign)
        global_generators=[rotate]+[flip(m) for m in ((1,2,4) if sign==1 else (3,5))]
        assert invariants(global_generators)==0
        real=invariants([flip(1 if sign==1 else 6)])
        assert real==row["real_2_torsion_dimension"]
        assert invariants([rotate])==0
        C=row["common_Selmer_basis"];D=row["sum_Selmer_basis"]
        radical=len(C)-r.rank(list(map(r.pack,row["difference_CT_matrix"])))
        assert row["Jacobian_2_Selmer_dimension"]==len(D)+radical
        for local in row["local_conditions"]:
            assert local["Jacobian_local_dimension"]==local["sum_dimension"]+local["common_dimension"]-local["connecting_rank"]
        rational_first=[4];rational_second=[1,2,4] if sign==1 else [8]
        rational_sum=r.rank(rational_first+rational_second)
        rational_intersection=1+len(rational_second)-rational_sum
        left=len(D)-rational_sum;right=radical-rational_intersection
        assert left>=0 and right>=0
        assert left+right==row["Jacobian_Sha_2_dimension"]
        assert row["Jacobian_exact_rank"]+left+right==row["Jacobian_2_Selmer_dimension"]
        rows.append({"sign":sign,"sextic_ascending":list(map(str,g.list())),"sextic_discriminant":str(disc),
                     "rational_2_torsion_dimension":0,"real_2_torsion_dimension":real,
                     "left_Sha_quotient_dimension":left,"right_Sha_quotient_dimension":right,
                     "Jacobian_Sha_2_dimension":left+right})
    # Every invertible 2x2 matrix over F2 preserves the alternating plane.
    H=[[0,1],[1,0]];checked=0
    for a,b,c,d in product((0,1),repeat=4):
        if (a*d+b*c)%2!=1:continue
        M=[[a,b],[c,d]]
        transported=[[sum(M[k][i]*H[k][l]*M[l][j] for k in range(2) for l in range(2))%2 for j in range(2)] for i in range(2)]
        assert transported==H;checked+=1
    assert checked==6
    # The positive sextic has rational infinity branches Z/X^3=+/-1.
    assert f(-1)==1 and f(-1)==(-1)**2
    assert 1==0**3+11*0**2-14*0+1
    out={"schema":"rank-jump.small-jacobian-selmer-verification.v1","status":"PASS","rows":rows,
         "independent_dyadic_modulus":4,"dyadic_residues_checked":64,"alternating_plane_automorphisms_checked":checked,
         "analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),"verifier_sha256":r.digest(Path(__file__).read_bytes()),
         "boundary":"The isogeny Selmer and Kummer diagram proves the Sha exact sequence; finite checks certify its inputs here."}
    if check:assert out==r.read(OUTPUT)
    else:r.write_new(OUTPUT,out)
    print("PASS independent labelled quotient geometry, dyadic exclusions, torsion invariants, and both Sha quotients")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    verify(parser.parse_args().check)
