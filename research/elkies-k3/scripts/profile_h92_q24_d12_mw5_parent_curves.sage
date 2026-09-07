#!/usr/bin/env sage -python
"""
Profile the five MW generators of the H3 q24 child D12 fibration as actual
curves in the old D13 parent.

Key geometric choice:
  F_new = D24.
The original q24 chamber gives two explicit old-fibre components with
F_new-degree one:
  R3 and D13_affine.
Choose R3 as the NEW zero section.  Then D13_affine is an explicit rational
section of the D12 fibration.

Build the child frame with this EFFECTIVE zero, root-adapt its D12 lattice,
LLL-reduce the MW5 quotient, and reconstruct an effective (-2)-section for
each of the five quotient basis vectors by exact closest-vector reduction in
D12.

For every generator report:
  * new MW vector and height;
  * D12 local correction and P.O_new;
  * its parent NS class;
  * degree on the OLD D13 fibre;
  * intersection with the OLD D13 zero;
  * whether it equals the explicit old affine component.

This tells us which MW5 curves are cheapest to realize in the existing D13
equation and therefore which should anchor the q24 -> D12 equation compiler.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix,
    lcm, matrix, pari, vector
)


Q24_WITNESS = vector(ZZ, (
    0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1,
))


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
        h/"projects"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def roots_and_data(gram):
    qf=pari(gram).qfminim(2)
    count=ZZ(qf[0])
    if not count:
        return (),matrix(ZZ,0,gram.nrows()),(0,0,1)
    half=[vector(ZZ,c) for c in matrix(ZZ,qf[2]).columns()]
    roots=tuple(half+[-r for r in half])
    basis=matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg=basis*gram*basis.transpose()
    return roots,basis,(basis.rank(),count,abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots,unused,data=roots_and_data(gram)
    rank=data[0]
    # Lexicographic positivity is deterministic and additive.
    positive=[
        r for r in roots
        if next(v for v in r if v)!=0 and next(v for v in r if v)>0
    ]
    pset={tuple(r) for r in positive}
    simple=[
        r for r in positive
        if not any(tuple(r-left) in pset for left in positive)
    ]
    M=matrix(ZZ,[list(r) for r in simple])
    assert M.nrows()==M.rank()==rank
    C=M*gram*M.transpose()
    assert all(C[i,i]==2 for i in range(rank))
    assert all(C[i,j] in (0,-1) for i in range(rank) for j in range(rank) if i!=j)
    return M,C


def child_frame_with_zero(ns,fibre,zero):
    assert fibre*ns*fibre==0
    assert zero*ns*zero==-2 and zero*ns*fibre==1
    mate=zero+fibre
    assert mate*ns*mate==0 and mate*ns*fibre==1
    orth=matrix(
        ZZ,[list(fibre*ns),list(mate*ns)]
    ).right_kernel_matrix()
    B=matrix(
        ZZ,[list(fibre),list(mate)]+[list(r) for r in orth.rows()]
    )
    assert abs(B.det())==1
    child=-(orth*ns*orth.transpose())
    assert B*ns*B.transpose()==block_diagonal_matrix(
        matrix(ZZ,((0,1),(1,0))),-child
    )
    return child,B


def root_adapt(child):
    roots,rb,data=roots_and_data(child)
    assert data==(12,264,4),data
    simple,cartan=deterministic_simple_roots(child)
    assert cartan.det()==4

    smith,left,right=rb.smith_form()
    assert smith==left*rb*right
    assert tuple(abs(smith[i,i]) for i in range(12))==(1,)*12

    completion=right.inverse()
    initial=simple.stack(completion[12:])
    assert abs(initial.det())==1

    G=initial*child*initial.transpose()
    root=G[:12,:12]
    coupling=G[:12,12:]
    tail=G[12:,12:]
    H=tail-coupling.transpose()*root.inverse()*coupling

    scale=lcm(QQ(v).denominator() for v in H.list())
    lll=matrix(ZZ,pari((scale*H).change_ring(ZZ)).qflllgram())
    assert abs(lll.det())==1
    change=block_diagonal_matrix(identity_matrix(ZZ,12),lll.transpose())
    A=change*initial
    G=A*child*A.transpose()

    root=G[:12,:12]
    coupling=G[:12,12:]
    tail=G[12:,12:]
    H=tail-coupling.transpose()*root.inverse()*coupling
    return A,G,H


def effective_section_for_z(G,H,z):
    z=vector(ZZ,z)
    root=G[:12,:12]
    base=vector(ZZ,[0]*12+list(z))
    pairing=vector(QQ,base*G[:,:12])
    dual=pairing*root.inverse()

    lattice=IntegralLattice(root)
    iterator=lattice.enumerate_close_vectors(-dual)

    # enumerate_close_vectors is ordered by distance.  The first integral lift
    # gives the minimal ADE correction in this quotient coset.
    shift=vector(ZZ,next(iterator))
    w=base+vector(ZZ,list(shift)+[0]*5)
    norm=ZZ(w*G*w)
    assert norm%2==0
    a=(norm-2)//2
    section=vector(ZZ,[a,1]+list(w))
    ns_child=block_diagonal_matrix(
        matrix(ZZ,((0,1),(1,0))),-G
    )
    assert section*ns_child*section==-2

    height=QQ(z*H*z)
    correction=QQ(norm)-height
    PdotO=ZZ(a-1)
    assert QQ(4+2*PdotO)-correction==height
    assert PdotO>=0
    return section,height,correction,PdotO,norm


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
FRAME=ROOT/"elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
LOCAL=ROOT/"artifacts/local/elkies-k3"
frame=load_gram(FRAME)
U2=matrix(ZZ,((0,1),(1,0)))
ns=block_diagonal_matrix(U2,-frame)

Fold=vector(ZZ,[1,0]+[0]*17)
Oold=vector(ZZ,[-1,1]+[0]*17)
Fnew=vector(ZZ,[12,2]+list(Q24_WITNESS))
assert Fnew*ns*Fnew==0 and Fnew*ns*Fold==2

# Old effective D13 simple components.
C=[]
for i in range(13):
    c=vector(ZZ,[0,0]+[-ZZ(j==i) for j in range(17)])
    assert c*ns*c==-2
    C.append(c)

# D13 affine.
root=frame[:13,:13]
qf=pari(root).qfminim(2)
half=matrix(ZZ,qf[2]).transpose().rows()
roots=[vector(ZZ,r) for r in half]+[-vector(ZZ,r) for r in half]
highest=max((r for r in roots if all(v>=0 for v in r)),key=lambda r:sum(r))
Caff=Fold+vector(ZZ,[0,0]+list(highest)+[0]*4)
assert Caff*ns*Caff==-2

pairings=[ZZ(Fnew*ns*c) for c in C]
assert pairings==[0,0,1,0,0,0,0,0,0,0,0,0,0]
assert Fnew*ns*Caff==1

# Choose the leaf R3 as effective new zero.
Onew=C[2]
assert Onew*ns*Fnew==1
assert Caff*ns*Fnew==1

child,Bzero=child_frame_with_zero(ns,Fnew,Onew)
A,G,H=root_adapt(child)
Badapt=block_diagonal_matrix(identity_matrix(ZZ,2),A)*Bzero
Gns=block_diagonal_matrix(U2,-G)
assert Badapt*ns*Badapt.transpose()==Gns

print(
    "Q24D12MW5_SETUP|"
    "new_zero=old_D13_R3|other_explicit_section=old_D13_affine|"
    f"root_data={roots_and_data(G)[2]}|"
    f"height_gram={';'.join(','.join(str(v) for v in row) for row in H.rows())}|"
    "status=PASS",
    flush=True,
)

# Coordinate/profile of explicit old affine section.
affcoords=vector(QQ,Caff)*Badapt.inverse()
assert all(v in ZZ for v in affcoords)
affcoords=vector(ZZ,affcoords)
assert affcoords[1]==1
zaff=vector(ZZ,affcoords[-5:])
haff=QQ(zaff*H*zaff)
POaff=ZZ(Caff*ns*Onew)
normaff=ZZ(vector(ZZ,affcoords[2:])*G*vector(ZZ,affcoords[2:]))
corraff=QQ(normaff)-haff
assert QQ(4+2*POaff)-corraff==haff

print(
    "Q24D12MW5_AFFINE|"
    f"mw={','.join(map(str,zaff))}|height={haff}|corr={corraff}|"
    f"PdotO={POaff}|old_degree={Caff*ns*Fold}|old_O={Caff*ns*Oold}|"
    "status=PASS_EXPLICIT_NEW_SECTION",
    flush=True,
)

profiles=[]
sections=[]
for i in range(5):
    z=vector(ZZ,[ZZ(j==i) for j in range(5)])
    section_child,height,corr,PO,norm=effective_section_for_z(G,H,z)
    parent=vector(ZZ,section_child)*Badapt
    assert parent*ns*parent==-2
    assert parent*ns*Fnew==1
    sections.append(parent)

    old_degree=ZZ(parent*ns*Fold)
    old_O=ZZ(parent*ns*Oold)
    equal_aff=bool(parent==Caff)
    equal_neg_aff=False  # class negation is not an effective section operation

    rec={
        "index":i+1,
        "mw":list(map(int,z)),
        "height":str(height),
        "D12_local_correction":str(corr),
        "P_dot_O_new":int(PO),
        "frame_norm":int(norm),
        "parent_class":list(map(int,parent)),
        "old_D13_fibre_degree":int(old_degree),
        "old_D13_zero_intersection":int(old_O),
        "equals_old_D13_affine":equal_aff,
    }
    profiles.append(rec)

    print(
        "Q24D12MW5_GEN|"
        f"i={i+1}|height={height}|corr={corr}|PdotO={PO}|"
        f"old_degree={old_degree}|old_O={old_O}|"
        f"equals_old_affine={int(equal_aff)}|"
        f"class={','.join(map(str,parent))}|status=PASS",
        flush=True,
    )

# Confirm MW independence and identify the explicit affine vector in the basis.
assert H.rank()==5
assert all(s*ns*Fnew==1 for s in sections)

# Search small MW words for the exact old affine section, since LLL sign/basis
# choices need not make it literally e1.
found_aff=None
for bound in range(1,5):
    from itertools import product
    for coeffs in product(range(-bound,bound+1),repeat=5):
        if not any(coeffs):
            continue
        z=vector(ZZ,coeffs)
        # only inspect words at least as short as a small explicit component.
        if QQ(z*H*z)>12:
            continue
        try:
            section_child,height,corr,PO,norm=effective_section_for_z(G,H,z)
        except Exception:
            continue
        parent=vector(ZZ,section_child)*Badapt
        if parent==Caff:
            found_aff=(
                list(map(int,z)),str(height),str(corr),int(PO)
            )
            break
    if found_aff:
        break
assert found_aff is not None

print(
    "Q24D12MW5_AFFINE_WORD|"
    f"mw={','.join(map(str,found_aff[0]))}|height={found_aff[1]}|"
    f"corr={found_aff[2]}|PdotO={found_aff[3]}|status=PASS",
    flush=True,
)

profiles_sorted=sorted(
    profiles,
    key=lambda r:(r["old_D13_fibre_degree"],r["P_dot_O_new"],r["index"])
)
print(
    "Q24D12MW5_CHEAPEST|"
    + "|".join(
        f"G{r['index']}:olddeg={r['old_D13_fibre_degree']},"
        f"h={r['height']},PO={r['P_dot_O_new']}"
        for r in profiles_sorted
    )
    + "|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-d12-mw5-parent-profiles.v1",
    "status":"PASS_EXACT_D12_MW5_PARENT_PROFILES",
    "new_fibre":list(map(int,Fnew)),
    "new_zero":{
        "curve":"old_D13_R3",
        "class":list(map(int,Onew)),
    },
    "explicit_other_section":{
        "curve":"old_D13_affine",
        "class":list(map(int,Caff)),
        "mw_coordinates":list(map(int,zaff)),
        "height":str(haff),
        "D12_local_correction":str(corraff),
        "P_dot_O":int(POaff),
        "basis_word":found_aff[0],
    },
    "MW_height_gram":[[str(v) for v in row] for row in H.rows()],
    "generators":profiles,
    "next":(
        "Use the lowest old-D13-degree new sections as equation-level anchors. "
        "In particular old_R3 and old_affine are already explicit curves, so "
        "the D12 model has an effective zero and at least one rational MW "
        "section before any new section search."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/"q24-d12-mw5-parent-profiles.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D12MW5_RESULT|rank=5|explicit_zero=old_R3|"
    "explicit_section=old_affine|status=PASS_EXACT_D12_MW5_PARENT_PROFILES",
    flush=True,
)
