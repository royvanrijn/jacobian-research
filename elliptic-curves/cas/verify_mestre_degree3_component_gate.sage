#!/usr/bin/env sage-python
"""Exact split-I4 component gate on six certified determinant468 parents.

Checks the rational affine-A3 configuration and all20 nonnegative degree3
intersection profiles. A triangle supplies A2, and an untouched old I4
component supplies an orthogonal A1. The argument applies to every old
section triangle, not a bounded word dictionary. One worker,30 seconds.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
from sage.all import ZZ,matrix,vector,block_diagonal_matrix

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
SOURCE=ART/'mestre_rational_ns_gram_v2.json'
OUT=ART/'mestre_degree3_component_gate_v1.json'

def rows(m):return [list(map(int,r)) for r in m.rows()]

def build():
    source=json.loads(SOURCE.read_text()); records=[]
    profiles=[]
    for a in range(4):
        for b in range(4-a):
            for c in range(4-a-b):
                profile=[a,b,c,3-a-b-c]
                zero=next(i for i,v in enumerate(profile) if v==0)
                profiles.append(dict(intersections=profile,untouched_component=zero))
    assert len(profiles)==20 and len({tuple(r['intersections']) for r in profiles})==20
    triangle=matrix(ZZ,[[-2,1,1],[1,-2,1],[1,1,-2]])
    extended=block_diagonal_matrix(triangle,matrix(ZZ,[[-2]]))
    assert extended.rank()==3 and extended*vector(ZZ,[1,1,1,0])==0
    assert extended[:2,:2].det()==3 and abs(extended.matrix_from_rows_and_columns([0,1,3],[0,1,3]).det())==6
    for r in source['rows']:
        G=matrix(ZZ,r['rational_NS_Gram']); assert G.nrows()==18 and G.det()==-468
        assert r['basis'][:7]==['F','O','I2_minus','I2_plus','I4_1','I4_2','I4_3']
        e=matrix.identity(ZZ,18); F=e[0]
        components=matrix(ZZ,[F-e[4]-e[5]-e[6],e[4],e[5],e[6]])
        affine=components*G*components.transpose()
        assert affine==matrix(ZZ,[[-2,1,0,1],[1,-2,1,0],[0,1,-2,1],[1,0,1,-2]])
        assert sum(components.rows(),vector(ZZ,18))==F
        assert all(c*G*F==0 for c in components.rows())
        assert G[:2,:2]==matrix(ZZ,[[0,1],[1,-2]])
        records.append(dict(outer_u=r['outer_u'],old_fibre=list(map(int,F)),I4_component_classes=rows(components),
                            I4_intersection_matrix=rows(affine),arithmetic_NS_rank=18,
                            degree3_MW_rank_upper_bound=15,old_section_triangle_MW_rank_upper_bound=13))
    assert len(records)==6
    return dict(schema='mestre.degree3-component-gate.v1',status='ALL_OLD_SECTION_TRIANGLES_HAVE_MW_AT_MOST13_AND_DEGREE3_HAS_MW_AT_MOST15',
                input_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SOURCE]},
                degree3_component_profiles=profiles,triangle_and_untouched_component_gram=rows(extended),rows=records,
                proof='For a nef new fibre D with D.F=3, its nonnegative integral intersections with the four rational old I4 components sum to3, so one component C is vertical. C cannot be a multiple of D because F.C=0 and F.D=3. Hence at least one rational fibre-root direction and MW<=18-2-1=15. If D is the sum of three old rational sections meeting pairwise once, C is disjoint from each, and their Gram with C has rank3 and radical D. The vertical root space contains A2+A1, so MW<=18-2-3=13. Rationality and effectivity of the old I4 components and the full arithmetic NS rank are inputs from the existing full-NS proof.',
                boundary='The six certified determinant468 surfaces only. The degree3 bound allows MW14 and MW15 for nontriangle divisors; old degree>=4 and other surfaces remain open. No302 inverse or full MW coordinate basis is constructed.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args()
    signal.alarm(30);result=build()
    if args.build:
        assert not OUT.exists();OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    assert result==json.loads(OUT.read_text())
    print(result['status'],flush=True)
