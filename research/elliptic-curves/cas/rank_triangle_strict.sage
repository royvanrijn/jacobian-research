#!/usr/bin/env sage-python
"""Exact Kummer/strict/ideal dictionary in the fixed accessibility complement."""
import json,gzip,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,AA,PolynomialRing,matrix,vector,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'rank_triangle_v1'
sys.path.insert(0,str(Path(__file__).parent))
from research_runtime.local_kummer import LocalSquareclasses
import icarm_curve302 as curve
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(M):return [list(map(int,row)) for row in M.rows()]
def main():
    paths=[ART/'curve302_recovered_quotient_local_filtration_v1.json',ART/'curve302_descent_anatomy_v1.json',
           ART/'curve302_descent_remaining_class_v1.json',ART/'rank_accessibility_subsets_v1/inputs.json.gz',Path(__file__)]
    loc,anatomy,remaining,data=[read(p) for p in paths[:4]];hp=anatomy['half_ideal_packet'];R=PolynomialRing(QQ,'z');f=R(hp['cubic_ascending'])
    gammas=[R([QQ(p['a']),-QQ(p['d'])**2]) for p in hp['points']]
    assert len(gammas)==31
    for p,P in zip(hp['points'],curve.POINTS):
        a,b,d=map(QQ,[p['a'],p['b'],p['d']]);x,y=map(QQ,P)
        assert a/d**2==4*x and b/d**3==8*y+4*x+4
        assert f(a/d**2)*d**6==b*b
    good=[[] for _ in gammas]
    for block in loc['kummer_rank_certificate']['complete_split_character_blocks']:
        p=block['prime'];roots=block['roots'];assert len(roots)==3
        for root in roots:assert f.change_ring(GF(p))(root)==0
        for i,g in enumerate(gammas):
            values=[g.change_ring(GF(p))(root) for root in roots];assert all(values)
            good[i].extend(int(not v.is_square()) for v in values)
        assert matrix(GF(2),good).rank()==block['rank_after_block']
    assert matrix(GF(2),good).rank()==31
    S=[r['place'] for r in loc['local_places'] if isinstance(r['place'],int)]
    nf=pari.nfinit([pari(f),S]);assert pari.nfcertify(nf)==[]
    betas=[pari.Mod(pari(g),pari(f)) for g in gammas];signatures=[[] for g in gammas]
    for p in S:
        L=LocalSquareclasses(nf,p)
        for row,g in zip(signatures,betas):row.extend(map(int,L.signature(g)))
    roots=f.roots(AA,multiplicities=False)
    for row,g in zip(signatures,gammas):row.extend(int(g(r)<0) for r in roots)
    L=matrix(GF(2),signatures);V=L.left_kernel();assert L.rank()==21 and V.dimension()==10
    J=matrix(GF(2),hp['words']);assert V==J.row_space()
    assert hp['words']==loc['local_filtration_mod_2']['strict_kernel_public_words']
    selected=next(r for r in data['curves'] if r['id']=='302');B=matrix(ZZ,selected['basis_words_in_D'])
    assert abs(B.det())==1
    B2=B.change_ring(GF(2));C=J*B2.inverse();E=C[:,17:];assert E.rank()==10
    generic_local=B2[:17,:]*L;assert generic_local.rank()==17
    quotient_checks=E.right_kernel().basis_matrix();assert quotient_checks.nrows()==4
    entries=[]
    for i in range(10):
        entries.append({'ideal_index':i,'word_in_public_D':hp['words'][i],
                        'generic_correction_mod2':list(map(int,C.row(i)[:17])),
                        'exceptional_word_mod2':list(map(int,C.row(i)[17:])),
                        'half_ideal_HNF':hp['half_ideals'][i],
                        'beta':'product of gamma_i for the displayed public word; (beta)=J_i^2'})
    assert remaining['ordinary_half_ideal_rank_interval']==[10,10]
    pure=[]
    for j in range(14):
        e=vector(GF(2),[int(k==j) for k in range(14)])
        q=quotient_checks*e;coeff=None;correction=None
        if not q:
            coeff=E.solve_left(e);correction=coeff*C[:,:17]
            assert (vector(GF(2),list(correction)+list(e))*B2*L).is_zero()
        pure.append({'target':j+1,'quotient_local_class':list(map(int,q)),
                     'strict_after_generic_correction':bool(not q),
                     'J_coordinates':list(map(int,coeff)) if coeff is not None else None,
                     'generic_correction':list(map(int,correction)) if correction is not None else None})
    out={'schema':'rank-triangle.strict-dictionary.v1','status':'PASS_EXACT_KUMMER_DICTIONARY',
         'bindings':{str(p):sha(p) for p in paths},'cubic':hp['cubic_ascending'],
         'kummer_generators':[list(map(str,g.list())) for g in gammas],
         'good_character_matrix':rows(matrix(GF(2),good)),'local_signature_matrix':rows(L),
         'known_Kummer_dimension':31,'generic_dimension':17,'strict_dimension':10,
         'strict_public_basis':rows(J),'strict_basis_in_M17_E':rows(C),'quotient_local_check_matrix':rows(quotient_checks),
         'strict_to_half_ideals':entries,'individual_exceptional_classes':pure,
         'ordinary_half_ideal_image_dimension':10,'strict_unit_kernel_dimension':0,
         'basis_independent_description':'V=kernel of full S-localization restricted to the known rational Kummer image W; V injects into W/G because G intersects V trivially. Coordinates, individual J_i, and a complement splitting are not canonical.',
         'full_Selmer':'dimension=21+c_S; c_S>=10 remains UNKNOWN above',
         'boundary':'The whole strict Selmer kernel is not asserted to equal this known ten-dimensional V. Ideal-square/Artin independence is imported from the sealed anatomy certificates; local signatures and good-character independence are freshly recomputed.'}
    path=OUT/'strict.json'
    if path.exists():assert read(path)==out
    else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('STRICT PASS: W31, G17, V10, local quotient4; exact J dictionary',flush=True)
if __name__=='__main__':main()
