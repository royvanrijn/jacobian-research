#!/usr/bin/env python3
"""Exact LDL and a separate rational group law for the relation census."""
import argparse
from itertools import combinations
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,pari
import retrospective as r
from fibre_discrimination import hash_file
from verify_native_intersection_solubility import short_vectors,add

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
COSETS=r.OUT/'rank_jump_degree_one_relation_cosets_v1.json'
PANEL=r.OUT/'rank_jump_degree_one_relation_panel_v1.json'
TRIPLES=r.OUT/'rank_jump_trace_zero_triple_panel_v1.json'
OUTPUT=r.OUT/'rank_jump_degree_one_relation_panel_verification_v1.json'


def mul(A,n,P):
    if n<0:return mul(A,-n,(P[0],-P[1])) if P else None
    out=None
    while n:
        if n&1:out=add(A,out,P)
        n//=2
        if n:P=add(A,P,P)
    return out


def row_rank(rows,n):
    piv={}
    for row in rows:
        v=[QQ(x) for x in row]
        for k in sorted(piv):
            if v[k]:
                a=v[k];v=[v[j]-a*piv[k][j] for j in range(n)]
        inds=[i for i,x in enumerate(v) if x]
        if inds:
            k=inds[0];a=v[k];piv[k]=[x/a for x in v]
    return len(piv)


def verify():
    inp=r.read(INPUT);co=r.read(COSETS);panel=r.read(PANEL);triples=r.read(TRIPLES)
    for doc in (inp,co,panel,triples):
        for path,sha in doc['bindings'].items():assert hash_file(r.ROOT/path)==sha
    R=PolynomialRing(QQ,'t');G=matrix(ZZ,inp['gram'])
    assert G.det()==948 and G.is_positive_definite() and G[0,0]==4
    U=matrix(ZZ,pari(G).qflllgram());assert abs(U.det())==1
    low,nodes=short_vectors(U.transpose()*G*U,3);assert low==[]
    short_checks=[]
    for old in co['cosets']:
        V=matrix(ZZ,old['reduced_basis']);w=vector(ZZ,old['representative'])
        assert abs(V.det())==2**16
        for v in V.rows():assert all(x%2==0 for x in v) or all((v[i]-w[i])%2==0 for i in range(17))
        vs,count=short_vectors(V*G*V.transpose(),6)
        actual={tuple(vector(ZZ,v)*V) for v in vs};expected=set()
        for item in old['short_representatives_up_to_sign']:
            z=vector(ZZ,item['vector']);assert z*G*z==item['norm']
            expected.add(tuple(z));expected.add(tuple(-z))
        assert actual==expected
        assert old['has_norm_six']==any(x['norm']==6 for x in old['short_representatives_up_to_sign'])
        short_checks.append({'parity':old['parity'],'signed_count':len(actual),'nodes':count})
    cosets={x['parity']:x for x in co['cosets']}
    for pair in co['pairs']:
        w=[sum(inp['covers'][label]['trace'][i] for label in pair['labels']) for i in range(17)]
        parity=r.pack([x%2 for x in w]);assert parity==pair['parity']
        assert pair['has_degree_one_relation']==cosets[parity]['has_norm_six']
        if pair['has_degree_one_relation']:
            z=cosets[parity]['short_representatives_up_to_sign'][0]['vector']
            assert pair['conjugate_generic_words']==[[(w[i]+sign*z[i])//2 for i in range(17)] for sign in (1,-1)]
    A,B=R(inp['A']),R(inp['B']);cover_functions={}
    for label,c in inp['covers'].items():
        q=R(c['q']);x0,x1,y0,y1=[R(c['lift'][key+'_coefficients']) for key in ('x0','x1','y0','y1')]
        assert q.degree()==2 and q.is_squarefree()
        assert y0*y0+q*y1*y1==x0**3+3*q*x0*x1*x1+A*x0+B
        assert 2*y0*y1==3*x0*x0*x1+q*x1**3+A*x1
        assert vector(ZZ,c['trace'])*G*vector(ZZ,c['trace'])==10
        cover_functions[label]=(q,x0,x1,y0,y1)
    pair_rules={tuple(x['labels']):x for x in co['pairs']};rows=[];total_edges=0
    expected_keys=sorted(inp['blocks']);assert sorted(x['block_key'] for x in panel['rows'])==expected_keys
    for old in panel['rows']:
        block=inp['blocks'][old['block_key']];labels=block['labels'];n=len(labels)
        assert old['compatible_cover_count']==n and old['labels']==labels
        t=QQ(block['parameter']);aa,bb=A(t),B(t);assert 4*aa**3+27*bb**2
        basis=[]
        for c in inp['sections']:
            x=R(c['x_coefficients_low_to_high'])(t)
            if 'y_coefficients_low_to_high' in c:y=R(c['y_coefficients_low_to_high'])(t)
            else:
                ch=c['chord'];P=basis[ch['reference_basis_index']]
                y=P[1]+R(ch['slope_coefficients_low_to_high'])(t)*(x-P[0])
            assert y*y==x**3+aa*x+bb;basis.append((x,y))
        generic_cache={}
        def generic(word):
            key=tuple(word)
            if key not in generic_cache:
                out=None
                for n,P in zip(word,basis,strict=True):out=add(aa,out,mul(aa,int(n),P))
                generic_cache[key]=out
            return generic_cache[key]
        points={}
        for label in labels:
            q,x0,x1,y0,y1=cover_functions[label];v=q(t);assert v>0 and v.is_square();u=v.sqrt()
            points[label]=[(x0(t)+sign*u*x1(t),y0(t)+sign*u*y1(t)) for sign in (1,-1)]
            assert points[label][0]!=points[label][1]
            assert add(aa,*points[label])==generic(inp['covers'][label]['trace'])
        checks=[];edges=[];relation_rows=[];eligible=0
        for label_pair in combinations(labels,2):
            pair=tuple(sorted(label_pair));rule=pair_rules[pair]
            if not rule['has_degree_one_relation']:continue
            eligible+=1;matches=[]
            for wi,word in enumerate(rule['conjugate_generic_words']):
                S=generic(word)
                for ai in range(2):
                    for bi in range(2):
                        if add(aa,points[pair[0]][ai],points[pair[1]][bi])==S:
                            matches.append({'word_index':wi,'branch_indices':[ai,bi]})
            checks.append({'labels':list(pair),'branch_matches':matches})
            assert len(matches) in (0,2)
            if matches:
                m=next(x for x in matches if x['word_index']==0);i,j=[labels.index(l) for l in pair]
                signs=[(-1)**a for a in m['branch_indices']]
                edges.append({'indices':[i,j],'labels':list(pair),'signs':signs,'generic_word':rule['conjugate_generic_words'][0]})
                v=[0]*n;v[i],v[j]=signs;relation_rows.append(v)
        assert checks==old['pair_checks'] and edges==old['signed_quotient_edges']
        assert eligible==old['globally_degree_one_pairs_among_compatible_covers']
        rank=row_rank(relation_rows,n)
        assert rank==old['signed_relation_rank_over_Q']
        assert old['native_quotient_rank_upper_bound_from_pair_relations']==n-rank
        assert len(edges)==old['rational_degree_one_components_at_this_parameter_mod_conjugation']
        rows.append({'block_key':old['block_key'],'covers':n,'eligible_degree_one_pairs':eligible,
                     'incidences':len(edges),'signed_relation_rank':rank,'native_quotient_upper_bound':n-rank})
        total_edges+=len(edges)
    triple_rows={x['block_key']:x for x in triples['rows']};triple_count=0
    for key,block in inp['blocks'].items():
        expected=[];count=0
        for labels in combinations(block['labels'],3):
            count+=1
            sums=[sum(inp['covers'][label]['trace'][i] for label in labels) for i in range(17)]
            if all(s%2==0 for s in sums):expected.append({'labels':list(labels),'generic_word':[s//2 for s in sums]})
        assert triple_rows[key]['trace_even_triples']==expected
        assert triple_rows[key]['co_split_triples']==count;triple_count+=count
    assert triple_count==triples['triple_tests']==2853 and triples['trace_even_count']==0
    assert total_edges==18 and sum(x['has_degree_one_relation'] for x in co['pairs'])==131
    return {'schema':'rank-jump.degree-one-relation-panel-verification.v1','status':'PASS','rows':rows,
            'coset_checks':short_checks,'global_minimum_enumeration_nodes':nodes,
            'generic_lift_identities_checked':len(cover_functions),'total_degree_one_incidences':total_edges,
            'eligible_degree_one_pairs':131,'trace_even_triple_count':0,'triple_parities_checked':triple_count,
            'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,COSETS,PANEL,TRIPLES,Path(__file__),HERE/'verify_native_intersection_solubility.py',HERE/'fibre_discrimination.py',HERE/'retrospective.py')},
            'boundary':'Independent exact short-coset enumeration and specialized group arithmetic. Generic trace-class identification is inherited from the pinned atlas; lift polynomial identities are replayed here. Degree-one incidences constrain quotient directions but are not full rank certificates.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=verify()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS 519 cosets; 118 generic lift identities; 131 degree-one pairs; 18 incidences; 2853 triple parities with no zero class',flush=True)
