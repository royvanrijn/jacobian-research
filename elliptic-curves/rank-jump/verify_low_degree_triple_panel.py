#!/usr/bin/env python3
"""Independent exact LDL enumeration and rational group-law replay."""
import argparse
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector
import retrospective as r
from fibre_discrimination import hash_file
from verify_native_intersection_solubility import short_vectors,add
from verify_degree_one_relation_panel import mul,row_rank

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
COSETS=r.OUT/'rank_jump_low_degree_triple_cosets_v1.json'
PANEL=r.OUT/'rank_jump_low_degree_triple_panel_v1.json'
PAIRS=r.OUT/'rank_jump_degree_one_relation_panel_v1.json'
PAIR_VERIFIED=r.OUT/'rank_jump_degree_one_relation_panel_verification_v1.json'
BARRIER=r.OUT/'rank_jump_triple_degree_barrier_v1.json'
OUTPUT=r.OUT/'rank_jump_low_degree_triple_panel_verification_v1.json'


def verify():
    inp=r.read(INPUT);co=r.read(COSETS);panel=r.read(PANEL);pairs=r.read(PAIRS)
    pv=r.read(PAIR_VERIFIED);barrier=r.read(BARRIER)
    assert co['status']==pv['status']==barrier['status']=='PASS'
    for doc in (inp,co,panel,pairs,pv,barrier):
        for path,sha in doc['bindings'].items():assert hash_file(r.ROOT/path)==sha
    G=matrix(ZZ,inp['gram']);assert G.det()==948 and G.is_positive_definite()
    checks=[];cosets={};triple_set=set();expected_keys=set()
    for block in inp['blocks'].values():
        for labels in combinations(sorted(block['labels']),3):
            triple_set.add(labels)
            w=[sum(inp['covers'][label]['trace'][i] for label in labels) for i in range(17)]
            parity=r.pack([x%2 for x in w]);assert parity;expected_keys.add(parity)
    assert {x['parity'] for x in co['cosets']}==expected_keys
    assert len(co['cosets'])==len(expected_keys)
    for old in co['cosets']:
        V=matrix(ZZ,old['reduced_basis']);w=vector(ZZ,old['representative'])
        assert abs(V.det())==2**16 and r.pack([int(x)%2 for x in w])==old['parity']
        for v in V.rows():assert all(x%2==0 for x in v) or all((v[i]-w[i])%2==0 for i in range(17))
        vs,nodes=short_vectors(V*G*V.transpose(),6)
        actual={tuple(vector(ZZ,v)*V) for v in vs};expected=set()
        for item in old['short_representatives_up_to_sign']:
            z=vector(ZZ,item['vector']);assert z*G*z==item['norm'] and item['norm'] in (4,6)
            assert all((z[i]-w[i])%2==0 for i in range(17))
            expected.add(tuple(z));expected.add(tuple(-z))
        assert actual==expected and len(actual)<=2
        cosets[old['parity']]=old
        checks.append({'parity':old['parity'],'signed_short_count':len(actual),'enumeration_nodes':nodes})
        if len(checks)%500==0:print('verified cosets',len(checks),flush=True)
    rules={}
    assert len(co['triples'])==len(triple_set)
    for rule in co['triples']:
        labels=tuple(rule['labels']);assert labels in triple_set and labels not in rules
        w=[sum(inp['covers'][label]['trace'][i] for label in labels) for i in range(17)]
        parity=r.pack([x%2 for x in w]);assert rule['parity']==parity
        found=cosets[parity]['short_representatives_up_to_sign'];assert rule['eligible']==bool(found)
        if found:
            z=found[0]['vector'];assert rule['total_intersection_degree']==found[0]['norm']+2
            assert rule['conjugate_generic_words']==[[(w[i]+s*z[i])//2 for i in range(17)] for s in (1,-1)]
        rules[labels]=rule
    R=PolynomialRing(QQ,'t');A=R(inp['A']);B=R(inp['B'])
    ps={x['block_key']:x for x in pairs['rows']};rows=[]
    assert sorted(x['block_key'] for x in panel['rows'])==sorted(inp['blocks'])
    for old in panel['rows']:
        block=inp['blocks'][old['block_key']];labels=block['labels'];n=len(labels)
        assert old['parameter']==block['parameter'] and old['labels']==labels and old['compatible_cover_count']==n
        tr_rows=[];pair_rows=[];eligible={'6':0,'8':0};incidences={'6':0,'8':0};tr_checks=[];relations=[];count=0
        for edge in ps[old['block_key']]['signed_quotient_edges']:
            rr=[0]*n
            for i,s in zip(edge['indices'],edge['signs']):rr[i]=s
            pair_rows.append(rr)
        if n>=3:
            t=QQ(block['parameter']);aa,bb=A(t),B(t);assert 4*aa**3+27*bb**2
            basis=[]
            for c in inp['sections']:
                x=R(c['x_coefficients_low_to_high'])(t)
                if 'y_coefficients_low_to_high' in c:y=R(c['y_coefficients_low_to_high'])(t)
                else:
                    ch=c['chord'];P=basis[ch['reference_basis_index']]
                    y=P[1]+R(ch['slope_coefficients_low_to_high'])(t)*(x-P[0])
                assert y*y==x**3+aa*x+bb;basis.append((x,y))
            @lru_cache(None)
            def generic(word):
                out=None
                for a,P in zip(word,basis,strict=True):out=add(aa,out,mul(aa,int(a),P))
                return out
            points={}
            for label in labels:
                c=inp['covers'][label];q=R(c['q'])(t);assert q>0 and q.is_square();u=q.sqrt()
                x0,x1,y0,y1=[R(c['lift'][s+'_coefficients'])(t) for s in ('x0','x1','y0','y1')]
                points[label]=[(x0+s*u*x1,y0+s*u*y1) for s in (1,-1)]
                for x,y in points[label]:assert y*y==x**3+aa*x+bb
                assert points[label][0]!=points[label][1]
                assert add(aa,*points[label])==generic(tuple(c['trace']))
            @lru_cache(None)
            def pair_sum(a,sa,b,sb):return add(aa,points[a][sa],points[b][sb])
            for triple in combinations(sorted(labels),3):
                count+=1;rule=rules[triple]
                if not rule['eligible']:continue
                degree=rule['total_intersection_degree'];eligible[str(degree)]+=1;matches=[]
                for wi,word in enumerate(rule['conjugate_generic_words']):
                    S=generic(tuple(word))
                    for bs in product(range(2),repeat=3):
                        if add(aa,pair_sum(triple[0],bs[0],triple[1],bs[1]),points[triple[2]][bs[2]])==S:
                            matches.append({'word_index':wi,'branch_indices':list(bs)})
                ms={(m['word_index'],tuple(m['branch_indices'])) for m in matches}
                assert all((1-wi,tuple(1-x for x in bs)) in ms for wi,bs in ms)
                tr_checks.append({'labels':list(triple),'total_intersection_degree':degree,'branch_matches':matches})
                for m in matches:
                    if m['word_index']:continue
                    signs=[1-2*x for x in m['branch_indices']];word=rule['conjugate_generic_words'][0].copy()
                    for label,sign in zip(triple,signs):
                        if sign<0:word=[a-b for a,b in zip(word,inp['covers'][label]['trace'])]
                    signed_sum=None
                    for label,sign in zip(triple,signs):signed_sum=add(aa,signed_sum,mul(aa,sign,points[label][0]))
                    assert signed_sum==generic(tuple(word))
                    indices=[labels.index(label) for label in triple];rr=[0]*n
                    for i,s in zip(indices,signs):rr[i]=s
                    tr_rows.append(rr);incidences[str(degree)]+=1
                    relations.append({'labels':list(triple),'indices':indices,'signs':signs,'total_intersection_degree':degree,
                        'branch_indices':m['branch_indices'],'branch_sum_translate_word':rule['conjugate_generic_words'][0],
                        'canonical_signed_relation_generic_word':word})
        pr=row_rank(pair_rows,n);tr=row_rank(tr_rows,n);cr=row_rank(pair_rows+tr_rows,n)
        assert old['co_split_triples']==count and old['eligible_by_total_degree']==eligible
        assert old['incidences_by_total_degree']==incidences and old['triple_checks']==tr_checks
        assert old['canonical_signed_relations']==relations
        assert (old['pair_relation_rank'],old['triple_relation_rank'],old['combined_relation_rank'])==(pr,tr,cr)
        assert old['extra_relation_rank_beyond_pairs']==cr-pr and old['native_quotient_upper_bound']==n-cr
        rows.append({'block_key':old['block_key'],'co_split_triples':count,'eligible_by_total_degree':eligible,
            'incidences_by_total_degree':incidences,'pair_relation_rank':pr,'triple_relation_rank':tr,
            'combined_relation_rank':cr,'extra_relation_rank_beyond_pairs':cr-pr,'native_quotient_upper_bound':n-cr})
        if n>=3:print('verified',old['block_key'],eligible,incidences,'extra relation rank',cr-pr,flush=True)
    assert sum(x['co_split_triples'] for x in rows)==2853
    return {'schema':'rank-jump.low-degree-triple-panel-verification.v1','status':'PASS','rows':rows,'coset_checks':checks,
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,COSETS,PANEL,PAIRS,PAIR_VERIFIED,BARRIER,
            Path(__file__),HERE/'verify_native_intersection_solubility.py',HERE/'verify_degree_one_relation_panel.py',
            HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
        'boundary':'Independent complete norm<=6 enumeration and all specialized branch equations replayed. Degree formula and generic minimum are inherited from the pinned verified barrier. Relation constraints bound only the native quotient; no full rank or new rational parameter is asserted.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=verify()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS',len(d['coset_checks']),'cosets;',len(d['rows']),'addresses',flush=True)
