#!/usr/bin/env python3
"""Uniform degree-six/eight relations at frozen addresses; no point search."""
import argparse
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import subprocess
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'LOW_DEGREE_TRIPLE_PANEL_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
BARRIER=r.OUT/'rank_jump_triple_degree_barrier_v1.json'
PAIRS=r.OUT/'rank_jump_degree_one_relation_panel_v1.json'
COSETS=r.OUT/'rank_jump_low_degree_triple_cosets_v1.json'
OUTPUT=r.OUT/'rank_jump_low_degree_triple_panel_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-low-degree-triple-panel-v1'


def triples(inp):
    return sorted({tuple(sorted(labels)) for b in inp['blocks'].values() for labels in combinations(b['labels'],3)})


def trace_sum(inp,labels):
    return [sum(inp['covers'][label]['trace'][i] for label in labels) for i in range(17)]


def bindings(paths):
    return {str(p.relative_to(r.ROOT)):hash_file(p) for p in paths}


def coset_worker():
    from sage.all import ZZ,matrix,vector,identity_matrix,pari
    inp=r.read(INPUT);barrier=r.read(BARRIER);assert barrier['status']=='PASS'
    G=matrix(ZZ,inp['gram']);assert G.det()==948 and G.is_positive_definite()
    assert int(pari(G).qfminim(3)[0])==0
    classes={}
    for labels in triples(inp):
        w=trace_sum(inp,labels);key=r.pack([x%2 for x in w]);assert key
        classes.setdefault(key,w)
    rows=[];folder=WORK/'cosets';folder.mkdir(parents=True,exist_ok=True)
    for key,ww in sorted(classes.items()):
        path=folder/f'{key}.json'
        if not path.exists():
            w=vector(ZZ,ww)
            H=(2*identity_matrix(ZZ,17)).stack(matrix(ZZ,[w])).row_module(ZZ).basis_matrix()
            U=matrix(ZZ,pari(H*G*H.transpose()).qflllgram());assert abs(U.det())==1
            V=U.transpose()*H;raw=pari(V*G*V.transpose()).qfminim(6)
            found=[]
            for v in matrix(ZZ,raw[2]).columns():
                z=v*V;norm=int(z*G*z);assert norm in (4,6)
                assert all((z[i]-w[i])%2==0 for i in range(17))
                zz=tuple(map(int,z));zz=min(zz,tuple(-x for x in zz))
                found.append({'vector':list(zz),'norm':norm})
            assert int(raw[0])==2*len(found) and len(found)<=1
            r.write_new(path,{'parity':key,'representative':ww,'reduced_basis':[[int(x) for x in row] for row in V],
                              'short_representatives_up_to_sign':found})
        rows.append(r.read(path))
    bykey={x['parity']:x for x in rows};rules=[]
    for labels in triples(inp):
        w=trace_sum(inp,labels);key=r.pack([x%2 for x in w]);found=bykey[key]['short_representatives_up_to_sign']
        rule={'labels':list(labels),'parity':key,'eligible':bool(found)}
        if found:
            z=found[0]['vector'];rule['total_intersection_degree']=found[0]['norm']+2
            rule['conjugate_generic_words']=[[(w[i]+sign*z[i])//2 for i in range(17)] for sign in (1,-1)]
        rules.append(rule)
    r.write_new(COSETS,{'schema':'rank-jump.low-degree-triple-cosets.v1','status':'PASS','cosets':rows,'triples':rules,
        'bindings':bindings((INPUT,BARRIER,PROTOCOL,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py'))})
    print('COSETS',len(rows),'TRIPLES',len(rules),'ELIGIBLE',sum(x['eligible'] for x in rules),flush=True)


def empty_row(key,block,pair):
    n=len(block['labels'])
    return {'block_key':key,'parameter':block['parameter'],'labels':block['labels'],'compatible_cover_count':n,
            'co_split_triples':0,'eligible_by_total_degree':{'6':0,'8':0},'incidences_by_total_degree':{'6':0,'8':0},
            'triple_checks':[],'canonical_signed_relations':[],'triple_relation_rank':0,
            'pair_relation_rank':pair['signed_relation_rank_over_Q'],
            'combined_relation_rank':pair['signed_relation_rank_over_Q'],
            'extra_relation_rank_beyond_pairs':0,
            'native_quotient_upper_bound':n-pair['signed_relation_rank_over_Q']}


def fibre_worker(index):
    from sage.all import QQ,EllipticCurve,matrix
    inp=r.read(INPUT);key,block=sorted(inp['blocks'].items())[index]
    pair=next(p for p in r.read(PAIRS)['rows'] if p['block_key']==key)
    rules={tuple(x['labels']):x for x in r.read(COSETS)['triples']}
    t=QQ(block['parameter']);labels=block['labels'];n=len(labels)
    def ev(cs):
        out=QQ(0)
        for c in reversed(cs):out=out*t+QQ(c)
        return out
    A,B=ev(inp['A']),ev(inp['B']);assert 4*A**3+27*B**2
    E=EllipticCurve(QQ,[A,B]);basis=[]
    for c in inp['sections']:
        x=ev(c['x_coefficients_low_to_high'])
        if 'y_coefficients_low_to_high' in c:y=ev(c['y_coefficients_low_to_high'])
        else:
            ch=c['chord'];P=basis[ch['reference_basis_index']]
            y=P[1]+ev(ch['slope_coefficients_low_to_high'])*(x-P[0])
        basis.append(E(x,y))
    @lru_cache(None)
    def generic(word):return sum((int(a)*P for a,P in zip(word,basis,strict=True)),E(0))
    points={}
    for label in labels:
        c=inp['covers'][label];q=ev(c['q']);assert q>0 and q.is_square();u=q.sqrt()
        x0,x1,y0,y1=[ev(c['lift'][s+'_coefficients']) for s in ('x0','x1','y0','y1')]
        pp=[E(x0+sign*u*x1,y0+sign*u*y1) for sign in (1,-1)]
        assert pp[0]!=pp[1] and pp[0]+pp[1]==generic(tuple(c['trace']))
        points[label]=pp
    @lru_cache(None)
    def pair_sum(a,sa,b,sb):return points[a][sa]+points[b][sb]
    row=empty_row(key,block,pair);triple_rows=[]
    for triple in combinations(sorted(labels),3):
        row['co_split_triples']+=1;rule=rules[triple]
        if not rule['eligible']:continue
        degree=rule['total_intersection_degree'];row['eligible_by_total_degree'][str(degree)]+=1
        matches=[]
        for wi,word in enumerate(rule['conjugate_generic_words']):
            S=generic(tuple(word))
            for branches in product(range(2),repeat=3):
                if pair_sum(triple[0],branches[0],triple[1],branches[1])+points[triple[2]][branches[2]]==S:
                    matches.append({'word_index':wi,'branch_indices':list(branches)})
        matchset={(m['word_index'],tuple(m['branch_indices'])) for m in matches}
        assert all((1-wi,tuple(1-x for x in bs)) in matchset for wi,bs in matchset)
        # Each conjugation orbit has exactly one word-index-zero labelled match.
        for m in matches:
            if m['word_index']!=0:continue
            signs=[(-1)**x for x in m['branch_indices']];word=rule['conjugate_generic_words'][0].copy()
            for label,sign in zip(triple,signs):
                if sign==-1:word=[a-b for a,b in zip(word,inp['covers'][label]['trace'])]
            indices=[labels.index(x) for x in triple];rr=[0]*n
            for i,s in zip(indices,signs):rr[i]=s
            triple_rows.append(rr)
            row['canonical_signed_relations'].append({'labels':list(triple),'indices':indices,'signs':signs,
                'total_intersection_degree':degree,'branch_indices':m['branch_indices'],
                'branch_sum_translate_word':rule['conjugate_generic_words'][0],
                'canonical_signed_relation_generic_word':word})
            row['incidences_by_total_degree'][str(degree)]+=1
        row['triple_checks'].append({'labels':list(triple),'total_intersection_degree':degree,'branch_matches':matches})
    pair_rows=[]
    for edge in pair['signed_quotient_edges']:
        rr=[0]*n
        for i,s in zip(edge['indices'],edge['signs']):rr[i]=s
        pair_rows.append(rr)
    def rank(rows):return int(matrix(QQ,rows).rank()) if rows else 0
    assert rank(pair_rows)==row['pair_relation_rank']
    row['triple_relation_rank']=rank(triple_rows)
    row['combined_relation_rank']=rank(pair_rows+triple_rows)
    row['extra_relation_rank_beyond_pairs']=row['combined_relation_rank']-row['pair_relation_rank']
    row['native_quotient_upper_bound']=n-row['combined_relation_rank']
    r.write_new(WORK/'fibres'/str(index)/'result.json',row)
    print(key,'eligible',row['eligible_by_total_degree'],'incidences',row['incidences_by_total_degree'],
          'extra rank',row['extra_relation_rank_beyond_pairs'],flush=True)


def supervise(command,folder,seconds):
    folder.mkdir(parents=True,exist_ok=True);log=folder/'worker.log';state=folder/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),*command],
                                 stdout=out,stderr=out,timeout=seconds)
                d={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:d={'status':'TIMEOUT'}
        r.write_new(state,d)
    return r.read(state),log.read_text()


def run():
    state,log=supervise(['cosets'],WORK/'coset_execution',180);print(state,log,flush=True)
    if state['status']!='COMPLETE':return
    inp=r.read(INPUT);pairs={x['block_key']:x for x in r.read(PAIRS)['rows']};rows=[]
    for i,(key,block) in enumerate(sorted(inp['blocks'].items())):
        if len(block['labels'])<3:rows.append(empty_row(key,block,pairs[key]));continue
        folder=WORK/'fibres'/str(i);state,log=supervise(['fibre','--index',str(i)],folder,60)
        print(state,log,flush=True)
        result=folder/'result.json';rows.append(r.read(result) if result.exists() else {'block_key':key,'status':'UNKNOWN','execution':state,'log':log})
    r.write_new(OUTPUT,{'schema':'rank-jump.low-degree-triple-panel.v1','rows':rows,
        'bindings':bindings((INPUT,COSETS,PAIRS,BARRIER,PROTOCOL,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py')),
        'boundary':r.read(PROTOCOL)['failure_semantics']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','cosets','fibre']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.mode=='run':run()
    elif a.mode=='cosets':coset_worker()
    else:fibre_worker(a.index)
