#!/usr/bin/env python3
"""Family-wide degree-one rule, evaluated only on fixed retrospective fibres."""
import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
import subprocess
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'DEGREE_ONE_RELATION_PANEL_PROTOCOL.json'
CENSUS=r.OUT/'rank_jump_fibre_discrimination_v1.json'
BASE=r.OUT/'rank_jump_native_pair_collapse_locus_inputs_v1.json'
LATTICE=r.OUT/'rank_jump_norm_six_carrier_solubility_inputs_v1.json'
ATLAS=r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
COSETS=r.OUT/'rank_jump_degree_one_relation_cosets_v1.json'
OUTPUT=r.OUT/'rank_jump_degree_one_relation_panel_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-degree-one-relation-panel-v1'


def capture():
    source=r.read(CENSUS)
    blocks={k:b for k,b in source['blocks'].items() if b['dictionary']=='published-R17'}
    labels={c['label'] for b in blocks.values() for c in b['compatible_forms']}
    pairs=sorted({tuple(sorted((a['label'],b['label']))) for block in blocks.values()
                  for a,b in combinations(block['compatible_forms'],2)})
    assert len(blocks)==165 and len(labels)==118 and len(pairs)==519
    assert hash_file(ATLAS)=='78e037dc4170955b8f79ddce4d1d3e0c0d3e9bb8f9614644c59ccc7d605226c4'
    covers={};buf=None
    with ATLAS.open() as f:
        for line in f:
            if line=='    {\n':buf=[]
            if buf is not None:buf.append(line)
            if buf is not None and line in ('    },\n','    }\n'):
                label=next((s.split('"')[3] for s in buf if s.startswith('      "label": ')),None)
                if label in labels:
                    c=json.loads(''.join(buf).rstrip().rstrip(','))
                    covers[label]={'trace':c['published_basis_w'],'q':c['residual_chord']['q_coefficients'],'lift':c['lifted_section']}
                buf=None
    assert set(covers)==labels
    base=r.read(BASE)
    data={'schema':'rank-jump.degree-one-relation-panel-inputs.v1','covers':covers,'pairs':[list(p) for p in pairs],
          'blocks':{k:{'parameter':b['native_parameter'],'labels':[x['label'] for x in b['compatible_forms']]} for k,b in blocks.items()},
          'gram':r.read(LATTICE)['gram'],**{k:base[k] for k in ('A','B','sections')},
          'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (PROTOCOL,CENSUS,BASE,LATTICE,ATLAS)},
          'boundary':r.read(PROTOCOL)['failure_semantics']}
    r.write_new(INPUT,data)
    print('CAPTURED 165 addresses, 118 generic maps, 519 co-split pairs',flush=True)


def coset_worker():
    from sage.all import ZZ,matrix,vector,identity_matrix,pari
    inp=r.read(INPUT);G=matrix(ZZ,inp['gram']);assert G.det()==948 and G.is_positive_definite()
    assert int(pari(G).qfminim(3)[0])==0
    for c in inp['covers'].values():
        w=vector(ZZ,c['trace']);assert w*G*w==10
    classes={}
    for pair in inp['pairs']:
        w=[sum(inp['covers'][label]['trace'][i] for label in pair) for i in range(17)]
        key=r.pack([v%2 for v in w]);assert key
        classes.setdefault(key,{'representative':w,'pairs':[]})['pairs'].append(pair)
    rows=[];wd=WORK/'cosets';wd.mkdir(parents=True,exist_ok=True)
    for key,old in sorted(classes.items()):
        path=wd/f'{key}.json'
        if not path.exists():
            w=vector(ZZ,old['representative'])
            H=(2*identity_matrix(ZZ,17)).stack(matrix(ZZ,[w])).row_module(ZZ).basis_matrix()
            U=matrix(ZZ,pari(H*G*H.transpose()).qflllgram());assert abs(U.det())==1
            V=U.transpose()*H;raw=pari(V*G*V.transpose()).qfminim(6)
            vs=matrix(ZZ,raw[2]).columns();found=[]
            for v in vs:
                z=v*V;norm=int(z*G*z);assert norm in (4,6)
                assert all((z[i]-w[i])%2==0 for i in range(17))
                # Canonical sign for a reproducible pair of conjugate translates.
                zz=tuple(map(int,z));zz=min(zz,tuple(-x for x in zz))
                found.append({'vector':list(zz),'norm':norm})
            assert int(raw[0])==2*len(found) and len(found)<=1
            r.write_new(path,{'parity':key,'representative':list(map(int,w)),
                'reduced_basis':[[int(x) for x in row] for row in V],
                'short_representatives_up_to_sign':found,'has_norm_six':any(x['norm']==6 for x in found)})
        rows.append(r.read(path))
    bykey={x['parity']:x for x in rows};pairs=[]
    for pair in inp['pairs']:
        w=[sum(inp['covers'][label]['trace'][i] for label in pair) for i in range(17)]
        key=r.pack([v%2 for v in w]);row=bykey[key]
        d={'labels':pair,'parity':key,'has_degree_one_relation':row['has_norm_six']}
        if row['has_norm_six']:
            z=row['short_representatives_up_to_sign'][0]['vector']
            d['conjugate_generic_words']=[[(w[i]+sign*z[i])//2 for i in range(17)] for sign in (1,-1)]
        pairs.append(d)
    r.write_new(COSETS,{'schema':'rank-jump.degree-one-relation-cosets.v1','status':'PASS','cosets':rows,'pairs':pairs,
        'norm_below_four_vector_count':0,
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':'Complete norm<=6 coset enumeration on all co-split pairs. The norm-six relation exists globally but its incidence at a particular panel parameter still needs testing.'})
    print('COSETS',len(rows),'degree-one pairs',sum(x['has_degree_one_relation'] for x in pairs),flush=True)


def fibre_worker(index):
    from sage.all import QQ,EllipticCurve
    inp=r.read(INPUT);cosets=r.read(COSETS);key,old=sorted(inp['blocks'].items())[index]
    t=QQ(old['parameter'])
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
    def generic(word):return sum((int(a)*P for a,P in zip(word,basis,strict=True)),E(0))
    points={}
    for label in old['labels']:
        c=inp['covers'][label];q=ev(c['q']);assert q>0 and q.is_square();u=q.sqrt()
        x0,x1,y0,y1=[ev(c['lift'][s+'_coefficients']) for s in ('x0','x1','y0','y1')]
        pp=[E(x0+sign*u*x1,y0+sign*u*y1) for sign in (1,-1)]
        assert pp[0]!=pp[1] and pp[0]+pp[1]==generic(c['trace'])
        points[label]=pp
    rules={tuple(p['labels']):p for p in cosets['pairs']}
    edges=[];eligible=0;pair_checks=[]
    for ii,jj in combinations(range(len(old['labels'])),2):
        labels=tuple(sorted((old['labels'][ii],old['labels'][jj])));rule=rules[labels]
        if not rule['has_degree_one_relation']:continue
        eligible+=1;matches=[]
        for word_index,word in enumerate(rule['conjugate_generic_words']):
            S=generic(word)
            for a in range(2):
                for b in range(2):
                    if points[labels[0]][a]+points[labels[1]][b]==S:
                        matches.append({'word_index':word_index,'branch_indices':[a,b]})
        assert len(matches) in (0,2)
        if matches:
            assert {m['word_index'] for m in matches}=={0,1}
            m=next(m for m in matches if m['word_index']==0)
            i,j=[old['labels'].index(label) for label in labels]
            edges.append({'indices':[i,j],'labels':list(labels),'signs':[(-1)**a for a in m['branch_indices']],
                          'generic_word':rule['conjugate_generic_words'][0]})
        pair_checks.append({'labels':list(labels),'branch_matches':matches})
    from sage.all import matrix
    n=len(old['labels']);M=matrix(QQ,len(edges),n)
    for j,e in enumerate(edges):
        for i,s in zip(e['indices'],e['signs']):M[j,i]=s
    rank=int(M.rank());adj=[set() for _ in range(n)]
    for e in edges:
        i,j=e['indices'];adj[i].add(j);adj[j].add(i)
    unseen=set(range(n));components=[]
    while unseen:
        todo=[min(unseen)];seen=set(todo)
        while todo:
            i=todo.pop()
            for j in adj[i]-seen:seen.add(j);todo.append(j)
        unseen-=seen;components.append(sorted(seen))
    row={'block_key':key,'parameter':old['parameter'],'labels':old['labels'],'compatible_cover_count':n,
         'globally_degree_one_pairs_among_compatible_covers':eligible,
         'rational_degree_one_components_at_this_parameter_mod_conjugation':len(edges),
         'pair_checks':pair_checks,'signed_quotient_edges':edges,'signed_relation_rank_over_Q':rank,
         'native_quotient_rank_upper_bound_from_pair_relations':n-rank,'graph_components':components}
    r.write_new(WORK/'fibres'/str(index)/'result.json',row)
    print(key,'covers',n,'eligible',eligible,'incidences',len(edges),'quotient upper',n-rank,flush=True)


def supervise(command,folder,seconds):
    folder.mkdir(parents=True,exist_ok=True);log=folder/'worker.log';state=folder/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),*command],stdout=out,stderr=out,timeout=seconds)
                d={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:d={'status':'TIMEOUT'}
        r.write_new(state,d)
    return r.read(state),log.read_text()


def run():
    state,log=supervise(['cosets'],WORK/'coset_execution',120);print(state,log,flush=True)
    if state['status']!='COMPLETE':return
    inp=r.read(INPUT);rows=[]
    for i,(key,block) in enumerate(sorted(inp['blocks'].items())):
        if len(block['labels'])<2:
            rows.append({'block_key':key,'parameter':block['parameter'],'labels':block['labels'],
                         'compatible_cover_count':len(block['labels']),'globally_degree_one_pairs_among_compatible_covers':0,
                         'rational_degree_one_components_at_this_parameter_mod_conjugation':0,'pair_checks':[],
                         'signed_quotient_edges':[],'signed_relation_rank_over_Q':0,
                         'native_quotient_rank_upper_bound_from_pair_relations':len(block['labels']),
                         'graph_components':[[j] for j in range(len(block['labels']))]})
            continue
        folder=WORK/'fibres'/str(i);state,log=supervise(['fibre','--index',str(i)],folder,30)
        print(state,log,flush=True)
        result=folder/'result.json';rows.append(r.read(result) if result.exists() else {'block_key':key,'status':'UNKNOWN','execution':state,'log':log})
    r.write_new(OUTPUT,{'schema':'rank-jump.degree-one-relation-panel.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,COSETS,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['failure_semantics']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','cosets','fibre']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    elif a.mode=='cosets':coset_worker()
    else:fibre_worker(a.index)
