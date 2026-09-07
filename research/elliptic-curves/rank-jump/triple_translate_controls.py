#!/usr/bin/env python3
"""Three lattice-selected translate controls for one fixed native triple."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'TRIPLE_TRANSLATE_CONTROL_PROTOCOL.json'
BASE=r.OUT/'rank_jump_native_triple_intersection_inputs_v1.json'
OLD=r.OUT/'rank_jump_native_triple_intersection_v1.json'
SELECTION=r.OUT/'rank_jump_triple_translate_selection_v1.json'
OUTPUT=r.OUT/'rank_jump_triple_translate_controls_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-triple-translate-controls-v1'


def casepath(i):return r.OUT/f'rank_jump_triple_translate_control_{i}_inputs_v1.json'
def resultpath(i):return r.OUT/f'rank_jump_triple_translate_control_{i}_v1.json'


def select(gram,traces):
    """Intentionally accepts no oracle relation, parameter, point or rank."""
    from sage.all import ZZ,matrix,vector,identity_matrix,pari
    G=matrix(ZZ,gram);tau=[vector(ZZ,x) for x in traces];w=tau[0]-tau[1]+tau[2]
    H=(2*identity_matrix(ZZ,17)).stack(matrix(ZZ,[w])).row_module(ZZ).basis_matrix()
    assert H.nrows()==17 and abs(H.det())==2**16
    U=matrix(ZZ,pari(H*G*H.transpose()).qflllgram());V=U.transpose()*H;assert abs(V.det())==2**16
    raw=pari(V*G*V.transpose()).qfminim(10);vecs=matrix(ZZ,raw[2]).columns();signed=[]
    for v in vecs:
        rr=v*V;norm=int(rr*G*rr);assert 0<norm<=10 and all((rr[i]-w[i])%2==0 for i in range(17))
        signed.extend((rr,-rr))
    assert len(signed)==int(raw[0])
    def key(v):return sum(abs(int(x)) for x in v),max(abs(int(x)) for x in v),tuple(map(int,v))
    eligible={tuple(map(int,(rr+w)/2)) for rr in signed if rr*G*rr==10 and ((rr+w)/2)*G*((rr+w)/2)==10}
    canonical={min([s]+([tuple(map(int,w-vector(ZZ,s)))] if tuple(w-vector(ZZ,s)) in eligible else []),key=key) for s in eligible}
    chosen=sorted(canonical,key=key)[:3]
    return {'trace_sum':list(map(int,w)),'reduced_coset_basis':[[int(x) for x in row] for row in V],
            'signed_vectors':[{'vector':list(map(int,v)),'norm':int(v*G*v)} for v in sorted(signed,key=lambda v:tuple(v))],
            'signed_short_vector_count':len(signed),'eligible_words':list(map(list,sorted(eligible,key=key))),
            'conjugacy_class_count':len(canonical),'selected_words':list(map(list,chosen))}


def capture():
    base=r.read(BASE);geometry={k:base[k] for k in ('A','B','sections','gram','covers')}
    selection=select(geometry['gram'],[c['published_basis_w'] for c in geometry['covers']])
    # Oracle comparison and cohort values enter only after deterministic selection.
    from sage.all import ZZ,vector
    w=vector(ZZ,selection['trace_sum']);oracle=vector(ZZ,[ZZ(x) for x in base['generic_word']])
    selection['matches_oracle_conjugacy_class']=[vector(ZZ,x)==oracle or vector(ZZ,x)==w-oracle for x in selection['selected_words']]
    r.write_new(SELECTION,{'schema':'rank-jump.triple-translate-selection.v1',**selection,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,BASE,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']})
    for i,word in enumerate(selection['selected_words']):
        r.write_new(casepath(i),{'schema':'rank-jump.triple-translate-control-input.v1',**geometry,
            'generic_word':list(map(str,word)),'parameters':base['parameters'],'selection_index':i,
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,BASE,SELECTION)}})
    print('Selected',len(selection['selected_words']),'from',selection['conjugacy_class_count'],'classes; oracle matches',selection['matches_oracle_conjugacy_class'])


def worker(i,check=False):
    import native_triple_intersection as t
    t.INPUT=casepath(i);t.PROTOCOL=PROTOCOL
    data=t.compute()
    if check:assert r.read(resultpath(i))==data
    else:r.write_new(resultpath(i),data)
    print('Control',i,'degree',data['intersection_degree'],'factor degrees',[len(x['coefficients'])-1 for x in data['factorization']],flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i,word in enumerate(r.read(SELECTION)['selected_words']):
        log=WORK/f'{i}.log';execution=WORK/f'{i}_execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i)],stdout=out,stderr=out,timeout=60)
                    status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
            r.write_new(execution,status)
        row={'index':i,'execution':r.read(execution),'log':log.read_text()}
        if resultpath(i).exists():
            data=r.read(resultpath(i));row.update({'result':str(resultpath(i).relative_to(r.ROOT)),
                'intersection_degree':data['intersection_degree'],'factor_degrees':[len(f['coefficients'])-1 for f in data['factorization']],
                'frozen_hits':[x['source_id'] for x in data['parameter_checks'] if x['on_relation_locus']]})
        rows.append(row);print(row,flush=True)
    paths=[PROTOCOL,BASE,OLD,SELECTION,Path(__file__),HERE/'native_triple_intersection.py',HERE/'retrospective.py']
    for i in range(len(rows)):
        paths.append(casepath(i))
        if resultpath(i).exists():paths.append(resultpath(i))
    r.write_new(OUTPUT,{'schema':'rank-jump.triple-translate-controls.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker','check']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:worker(a.case,check=a.mode=='check')
