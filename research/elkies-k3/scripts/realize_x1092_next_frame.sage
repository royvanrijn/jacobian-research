#!/usr/bin/env sage-python
"""Bounded exact realization of the next prospectively ranked X1092 frame.

Reuse the certified class1 RR/section compiler with new lattice inputs.
No exceptional specialization data enter discovery or realization.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import time
import zipfile
from sage.all import *
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
BASE=ROOT/'elkies-k3/scripts/realize_x1092_class1.sage'
SOURCE=ART/'curve302_recovered_mw17_parent_v1.json'
CENSUS=ART/'det1092_pruned_rootless_j2_census_v1.json'
PRIORITY=ART/'det1092_frame_realization_priority_v1.json'
ORBITS=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
PACKETS=ART/'det1092_pruned_anchor_packets_v1.zip'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(m):return [list(map(int,r)) for r in m.rows()]
def read(p):return json.loads(p.read_text())
def emit(stage,payload):
    payload.update(schema='x1092.frame-realization.'+stage+'.v1',class_index=args.class_index,sage_version=SAGE_VERSION)
    if payload.get('family')=='x1092-j2-class1':payload['family']='x1092-j2-class'+str(args.class_index)
    p=ART/(PREFIX+'_'+stage+'_v1.json')
    if p.exists():
        raise RuntimeError('Refusing to overwrite retained proof packet: '+str(p))
    # The inherited compiler emits incremental progress repeatedly; use indexed snapshots.
    p.write_text(json.dumps(payload,indent=2)+'\n')
    print(stage,payload.get('status'),flush=True)

def compiler_emit(stage,payload):
    if stage=='sections_progress':
        stage += '_'+str(len(payload['sections']))
    emit(stage,payload)

def discover():
    G=matrix(ZZ,read(SOURCE)['generic_height_gram'])
    target=matrix(ZZ,read(CENSUS)['rootless_classes'][args.class_index-1]['gram'])
    priority=next(r for r in read(PRIORITY)['ranked_new_types'] if r['class_index']==args.class_index)
    candidates=[r for r in csv.DictReader(ORBITS.open(),delimiter='\t') if int(r['minimum_norm'])==12]
    # Coordinate complexity and mask only: no parameter/search data.
    candidates.sort(key=lambda r:(sum(abs(int(c)) for c in r['parent_MW17_w'].split()),int(r['orbit_mask'])))
    J=matrix(ZZ,[[0,1],[1,0]]);N=block_diagonal_matrix(J,-G)
    checked=[]
    for r in candidates:
        w=list(map(int,r['parent_MW17_w'].split()))
        D=vector(ZZ,[3,2]+w);O=vector(ZZ,[-1,1]+[0]*17)
        U=matrix(ZZ,[D,D+O]);W=(U*N).right_kernel_matrix();H=-W*N*W.T
        iso=pari(target).qfisom(pari(H))
        checked.append({'mask':int(r['orbit_mask']),'matches_target':bool(iso)})
        if iso:
            emit('discovery',{'status':'PASS_TARGET_FRAME_IN_GENERIC_DEGREE2_WINDOW',
                'orbit_mask':int(r['orbit_mask']),'word':w,'checked':checked,
                'prospective_priority':priority,'selection_rule':'First exact target isometry in norm12 representatives ordered by L1 coordinate complexity, then mask.',
                'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,CENSUS,PRIORITY,ORBITS,BASE)},
                'specialization_inputs':0})
            return
    emit('discovery',{'status':'UNKNOWN_NO_MATCH_IN_RETAINED_DEGREE2_WINDOW','checked':checked})

def marking():
    discovery=read(ART/(PREFIX+'_discovery_v1.json'))
    assert discovery['status']=='PASS_TARGET_FRAME_IN_GENERIC_DEGREE2_WINDOW'
    for p,digest in discovery['inputs'].items():assert sha(ROOT/p)==digest
    G=matrix(ZZ,read(SOURCE)['generic_height_gram'])
    target=matrix(ZZ,read(CENSUS)['rootless_classes'][args.class_index-1]['gram'])
    witness=discovery['prospective_priority']['witness']
    with zipfile.ZipFile(PACKETS) as z:packet=json.loads(z.read('anchor-16.json'))
    left,right=[packet['embeddings'][witness[k]] for k in ('known_embedding_index','new_embedding_index')]
    assert left['sixth_index']==right['sixth_index']==witness['shared_sixth_index']
    C=matrix(ZZ,witness['common_core_basis_in_niemeier'])
    for e in (left,right):
        B=matrix(ZZ,e['complement_basis_in_ambient']);coords=matrix(ZZ,B.T.solve_right(C.T).T)
        assert coords*matrix(ZZ,e['gram'])*coords.T==matrix(ZZ,witness['common_core_gram'])
        assert coords.row_module().saturation()==coords.row_module()
    assert matrix(ZZ,witness['common_core_gram']).det()==witness['common_core_determinant']
    WORD=discovery['word'];w=vector(ZZ,WORD);assert w*G*w==12
    J=matrix(ZZ,[[0,1],[1,0]]);N=block_diagonal_matrix(J,-G)
    D=vector(ZZ,[3,2]+WORD);O=vector(ZZ,[-1,1]+[0]*17)
    U=matrix(ZZ,[D,D+O]);W=(U*N).right_kernel_matrix();H=-W*N*W.T;T=U.stack(W)
    assert abs(T.det())==1 and T*N*T.T==block_diagonal_matrix(J,-H)
    iso=matrix(ZZ,pari(target).qfisom(pari(H)).sage())
    assert abs(iso.det())==1 and iso.T*H*iso==target
    requested=matrix(ZZ,right['gram']);reqiso=matrix(ZZ,pari(requested).qfisom(pari(H)).sage())
    assert abs(reqiso.det())==1 and reqiso.T*H*reqiso==requested
    requested_rows=reqiso.T*W
    assert requested_rows*N*requested_rows.T==-requested and U*N*requested_rows.T==0
    assert pari(H).qfminim(2)[0]==0
    assert not pari(G).qfisom(pari(H))
    class1=matrix(ZZ,read(CENSUS)['rootless_classes'][0]['gram'])
    assert not pari(class1).qfisom(pari(H))
    walls,shells=compiler['negative_walls'](D,G);assert not walls
    emit('marking',{'status':'PASS_EXACT_RATIONAL_MARKED_NEF_U_AND_FRAME_TRANSPORT',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,CENSUS,PRIORITY,PACKETS,ORBITS,BASE)},
        'requested_shared_core_witness':witness,
        'transport_scope':'Same requested J2 frame; not asserted to extend the particular shared-core identification.',
        'source_NS_gram':rows(N),'fibre_D':list(map(int,D)),'rational_zero_O':list(map(int,O)),
        'orbit_mask':discovery['orbit_mask'],'generic_trace_word':WORD,
        'transport_rows_D_D_plus_O_complement':rows(T),'transport_determinant':int(T.det()),
        'frame_gram':rows(H),'representative_columns_in_complement':rows(iso),
        'requested_embedding_columns_in_complement':rows(reqiso),
        'requested_embedding_frame_rows_in_source_NS':rows(requested_rows),
        'full_gram_identity_verified':True,'complement_integral_isometry_verified':True,
        'nonisometric_to_classes':[1,6],
        'nef_certificate':{'old_fibre_degree':2,'old_zero_intersection':1,'primitive':True,
            'effective':'RR and D.F_old=2>0 select D, since F_old is nef.',
            'old_vertical_root_rank':0,'complete_horizontal_wall_shells':shells,'negative_walls':walls,
            'completeness':'A negative-intersection irreducible curve is a (-2) component of effective D, hence has old fibre degree at most2. No vertical roots.'},
        'equation':'UNKNOWN','rational_sections':'UNKNOWN','parameter_search':'BLOCKED'})
    WORK.mkdir(parents=True,exist_ok=True);save((G,N,D,O,H,T),str(WORK/'marking.sobj'))

def section_plan():
    G,N,D,O,H,T=load(str(WORK/'marking.sobj'));inverse=T.inverse()
    discovery=read(ART/(PREFIX+'_discovery_v1.json'))
    short=matrix(ZZ,pari(G).qfminim(8)[2]);candidates=[]
    for col in short.columns()+[vector(ZZ,discovery['word'])]:
        for word in (col,-col):
            div=vector(ZZ,[(word*G*word-2)//2,1]+list(word))
            if div*N*D!=1:continue
            child=div*inverse;assert child[1]==1 and all(c in ZZ for c in child)
            candidates.append((word,vector(ZZ,child[2:])))
    candidates.sort(key=lambda x:(x[0]*G*x[0],sum(abs(c) for c in x[0]),tuple(x[0])))
    helper=runpy.run_path(str(ROOT/'elkies-k3/scripts/plan_r17_norm12_direct_section_basis.sage'))
    # Full rational rank is insufficient for the modulo-two production gate.
    # improve_basis, unlike reduce_basis, rejects a remaining nontrivial index.
    result=helper['improve_basis']([],candidates);glue_word=None;glue=[]
    if result is None:
        for row in csv.DictReader(ORBITS.open(),delimiter='\t'):
            if int(row['minimum_norm'])!=10:continue
            word=vector(ZZ,list(map(int,row['parent_MW17_w'].split())))
            for sign in (1,-1):
                v=sign*word;div=vector(ZZ,[2,2]+list(v))
                if div*N*D==1:
                    child=div*inverse
                    glue.append((sum(abs(c) for c in v),int(row['orbit_mask']),sign,v,vector(ZZ,child[2:])))
        glue.sort(key=lambda x:tuple(x[:3]))
        for _,mask,sign,v,child in glue:
            result=helper['improve_basis']([child],candidates)
            if result is not None:
                glue_word=list(map(int,v));break
    if result is None:
        emit('section_plan',{'status':'UNKNOWN_BOUNDED_SATURATED_SECTION_SPAN',
            'candidate_count':len(candidates),'bisection_candidate_count':len(glue)})
        return
    selected,B=result;assert abs(B.det())==1
    emit('section_plan',{'status':'PASS_EXACT_RANK17_SECTION_CLASS_PLAN','height_bound':12,
        'candidate_count':len(candidates),'selected_source_words':rows(matrix(ZZ,selected)),
        'glue_bisection_word':glue_word,'child_frame_coordinates':rows(B),'subgroup_index':1,
        'height_gram':rows(B*H*B.T),'height_determinant':int((B*H*B.T).det()),
        'rational_coordinate_recovery':'PENDING',
        'saturation_gate':'A rank17 but index>1 old-section subgroup triggers the retained norm10 bisection completion.'})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['discover','marking','trace','equation','section_plan','sections','normalize','compact'])
    p.add_argument('--class-index',type=int,default=3)
    args=p.parse_args();assert args.class_index not in (1,6)
    PREFIX='x1092_class'+str(args.class_index)+'_realization'
    WORK=ROOT/('artifacts/local/elkies-k3/x1092-class'+str(args.class_index)+'-realization-v1')
    WORK.mkdir(parents=True,exist_ok=True)
    compiler=runpy.run_path(str(BASE))
    if args.stage in ('discover','marking','section_plan'):globals()[args.stage]()
    else:
        discovery=read(ART/(PREFIX+'_discovery_v1.json'))
        assert discovery['status']=='PASS_TARGET_FRAME_IN_GENERIC_DEGREE2_WINDOW'
        env=compiler[args.stage].__globals__
        env.update(WORK=WORK,PREFIX=PREFIX,WORD=discovery['word'],MASK=discovery['orbit_mask'],emit=compiler_emit)
        compiler[args.stage]()
