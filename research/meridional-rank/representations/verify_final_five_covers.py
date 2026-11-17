#!/usr/bin/env python3
"""Independent finite replay for five source models; no producer imports.

The written proof separately realizes the manifolds, verifies canonical JSJ
and applies the imported classical rank2 necessary condition.
"""
from copy import deepcopy
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
from verify_compact_covers import replay,clean,reduce_word,abelianization,replace
from verify_port_assemblies import check_fragment,check_join
from verify_mutation_covers import regions_and_graph
from verify_three_piece_covers import all_cuts,canonical_ports,root_identity,check_root as check_old_root
from verify_four_piece_covers import H_normal,stack,inverse,subst,power,cyclic_key,commutator,det
from verify_double_transpositions import geometric_crossings,check_trace

HERE=Path(__file__).resolve().parent
NAMES=['13a_234','13a_242','13a_243','13a_244','13a_265']
SPECS={
    '13a_234':{'P_type':233,'m':2,'d':[2,-1],'Q_type':'root_right','conjugator':[-1],'k_h_sign':1,'outer_order':'ks','M':[[2,3],[-1,-2]],'D':389},
    '13a_242':{'P_type':239,'m':3,'d':[2,-1],'Q_type':'root_right','conjugator':[3],'k_h_sign':1,'outer_order':'sk','M':[[-1,-2],[0,1]],'D':359},
    '13a_243':{'P_type':239,'m':2,'d':[1,-2],'Q_type':'root_left','conjugator':[1],'k_h_sign':1,'outer_order':'sk','M':[[-1,-3],[1,2]],'D':337},
    '13a_244':{'P_type':239,'m':2,'d':[2,-1],'Q_type':'root_right','conjugator':[-1],'k_h_sign':1,'outer_order':'ks','M':[[-1,-3],[1,2]],'D':311},
    '13a_265':{'P_type':253,'m':2,'d':[2,-1,-1],'Q_type':'root_squared','conjugator':[],'k_h_sign':-1,'outer_order':'sk','M':[[-1,-3],[1,2]],'D':311}}


def check_source_rows(result,plan):
    assert result['status']=='COMPLETE_FROZEN_SELECTION'
    assert plan['targets']==NAMES and plan['controls']==['13a_233','13a_239','13a_253']
    assert [r['name'] for r in result['results']]==plan['controls']+plan['targets']
    panel={r['name']:r['pd'] for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    old={r['name']:r for r in json.loads((HERE/'remaining_cover_structure_results.json').read_text())['results']}
    checked=[]
    for row in result['results']:
        pd=row['pd'];assert pd==panel[row['name']]
        graph=regions_and_graph(pd,row['checkerboards']['selected'])
        raw=[reduce_word([x for x in (o+1,-a-1,o+1,-b-1) if abs(x)!=1]) for o,a,b in geometric_crossings(pd)]
        assert raw==row['full_raw_core_relators']
        assert row['compact_core']['initial']=={'generators':list(range(2,14)),'relators':clean(raw)}
        generators,relators=replay(row['compact_core']);assert len(generators)==len(relators)==3
        smith=abelianization(generators,relators);assert smith==[1,1,graph['spanning_trees']]
        cuts=all_cuts(pd);assert cuts==row['candidate_cuts'] and row['cut_status']=='ALL_CANDIDATES'
        assert len(cuts)==len(row['cuts'])
        for expected,cut in zip(cuts,row['cuts']):
            assert all(cut[k]==v for k,v in expected.items())
            for field,fragment in zip(('inside_crossings','outside_crossings'),cut['fragments']):
                check_fragment(fragment,[pd[i] for i in cut[field]])
            check_join(*cut['fragments'],pd)
        if row['name'] in plan['controls']:
            previous=old[row['name']]
            for field in ('pd','full_raw_core_relators','candidate_cuts','cuts'):
                if field=='cuts':
                    # Group elapsed times are observational, not mathematical inputs.
                    assert len(row[field])==len(previous[field])
                    for A,B in zip(row[field],previous[field]):
                        assert A['inside_crossings']==B['inside_crossings'] and A['boundary']==B['boundary']
                        for a,b in zip(A['fragments'],B['fragments']):
                            assert a['pd']==b['pd'] and a['arc_edge_classes']==b['arc_edge_classes']
                            assert a['formal_core_reduction']['final']==b['formal_core_reduction']['final']
                else:assert row[field]==previous[field]
            assert row['compact_core']['final']==previous['compact_core']['final']
        checked.append({'name':row['name'],'H1_smith':smith,'source_generators':3,'balanced_cuts':len(cuts),'fragment_reductions':2*len(cuts)})
    assert sum(r['fragment_reductions'] for r in checked)==32
    return checked


def filling_certificate(m,r,q,order):
    assert (m,abs(r)) in ((2,3),(3,2))
    identity=(0,1,2);trans=(1,0,2);cycle=(1,2,0)
    s,c=(trans,cycle) if m==2 else (cycle,trans)
    compose=lambda a,b:tuple(a[b[i]] for i in range(3))
    def exponent(p,n):
        out=identity
        for _ in range(n%6):out=compose(out,p)
        return out
    k=compose(c,exponent(s,-1)) if order=='ks' else compose(exponent(s,-1),c)
    f=exponent(s,m)
    actual_c=compose(k,s) if order=='ks' else compose(s,k)
    assert f==identity and actual_c==c
    assert compose(k,f)==compose(f,k)
    assert compose(exponent(actual_c,r),exponent(f,q))==identity
    assert compose(s,c)!=compose(c,s)
    seen={identity};queue=[identity]
    for a in queue:
        for b in (s,k):
            value=compose(a,b)
            if value not in seen:seen.add(value);queue.append(value)
    assert len(seen)==6
    return {'s':s,'k':k,'c':c,'fiber':f,'image_order':6,'noncommuting_words':['s','c'],'filling_relation_passes':True}


def check_model(row,source):
    name=row['name'];spec=SPECS[name];m=spec['m'];d=spec['d'];M=spec['M'];kind=spec['Q_type']
    for k in ('P_type','m','d','Q_type','conjugator','k_h_sign','outer_order'):assert row[k]==spec[k]
    assert row['P_to_C']==M and row['determinant']==spec['D']
    cut=source['cuts'][0];pieces,ports=canonical_ports(cut)
    assert row['boundary']==cut['boundary'] and row['ports']==ports
    assert row['pieces']==[{**p,'arc_words':{str(g):w for g,w in p['arc_words'].items()}} for p in pieces]
    assert pieces[0]['relator']==[-2,-1,-1,-2,1]
    tP=[1,1,2];assert cyclic_key(pieces[0]['relator'])==cyclic_key(tP+tP+[-1]*3)
    Pmap={1:[1],2:[-1,-1,2]};Qmap={1:[3],2:[1]} if kind=='root_left' else {1:[1],2:[3]}
    assert row['P_images']=={str(g):w for g,w in Pmap.items()} and row['Q_images']=={str(g):w for g,w in Qmap.items()}
    tQ=([2]*m+[1]) if kind=='root_right' else ([-1,-1,2] if kind=='root_left' else [2,2,1,1])
    base=2 if kind=='root_left' else 1
    assert cyclic_key(pieces[1]['relator'])==cyclic_key(tQ+tQ+[-base]*3)
    assert root_identity(subst(tQ,Qmap)+[-2],d,m)
    assert root_identity(subst(pieces[1]['relator'],Qmap),d,m)
    assert H_normal(subst(tP,Pmap))==H_normal([2])
    h=[1]*3;e=stack(Pmap[2]+inverse(h));k=power(h,spec['k_h_sign']);f=[3]*m
    c=(k+[3]) if spec['outer_order']=='ks' else ([3]+k)
    assert row['outer_c_f']==[c,f]
    Pcoords={233:[[0,1],[1,1],[-1,2],[2,0]],239:[[0,1],[1,0],[-1,1],[0,0]],253:[[-2,1],[1,0],[-1,1],[0,0]]}[spec['P_type']]
    Qcoords={'13a_234':[[0,0],[2,-1],[1,-1],[1,0]],'13a_242':[[0,0],[1,-1],[1,0],[2,-1]],
        '13a_243':[[0,0],[2,-1],[1,-1],[3,-2]],'13a_244':[[0,0],[2,-1],[1,-1],[3,-2]],'13a_265':[[0,0],[0,1],[-1,1],[1,0]]}[name]
    assert row['P_port_coordinates_e_h']==Pcoords and row['Q_port_coordinates_c_f']==Qcoords
    expectedP=[stack(power(e,A)+power(h,B)) for A,B in Pcoords]
    expectedQ=[stack(power(c,A)+power(f,B)) for A,B in Qcoords]
    assert row['P_ports_expected']==expectedP and row['Q_ports_expected']==expectedQ
    b=spec['conjugator']
    for port,P,Q in zip(ports,expectedP,expectedQ):
        assert H_normal(subst(port['P'],Pmap)+inverse(P))==(0,()),(name,port['edge'],'P')
        actual=inverse(b)+subst(port['Q'],Qmap)+b
        assert root_identity(actual+inverse(Q),d,m),(name,port['edge'],'Q')
    assert root_identity(commutator(c,f),d,m)
    assert not root_identity(commutator(h,[3]),d,m)
    N=[[0,1],[spec['k_h_sign'],-1 if kind=='root_squared' else 0]]
    assert row['C_to_Q']==N and abs(det(N))==abs(det(M))==1
    assert row['fiber_intersections']==[abs(M[0][1]),1]
    assert row['proposed_JSJ_bases']==['D(2,3)',f'A({m})','D(2,3)']
    base_port=0 if spec['P_type']==233 else 3
    z=[-M[0][1],-M[1][1]] if base_port==0 else Qcoords[3]
    assert row['z_in_c_f']==z
    assert all([sum(M[j][k]*p[k] for k in (0,1))+z[j] for j in (0,1)]==q for p,q in zip(Pcoords,Qcoords))
    diffs=[[p[j]-Pcoords[base_port][j] for j in (0,1)] for i,p in enumerate(Pcoords) if i!=base_port]
    assert any(abs(det(list(pair)))==1 for pair in combinations(diffs,2))
    assert abelianization(*replay(source['compact_core']))==[1,1,spec['D']]
    r,q=M[0][1],M[1][1];certificate=filling_certificate(m,r,q,spec['outer_order'])
    return {'name':name,'root_degree':m,'source_H1_order':spec['D'],'P_to_C':M,'C_to_Q':N,
        'central_outer_order':spec['outer_order'],'inner_fiber_sign':spec['k_h_sign'],
        'fiber_intersections':row['fiber_intersections'],'central_filling_slope_c_f':[r,q],
        'nonabelian_filling_certificate':certificate,'proposed_JSJ_pieces':3,'finite_inputs':'PASS'}


def main():
    started=time.monotonic();plan=json.loads((HERE/'final_five_peripheral_plan.json').read_text())
    structural=json.loads((HERE/'final_five_structure_results.json').read_text())
    splan=json.loads((HERE/'final_five_structure_plan.json').read_text())
    result=json.loads((HERE/'final_five_peripheral_results.json').read_text())
    assert result['status']=='COMPLETE_FROZEN_PERIPHERAL_INPUTS' and plan['targets']==NAMES
    for record in (structural,result):
        for p,digest in record['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest,p
    checked=check_source_rows(structural,splan);sources={r['name']:r for r in structural['results']}
    old=json.loads((HERE/'three_piece_cover_results.json').read_text())
    controls=[check_old_root(row,sources[row['name']]) for row in old['roots']]
    models=[check_model(row,sources[row['name']]) for row in result['targets']];assert [r['name'] for r in models]==NAMES
    previous=json.loads((HERE/'three_piece_cover_verification.json').read_text())['roots']
    assert json.loads(json.dumps(controls))==previous
    uppers={r['name']:r for r in json.loads((HERE/'joint_a5_results.json').read_text())['results']};upper=[]
    for name in NAMES:
        row=uppers[name];assert len(row['seeds'])==4
        check_trace(geometric_crossings(sources[name]['pd']),row['seeds'],row['seed_trace'])
        upper.append({'name':name,'seeds':row['seeds'],'propagation_steps':len(row['seed_trace']),'upper4_input_checked':True})
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    bad=deepcopy(result['targets'][0]);bad['Q_ports_expected'][1][0]*=-1
    reject('changed source port',lambda:check_model(bad,sources[bad['name']]))
    bad=deepcopy(result['targets'][0]);bad['outer_c_f'][0]=[3,1,1,1]
    reject('outer ks replaced by sk without its conjugation',lambda:check_model(bad,sources[bad['name']]))
    bad=deepcopy(result['targets'][1]);bad['conjugator']=[]
    reject('lost Q242 base-path conjugation',lambda:check_model(bad,sources[bad['name']]))
    bad=deepcopy(result['targets'][-1]);bad['C_to_Q'][1][0]=1
    reject('wrong Q265 inner fiber sign',lambda:check_model(bad,sources[bad['name']]))
    bad=deepcopy(structural);bad['results'][3]['cuts'][0]['fragments'][0]['arc_edge_classes'][0].append(9)
    reject('wrong physical source arc',lambda:check_source_rows(bad,splan))
    reject('genuine meridian slope r1 has no claimed nonabelian filling obstruction',lambda:filling_certificate(2,1,0,'sk'))
    relation=commutator([1],[2,2]);assert not subst(relation,{1:[],2:[2]}) and not subst(relation,{1:[-2],2:[2]})
    assert time.monotonic()-started<plan['limits']['maximum_seconds']
    files=['verify_final_five_covers.py','final_five_structure_plan.json','final_five_structure_results.json',
        'final_five_peripheral_plan.json','final_five_peripheral_results.json','parabolic_panel.json',
        'verify_compact_covers.py','verify_port_assemblies.py','verify_mutation_covers.py','verify_three_piece_covers.py',
        'verify_four_piece_covers.py','verify_double_transpositions.py','parabolic_lifts.py',
        'three_piece_cover_results.json','three_piece_cover_verification.json','three_piece_source_audit.json',
        'joint_a5_results.json','joint_a5_verification.json']
    receipt={'status':'PASS independent final-five source, peripheral, filling and upper4 inputs',
        'source_audit':checked,'positive_controls':controls,'targets':models,'upper4':upper,'corruptions_rejected':rejected,
        'cyclic_meridian_filling_control':'k=1 or sk=1 leaves <s> with no relation',
        'actual_rank_transfers':0,
        'assurance':'All eight literal sources checked directly; all5new lower inputs and individual uppers pass. No mutation or cross-knot rank transfer. Finite normal forms, primitive slopes and nonabelian filling quotients are independently replayed without producer imports; legacy parser retained only for index transport. Manifold realization, canonical JSJ and imported BW05 case6 necessity remain the separate written uniform proof, unformalized and externally unreviewed.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'final_five_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'new_targets':len(models),'positive_controls':len(controls),
        'source_groups':len(checked),'fragment_reductions':32,'upper4':len(upper),'corruptions':len(rejected)}))


if __name__=='__main__':main()
