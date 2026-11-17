#!/usr/bin/env python3
"""Independent exact source, normal-form, slope, mutation and upper replay.

No presentation simplifier, graph producer or search routine is imported.
Uniform topology and the classical rank theorem are companion-proof inputs.
"""
from collections import defaultdict
from copy import deepcopy
import hashlib
from itertools import combinations
import json
from pathlib import Path
from verify_compact_covers import replay as presentation_replay,clean,reduce_word,abelianization,replace
from verify_port_assemblies import check_fragment,check_join
from verify_mutation_covers import regions_and_graph,replay as graph_replay
from verify_four_piece_covers import H_normal,stack,inverse,subst,power,cyclic_key,commutator,det
from verify_double_transpositions import geometric_crossings,check_trace

HERE=Path(__file__).resolve().parent


def all_cuts(pd):
    # Enumerate masks rather than the producer's combinations of crossing lists.
    edges=defaultdict(set)
    for c,row in enumerate(pd):
        for e in row:edges[e].add(c)
    def connected(S):
        reached={min(S)}
        while True:
            enlarged=reached|{v for ends in edges.values() if ends&reached for v in ends&S}
            if enlarged==reached:return reached==S
            reached=enlarged
    cuts=[]
    for mask in range(1,1<<13,2):
        A={i for i in range(13) if mask>>i&1};B=set(range(13))-A
        if len(A) not in (6,7):continue
        boundary=sorted(e for e,ends in edges.items() if ends&A and ends&B)
        if len(boundary)==4 and connected(A) and connected(B):
            cuts.append({'inside_crossings':sorted(A),'outside_crossings':sorted(B),'boundary':boundary})
    return sorted(cuts,key=lambda r:(r['inside_crossings']!=list(range(7)),r['inside_crossings']))


def root_identity(word,d,m):
    """Complete H*_(d=s^m)<s> normal form, H=<a,t|a^3=t^2>."""
    exponent=lambda w:sum((2 if abs(x)==1 else 3)*(1 if x>0 else -1) for x in w)
    slope=exponent(d);assert abs(slope)==1 and m in (2,3)
    def member(w):
        k=exponent(w)//slope
        assert abs(k)*len(d)<10000
        return k if H_normal(w)==H_normal(power(d,k)) else None
    pieces=[]
    for x in word:
        kind='S' if abs(x)==3 else 'H';v=(1 if x>0 else -1) if kind=='S' else [x]
        if pieces and pieces[-1][0]==kind:pieces[-1][1]+=v
        else:pieces.append([kind,v])
    for _ in range(10000):
        merged=[]
        for kind,v in pieces:
            if (kind=='S' and v==0) or (kind=='H' and H_normal(v)==(0,())):continue
            if merged and merged[-1][0]==kind:merged[-1][1]+=v
            else:merged.append([kind,v])
        pieces=merged
        for i,(kind,v) in enumerate(pieces):
            if kind=='S' and v%m==0:
                pieces[i]=['H',power(d,v//m)];break
        else:
            for i in range(1,len(pieces)-1):
                if pieces[i][0]=='H':
                    k=member(pieces[i][1])
                    if k is not None:pieces[i]=['S',m*k];break
            else:return not pieces or (len(pieces)==1 and pieces[0][0]=='H' and H_normal(pieces[0][1])==(0,()))
            continue
        continue
    raise AssertionError('normal form cap')


def check_small_audit(report,plan):
    assert report['status']=='COMPLETE_FROZEN_SELECTION'
    for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest,p
    panel={r['name']:r['pd'] for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    assert [r['name'] for r in report['results']]==plan['controls']+plan['representatives']
    checked=[]
    for row in report['results']:
        pd=panel[row['name']];assert pd==row['pd'];regions_and_graph(pd,row['checkerboards']['selected'])
        raw=[reduce_word([x for x in (o+1,-a-1,o+1,-b-1) if abs(x)!=1]) for o,a,b in geometric_crossings(pd)]
        assert raw==row['full_raw_core_relators']
        assert row['compact_core']['initial']=={'generators':list(range(2,14)),'relators':clean(raw)}
        gens,rels=presentation_replay(row['compact_core']);assert len(gens)==len(rels)==3
        smith=abelianization(gens,rels)
        assert row['candidate_cuts']==all_cuts(pd) and row['cut_status']=='ALL_CANDIDATES'
        assert len(row['cuts'])==len(row['candidate_cuts'])
        for cut,expected in zip(row['cuts'],row['candidate_cuts']):
            assert all(cut[k]==v for k,v in expected.items())
            for field,fragment in zip(('inside_crossings','outside_crossings'),cut['fragments']):
                check_fragment(fragment,[pd[i] for i in cut[field]])
            check_join(*cut['fragments'],pd)
        checked.append({'name':row['name'],'core_generators':3,'H1_smith':smith,'cuts':len(row['cuts']),
            'fragment_reductions':2*len(row['cuts'])})
    old={r['name']:r['final'] for r in json.loads((HERE/'compact_cover_results.json').read_text())['targets']}
    for row in report['results'][:2]:assert row['compact_core']['final']==old[row['name']]
    return checked


def canonical_ports(cut):
    pieces=[];ports=[]
    for fragment in cut['fragments']:
        reduction=fragment['formal_core_reduction'];rename={g:[i+1] for i,g in enumerate(reduction['final']['generators'])}
        words={int(g):replace(w,rename) for g,w in reduction['original_generators_in_final'].items()};words[1]=[]
        pieces.append({'relator':replace(reduction['final']['relators'][0],rename),'arc_words':words})
    for edge in cut['boundary']:
        words=[]
        for frag,piece in zip(cut['fragments'],pieces):
            j=next(i+1 for i,g in enumerate(frag['arc_edge_classes']) if edge in g);words.append(piece['arc_words'][j])
        ports.append({'edge':edge,'P':words[0],'Q':words[1]})
    return pieces,ports


def check_root(row,source):
    name=row['name'];i=['13a_233','13a_239','13a_253'].index(name)
    cut=source['cuts'][0];pieces,ports=canonical_ports(cut)
    expected_pieces=[{**p,'arc_words':{str(g):w for g,w in p['arc_words'].items()}} for p in pieces]
    assert row['pieces']==expected_pieces and row['port_words']==ports and row['boundary']==cut['boundary']
    assert pieces[0]['relator']==[-2,-1,-1,-2,1]
    tP=[1,1,2];assert cyclic_key(pieces[0]['relator'])==cyclic_key(tP+tP+[-1]*3)
    Qtarget=[[-1,-1,2],[1,1,1,2],[2,2,1,1]][i];qbase=2 if i<2 else 1
    assert cyclic_key(pieces[1]['relator'])==cyclic_key(Qtarget+Qtarget+[-qbase]*3)
    m=[2,3,2][i];d=[[1,-2],[2,-1],[2,-1,-1]][i]
    Pmap={1:[1],2:[-1,-1,2]};Qmap=({1:[3],2:[1]} if i<2 else {1:[1],2:[3]})
    conjugator=[1] if i==0 else []
    assert row['P_images']=={str(g):w for g,w in Pmap.items()} and row['Q_images']=={str(g):w for g,w in Qmap.items()}
    assert row['root_degree']==m and row['root_word_d']==d and row['Q_conjugator']==conjugator
    # The definitions of t are inverse eliminations, not only homomorphisms.
    assert root_identity(subst(Qtarget,Qmap)+[-2],d,m)
    assert root_identity(subst(pieces[1]['relator'],Qmap),d,m)
    assert H_normal(subst(tP,Pmap))==H_normal([2])
    h=[1]*3;e=stack(Pmap[2]+inverse(h));f=[3]*m;eta=[1,-1,1][i];a=power([3],eta)+h
    assert row['central_outer_eta']==eta and row['central_outer_a_f']==[a,f]
    expectedP=[[h,stack(e+h),stack(h+h+inverse(e)),e+e],
        [h,e,stack(h+inverse(e)),[]],[stack(h+power(e,-2)),e,stack(h+inverse(e)),[]]][i]
    expectedQ=[[[],stack(inverse(f)+a+a),stack(inverse(f)+a),a],
        [[],stack(f+a),a,stack(f+a+a)],[[],f,stack(inverse(a)+f),a]][i]
    assert row['P_ports_expected']==expectedP and row['Q_ports_expected']==expectedQ
    for port,P,Q in zip(ports,expectedP,expectedQ):
        assert H_normal(subst(port['P'],Pmap)+inverse(P))==(0,())
        actual=inverse(conjugator)+subst(port['Q'],Qmap)+conjugator
        assert root_identity(actual+inverse(Q),d,m),(name,port['edge'])
    assert root_identity(commutator(a,f),d,m)
    assert not root_identity(commutator(h,[3]),d,m)
    matrices=[[[2,3],[-1,-2]],[[-1,-2],[0,-1]],[[-1,-3],[1,2]]]
    M=matrices[i];N=[[0,1],[1,-1 if i==2 else 0]]
    assert row['P_to_C_matrix_e_h_in_a_f']==M and row['C_to_Q_matrix_k_f_in_e_h']==N
    assert abs(det(M))==abs(det(N))==1
    assert row['fiber_intersections']==[abs(M[0][1]),1]
    pcoords=[[(0,1),(1,1),(-1,2),(2,0)],[(0,1),(1,0),(-1,1),(0,0)],[(-2,1),(1,0),(-1,1),(0,0)]][i]
    qcoords=[[(0,0),(2,-1),(1,-1),(1,0)],[(0,0),(1,1),(1,0),(2,1)],[(0,0),(0,1),(-1,1),(1,0)]][i]
    z=[(-3,2),(2,1),(1,0)][i]
    assert all(tuple(sum(M[j][k]*p[k] for k in (0,1))+z[j] for j in (0,1))==q for p,q in zip(pcoords,qcoords))
    # After eliminating z at its unit port, some two differences form a
    # unimodular basis in(e,h); hence all and only the asserted edge equalities.
    base=[0,3,3][i];differences=[tuple(p[j]-pcoords[base][j] for j in (0,1)) for k,p in enumerate(pcoords) if k!=base]
    assert any(abs(det([list(v) for v in pair]))==1 for pair in combinations(differences,2))
    D=[403,361,337][i];assert row['expected_H1_order']==D
    assert abelianization(*presentation_replay(source['compact_core']))==[1,1,D]
    r=M[0][1];s=M[1][1]
    # Explicit nonabelian S3 image of the central filling: f=1,sroot^m=1,
    # a^r=1.  k=sroot^-eta*a ensures a=sroot^eta*k.
    trans=(1,0,2);cycle=(1,2,0);root,outer=(trans,cycle) if m==2 else (cycle,trans)
    compose=lambda u,v:tuple(u[v[j]] for j in range(3))
    def perm_power(p,n):
        result=(0,1,2)
        for _ in range(n%6):result=compose(result,p)
        return result
    k=compose(perm_power(root,-eta),outer);fiber=perm_power(root,m)
    assert fiber==(0,1,2) and compose(perm_power(root,eta),k)==outer
    assert perm_power(outer,r)==(0,1,2) and compose(root,outer)!=compose(outer,root)
    return {'name':name,'root_degree':m,'P_to_C_matrix':M,'C_to_Q_matrix':N,'fiber_intersections':row['fiber_intersections'],
        'H1_order':D,'central_filling_slope_a_f':[r,s],
        'central_filling_nonabelian_S3':{'root':root,'k':k,'a':outer},'JSJ_pieces':3,'finite_inputs':'PASS'}


def main():
    result=json.loads((HERE/'three_piece_cover_results.json').read_text());plan=json.loads((HERE/'three_piece_cover_inputs.json').read_text())
    structure=json.loads((HERE/'remaining_cover_structure_results.json').read_text())
    smallplan=json.loads((HERE/'remaining_cover_structure_inputs.json').read_text())
    assert result['status']=='PARTIAL_TRANSFER_INPUTS'
    for p,digest in result['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest,p
    small=check_small_audit(structure,smallplan);by_name={r['name']:r for r in structure['results']}
    roots=[check_root(row,by_name[row['name']]) for row in result['roots']];assert [r['name'] for r in roots]==plan['roots']
    panel={r['name']:r['pd'] for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    diagrams={}
    for D,names in plan['blocks'].items():
        for name in names:
            data=result['diagrams'][name];assert data['pd']==panel[name]
            diagrams[name]=regions_and_graph(data['pd'],data['selected'])
        assert diagrams[names[0]]['spanning_trees']==int(D)
    expected=[(names[0],name) for names in plan['blocks'].values() for name in names[1:]]
    assert [(p['source'],p['target']) for p in result['pairs']]==expected
    pairs=[];unproved=[]
    for pair in result['pairs']:
        A=result['diagrams'][pair['source']]['selected'];B=result['diagrams'][pair['target']]['selected']
        if pair['certificate']['status'].startswith('CERTIFIED'):
            assert A['checkerboard_parity']==B['checkerboard_parity']
            pairs.append({'source':pair['source'],'target':pair['target'],**graph_replay(A,B,pair['certificate'])})
        else:
            assert pair['certificate']['status']=='NO_CERTIFICATE'
            unproved.append({'source':pair['source'],'target':pair['target'],
                'source_determinant':diagrams[pair['source']]['spanning_trees'],
                'target_determinant':diagrams[pair['target']]['spanning_trees'],
                'status':'NO_TRANSFER_CLAIM'})
    assert not pairs and len(unproved)==5
    sources={r['name']:r for r in json.loads((HERE/'joint_a5_results.json').read_text())['results']};upper=[]
    for name in diagrams:
        row=sources[name];assert len(row['seeds'])==4
        check_trace(geometric_crossings(panel[name]),row['seeds'],row['seed_trace'])
        upper.append({'name':name,'seeds':row['seeds'],'steps':len(row['seed_trace']),'upper4_input_checked':True})
    assert len(upper)==8
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    for label,field in [('changed port word','Q_ports_expected'),('altered primitive slope','P_to_C_matrix_e_h_in_a_f')]:
        bad=deepcopy(result['roots'][0]);bad[field][1 if field=='Q_ports_expected' else 0][0]+=1
        reject(label,lambda:check_root(bad,by_name[bad['name']]))
    bad=deepcopy(result['roots'][0]);bad['Q_conjugator']=[]
    reject('lost necessary peripheral conjugation',lambda:check_root(bad,by_name[bad['name']]))
    bad=deepcopy(structure);bad['results'][2]['cuts'][0]['fragments'][0]['arc_edge_classes'][0].append(9)
    reject('wrong physical fragment arc',lambda:check_small_audit(bad,smallplan))
    reject('missing graph certificate cannot transfer rank',lambda:graph_replay(result['diagrams']['13a_233']['selected'],result['diagrams']['13a_234']['selected'],result['pairs'][0]['certificate']))
    # The valid rank2 three-piece comparison is the A(2)=(2,4)-torus-link
    # exterior with its two genuine meridians k and s*k: their fiber
    # intersections are1, and filling either leaves infinite cyclic group.
    central_relation=commutator([1],[2,2])
    assert subst(central_relation,{1:[],2:[2]})==[]
    assert subst(central_relation,{1:[-2],2:[2]})==[]
    controls={'three_piece_count_alone_is_insufficient':True,'central_A2_genuine_meridian_intersections':[1,1],
        'central_A2_fill_k':{'k_image':[],'remaining_generators':['s'],'remaining_relators':[]},
        'central_A2_fill_sk':{'k_image':['s^-1'],'remaining_generators':['s'],'remaining_relators':[]},
        'old_four_piece_positive_roots':['13a_129','13a_135'],
        'S3_C2_C3_noncommuting':True}
    files=['verify_three_piece_covers.py', 'three_piece_cover_inputs.json', 'three_piece_cover_results.json', 'remaining_cover_structure_inputs.json', 'remaining_cover_structure_results.json', 'compact_cover_results.json', 'parabolic_panel.json', 'joint_a5_results.json', 'verify_compact_covers.py', 'verify_port_assemblies.py', 'verify_mutation_covers.py', 'verify_four_piece_covers.py', 'verify_double_transpositions.py', 'parabolic_lifts.py', 'THREE_PIECE_DOUBLE_COVERS_SHARP_FOUR.tex']
    receipt={'status':'PASS independent three-root peripheral and upper4 inputs; five transfers unproved','small_source_audit':small,
        'roots':roots,'diagrams':diagrams,'mutation_pairs':pairs,'unproved_transfers':unproved,'upper4':upper,'controls':controls,'corruptions_rejected':rejected,
        'assurance':'All finite source traces, selected ports, complete amalgam normal forms, slopes, S3 filling quotients and eight uppers independently checked without producer imports. Five proposed transfers have no certificate and transfer no rank conclusion; all five unequal-determinant pairs cannot have homeomorphic covers. Legacy source parser used only for arc-index transport. Canonical JSJ, realization, two-bridge unknot/filling fact and BW05 case6 exclusion are written uniform arguments and explicit classical imports, not proved by this finite receipt.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'three_piece_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'roots':len(roots),'mutation_pairs':len(pairs),'upper4':len(upper),'corruptions':len(rejected)}))


if __name__=='__main__':main()
