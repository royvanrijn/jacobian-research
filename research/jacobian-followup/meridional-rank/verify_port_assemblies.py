#!/usr/bin/env python3
"""Independent fragment geometry, exact moves, ports and full-core assembly."""
from collections import Counter,defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from verify_compact_covers import replay,clean,replace,backwards,reduce_word,keys,abelianization
from verify_double_transpositions import geometric_crossings

HERE=Path(__file__).resolve().parent


def arc_classes(pd):
    adjacent=defaultdict(set)
    for row in pd:
        for e in row:adjacent[e]
        adjacent[row[1]].add(row[3]);adjacent[row[3]].add(row[1])
    remaining=set(adjacent);groups=[]
    while remaining:
        reached={min(remaining)};queue=list(reached)
        for e in queue:
            for f in adjacent[e]:
                if f not in reached:reached.add(f);queue.append(f)
        remaining-=reached;groups.append(sorted(reached))
    return sorted(groups)


def check_fragment(data,actual_pd):
    assert data['pd']==actual_pd
    groups=arc_classes(actual_pd);assert groups==data['arc_edge_classes']
    counts=Counter(e for row in actual_pd for e in row)
    assert all(v in (1,2) for v in counts.values())
    boundary=sorted(e for e,v in counts.items() if v==1);assert boundary==data['boundary_labels'] and len(boundary)==4
    arcs={e:i+1 for i,g in enumerate(groups) for e in g}
    triples=[[arcs[row[1]],arcs[row[0]],arcs[row[2]]] for row in actual_pd]
    assert triples==data['crossing_arcs']
    rels=[reduce_word([x for x in (o,-a,o,-b) if abs(x)!=1]) for o,a,b in triples]
    assert rels==data['raw_relators'];reduction=data['formal_core_reduction']
    assert reduction['initial']=={'generators':list(range(2,len(groups)+1)),'relators':clean(rels)}
    gens,rels=replay(reduction);assert len(gens)==2 and len(rels)==1
    return groups


def full_arc_names(pd):
    triples=geometric_crossings(pd);names={}
    for row,(o,a,b) in zip(pd,triples):
        for edge,name in zip(row,(a,o,b,o)):
            assert edge not in names or names[edge]==name;names[edge]=name
    groups=arc_classes(pd);assert len(groups)==13
    assert len({next(iter({names[e] for e in g})) for g in groups})==13
    assert all(len({names[e] for e in g})==1 for g in groups)
    return names,triples


def check_join(I,O,pd):
    # The unnormalized fragment arc classes glue across four ports to the full arcs.
    pieces=[set(g) for g in I['arc_edge_classes']]+[set(g) for g in O['arc_edge_classes']]
    assert len(pieces)==17
    for edge in I['boundary_labels']:
        incident=[g for g in pieces if edge in g];assert len(incident)==2
        a,b=incident;pieces.remove(a);pieces.remove(b);pieces.append(a|b)
    assert len(pieces)==13 and sorted(sorted(g) for g in pieces)==arc_classes(pd)
    return {'fragment_arc_classes':17,'independent_port_identifications':4,'full_arc_classes':13}


def check_darts(A,B,certificate):
    orientation=certificate['cyclic_orientation'];assert orientation==-1
    mapping={int(c):v for c,v in certificate['crossing_rotations'].items()}
    assert set(mapping)==set(range(6)) and sorted(v[0] for v in mapping.values())==list(range(6))
    edge_image={};ports=[]
    countA=Counter(e for row in A for e in row);countB=Counter(e for row in B for e in row)
    for c,row in enumerate(A):
        d,k=mapping[c];assert k in (0,2)
        for h,e in enumerate(row):
            f=B[d][(orientation*h+k)%4]
            assert e not in edge_image or edge_image[e]==f;edge_image[e]=f
            assert countA[e]==countB[f]
            if countA[e]==1:ports.append([e,f])
    assert len(set(edge_image.values()))==len(edge_image) and ports==certificate['boundary_edge_map']


def check_assembly(row,I,O,pd):
    assert I['name']==O['name']==row['name'];check_join(I,O,pd)
    pieces=[]
    for data,gens in ((I,[1,2]),(O,[3,4])):
        reduction=data['formal_core_reduction'];rename={g:[h] for g,h in zip(reduction['final']['generators'],gens)}
        originals={int(g):replace(w,rename) for g,w in reduction['original_generators_in_final'].items()};originals[1]=[]
        pieces.append({'generators':gens,'relators':[replace(w,rename) for w in reduction['final']['relators']],
            'original_arc_words':originals})
    assert row['pieces']==[{**p,'original_arc_words':{str(k):v for k,v in p['original_arc_words'].items()}} for p in pieces]
    ports=[]
    for edge in I['boundary_labels']:
        i=next(j+1 for j,s in enumerate(I['arc_edge_classes']) if edge in s)
        o=next(j+1 for j,s in enumerate(O['arc_edge_classes']) if edge in s)
        u=pieces[0]['original_arc_words'][i];v=pieces[1]['original_arc_words'][o]
        ports.append({'edge':edge,'inner_arc':i,'outer_arc':o,'inner_word':u,'outer_word':v,
            'relation':reduce_word(u+[5]+backwards(v))})
    assert ports==row['ports']
    relations=pieces[0]['relators']+pieces[1]['relators']+[p['relation'] for p in ports]
    reduction=row['assembly_reduction']
    assert reduction['initial']=={'generators':[1,2,3,4,5],'relators':clean(relations)}
    gens,rels=replay(reduction);assert len(gens)==len(rels)==3
    return abelianization(gens,rels)


def check_transport(row,O,pd):
    names,triples=full_arc_names(pd)
    mapped=[next(iter({names[e] for e in group})) for group in O['arc_edge_classes']]
    assert mapped==row['outside_to_full_arcs'] and mapped[0]==row['outside_base_full_arc']==5
    original=keys(O['formal_core_reduction']['original_generators_in_final']);original[1]=[]
    initial=[{'outside_arc':i+1,'full_arc':arc,'word':original[i+1]} for i,arc in enumerate(mapped)]
    assert row['initial']==initial;known={}
    for entry in initial:known.setdefault(entry['full_arc'],entry['word'])
    for step in row['trace']:
        o,a,b=triples[step['crossing']];source,target=step['under_known'],step['under_new']
        assert {source,target}=={a,b} and o==step['over']
        assert o in known and source in known and target not in known
        known[target]=reduce_word(known[o]+backwards(known[source])+known[o])
    assert known==keys(row['known_words'])
    assert all(not (o in known and ((a in known)!=(b in known))) for o,a,b in triples)
    missing=sorted(set(range(13))-set(known));assert missing==row['unresolved_full_arcs']==[0,1,2,3,4]
    assert row['status']=='BOUNDARY_GENERATION_STALLED'
    return missing


def main():
    source=json.loads((HERE/'mutation_cover_results.json').read_text())['diagrams']
    inner=json.loads((HERE/'inner_tangle_results.json').read_text())
    outer=json.loads((HERE/'outer_tangle_results.json').read_text())
    assembly=json.loads((HERE/'port_assembly_results.json').read_text())
    transport=json.loads((HERE/'outer_tangle_transport.json').read_text())
    names=['13a_129','13a_135']
    assert [r['name'] for r in inner['fragments']]==[r['name'] for r in outer['fragments']]==names
    assert len(outer['outside_dart_isomorphisms'])==1
    check_darts(outer['fragments'][0]['pd'],outer['fragments'][1]['pd'],outer['outside_dart_isomorphisms'][0])
    targets=[]
    for k,name in enumerate(names):
        pd=source[name]['pd'];I=inner['fragments'][k];O=outer['fragments'][k];A=assembly['targets'][k]
        check_fragment(I,pd[:7]);check_fragment(O,pd[7:])
        assert I['boundary_labels']==O['boundary_labels']==[12,15,17,26]
        diagonal=check_assembly(A,I,O,pd);assert diagonal==[1,1,[391,401][k]]
        missing=check_transport(transport['results'][k],O,pd)
        targets.append({'name':name,'inner_relator_length':11,'outer_relator_length':7,
            'arc_join':check_join(I,O,pd),'assembly_abelianization_smith':diagonal,
            'port_words':A['ports'],'assembly_final':A['assembly_reduction']['final'],
            'unresolved_arcs_after_outer_generation':missing,'actual_rank_interval':[2,3]})
    assert inner['fragments'][0]['formal_core_reduction']['final']==inner['fragments'][1]['formal_core_reduction']['final']
    assert assembly['targets'][0]['pieces'][0]['relators']==assembly['targets'][1]['pieces'][0]['relators']
    assert assembly['targets'][0]['pieces'][1]['relators']==assembly['targets'][1]['pieces'][1]['relators']
    assert [len(w) for w in assembly['targets'][0]['pieces'][0]['relators']]==[11]
    assert [len(w) for w in assembly['targets'][0]['pieces'][1]['relators']]==[7]
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    I=inner['fragments'][0];O=outer['fragments'][0];pd=source[names[0]]['pd']
    bad=deepcopy(assembly['targets'][0]);bad['ports'][0]['relation'][0]*=-1
    reject('altered port relation',lambda:check_assembly(bad,I,O,pd))
    bad=deepcopy(O);bad['arc_edge_classes'][0].append(17)
    reject('wrong fragment arc identification',lambda:check_fragment(bad,pd[7:]))
    bad=deepcopy(outer['outside_dart_isomorphisms'][0]);bad['boundary_edge_map'][0][1]=26
    reject('wrong reflected port transport',lambda:check_darts(O['pd'],outer['fragments'][1]['pd'],bad))
    bad=deepcopy(transport['results'][0]);bad['unresolved_full_arcs']=[]
    reject('false closed2generator completion',lambda:check_transport(bad,O,pd))
    files=['verify_port_assemblies.py','verify_compact_covers.py','compact_cover_verification.json',
        'inner_tangle_plan.json','inner_tangle_results.json','outer_tangle_plan.json','outer_tangle_results.json',
        'outer_tangle_transport_plan.json','outer_tangle_transport.json','port_assembly_plan.json','port_assembly_results.json',
        'mutation_cover_results.json','verify_double_transpositions.py','parabolic_lifts.py']
    for report in (inner,outer,assembly,transport):
        for path,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest,path
    receipt={'status':'PASS independent four-port core assemblies','targets':targets,'corruptions_rejected':rejected,
        'common_formal_inner_and_outer_groups':True,'outer_ribbon_map_orientation':-1,
        'rank3_proved':False,'rank2_proved':False,'new_target_closures':0,
        'assurance':'Exact normalized core-group assembly by base change u_i*z and four generator identifications; all fragment and assembly moves replay. Helper arithmetic comes from independent compact verifier, not producer. Shared legacy arc-label parser disclosed. No injectivity/peripheral/Seifert/JSJ or filling-slope theorem is inferred. Uniform elementary assembly proof must accompany the finite inputs.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'port_assembly_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'targets':2,'fragment_presentations':4,'assembly_presentations':2,
        'corruptions':len(rejected),'actual_rank_interval':[2,3]}))


if __name__=='__main__':main()
