#!/usr/bin/env python3
"""Replay exact Nielsen/Tietze certificates; --snappy rechecks PD exteriors."""
import argparse
from copy import deepcopy
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPECTED = ['13a_373','13a_375','13a_391','13a_392',
            '13a_3264','13a_3266','13a_3275','13a_3276'] + ['13a_304', '13a_362', '13a_623', '13a_629', '13a_3232', '13a_3259', '13a_3374', '13a_3375', '13n_834', '13n_838', '13n_839', '13n_991', '13n_992', '13n_993', '13n_1367', '13n_1370', '13n_1371', '13n_1374', '13n_1375', '13n_1376', '13n_3814', '13n_3840', '13n_3911', '13n_3914']


def reduce_word(word):
    answer=[]
    for a in word:
        assert isinstance(a,int) and a
        if answer and answer[-1]==-a:answer.pop()
        else:answer.append(a)
    return answer


def inverse(word):return [-a for a in reversed(word)]


def cyclic_reduce(word):
    word=reduce_word(word)
    while len(word)>1 and word[0]==-word[-1]:word=word[1:-1]
    return word


def substitute(word,generator,replacement):
    other=inverse(replacement)
    return reduce_word([a for b in word
                        for a in (replacement if b==generator
                                  else other if b==-generator else [b])])


def expand(word,basis):
    return reduce_word([a for b in word
                        for a in (basis[b-1] if b>0 else inverse(basis[-b-1]))])


def nielsen_controls():
    for a in (1,-1):
        for b in (2,-2):
            sign=-1 if a*b<0 else 1
            factor=[2*sign]
            new_basis=[[1]+factor if a>0 else factor+[1],[2],[3]]
            inverse_map=[[1,-2*sign] if a>0 else [-2*sign,1],[2],[3]]
            assert all(expand(w,new_basis)==[i+1] for i,w in enumerate(inverse_map))
            assert all(expand(w,inverse_map)==[i+1] for i,w in enumerate(new_basis))
    return {'handle_slide_sign_cases':4,'free_inverse_compositions':8}


def replay(row):
    n=row['original_generator_count'];original_n=n
    relators=deepcopy(row['original_relators'])
    assert all(all(1<=abs(a)<=n for a in w) for w in relators)
    original_in_current=[[i] for i in range(1,n+1)]
    current_in_original=[[i] for i in range(1,n+1)]
    moves=list(row['word_moves']);position=0
    counts={'handle_slides':0,'inversions':0,'introductions':0,'unit_deletions':0}
    deletion_lengths=[]
    def change(generator,replacement):
        nonlocal relators,original_in_current
        relators=[cyclic_reduce(substitute(w,generator,replacement)) for w in relators]
        original_in_current=[substitute(w,generator,replacement) for w in original_in_current]
    while position<len(moves):
        a=moves[position];position+=1
        if a>n:
            assert a==n+1,'new generator index is not next'
            end=moves.index(a,position);definition=moves[position:end];position=end+1
            assert all(1<=abs(b)<=n for b in definition)
            relators.append([a]+inverse(definition))
            current_in_original.append(expand(definition,current_in_original))
            n+=1;counts['introductions']+=1
            continue
        assert 1<=abs(a)<=n and position<len(moves)
        b=moves[position];position+=1
        assert 1<=abs(b)<=n
        if a==b:
            assert a>0,'negative deletion index'
            units=[i for i,w in enumerate(relators) if cyclic_reduce(w) in ([a],[-a])]
            assert units,'deleted generator has no current unit relator'
            deletion_lengths.append(len(relators[units[0]]))
            relators.pop(units[0]);change(a,[])
            if a!=n:
                change(n,[a]);current_in_original[a-1]=current_in_original[-1]
            current_in_original.pop();n-=1;counts['unit_deletions']+=1
        elif a==-b:
            assert a>0,'negative inversion index'
            change(a,[-a]);current_in_original[a-1]=inverse(current_in_original[a-1])
            counts['inversions']+=1
        else:
            index=abs(a);other=abs(b)
            exponent=-1 if a*b<0 else 1
            A=current_in_original[index-1]
            B=current_in_original[other-1]
            if exponent<0:B=inverse(B)
            current_in_original[index-1]=reduce_word(A+B if a>0 else B+A)
            inverse_letter=-exponent*other
            replacement=[index,inverse_letter] if a>0 else [inverse_letter,index]
            change(index,replacement);counts['handle_slides']+=1
    assert n==3
    assert len(original_in_current)==original_n
    assert all(all(1<=abs(a)<=3 for a in w) for w in relators+original_in_current)
    if 'final_relators' in row:
        assert relators==row['final_relators']
        assert original_in_current==row['original_generators_in_final_three']
        assert current_in_original==row['final_three_in_original_generators']
    return {'name':row['name'],'final_generator_count':n,
            'final_relators':relators,
            'original_generators_in_final_three':original_in_current,
            'final_three_in_original_generators':current_in_original,
            'move_counts':counts,'deletion_witness_lengths':deletion_lengths}


def combinatorial_data(text):
    """Discard floating-point solution and shape data from SnapPea text."""
    lines=text.splitlines()
    # A SnapPea triangulation stores each tetrahedron in eight nonblank lines:
    # neighbors, face permutations, cusp indices, four peripheral rows, shape.
    blank_free=[line.split() for line in lines if line.strip()]
    cursor=next(i for i,row in enumerate(blank_free)
                if len(row)==2 and all(t.lstrip('-').isdigit() for t in row)
                and row[0]=='1' and row[1]=='0')
    assert blank_free[cursor+1][0]=='torus'
    assert all(Fraction(a)==0 for a in blank_free[cursor+1][1:3]),'filled snapshot'
    count=int(blank_free[cursor+2][0]);start=cursor+3
    tetrahedra=[]
    for i in range(count):
        block=blank_free[start+8*i:start+8*(i+1)]
        assert len(block)==8 and len(block[0])==len(block[1])==len(block[2])==4
        assert all(len(r)==16 for r in block[3:7])
        tetrahedra.append(block[:7])
    return tetrahedra


def snappy_replay(rows):
    import random
    import warnings
    warnings.filterwarnings('ignore',message='Plink failed to import tkinter.*')
    import snappy
    import spherogram
    from importlib.metadata import version
    from audit_diagram_identity import normalize,pd_walk
    assert snappy.__version__=='3.3.2'
    assert version('spherogram')=='2.4.1'
    reports=[]
    for row in rows:
        lookup='K'+row['name'].replace('_','')
        census_pd=spherogram.Link(lookup).PD_code()
        assert normalize(pd_walk(row['pd']),True)[0]==normalize(pd_walk(census_pd),True)[0]
        seed=row['constructor_seed'];random.seed(seed);snappy.set_rand_seed(seed)
        manifold=spherogram.Link(row['pd']).exterior()
        assert manifold.num_cusps()==1
        assert combinatorial_data(manifold._to_string())==combinatorial_data(row['snappea_triangulation'])
        frozen=snappy.Triangulation(row['snappea_triangulation'])
        assert frozen.num_tetrahedra()==row['tetrahedra']
        assert frozen.cusp_info('filling')==[(0.0,0.0)]
        assert frozen.triangulation_isosig()==row['triangulation_isosig']
        raw=frozen.fundamental_group(False,False,False)
        assert raw.num_generators()==row['original_generator_count']
        assert raw.relators(as_int_list=True)==row['original_relators']
        simplified=frozen.fundamental_group()
        assert simplified.num_generators()==3
        assert simplified._word_moves()==row['word_moves']
        assert simplified.generators_in_originals(as_int_list=True)==row['final_three_in_original_generators']
        reports.append({'name':row['name'],'literal_PD_to_snapshot_combinatorics':True,
                        'named_census_planar_map':True,'original_cell_presentation':True,
                        'numerical_geometry_used':False})
    return {'snappy_version':snappy.__version__,'spherogram_version':version('spherogram'),'cases':reports}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--snappy',action='store_true')
    args=parser.parse_args()
    controls=nielsen_controls()
    data=json.loads((HERE/'rank_three_knot_group_inputs.json').read_text())
    assert [r['name'] for r in data['cases']]==EXPECTED
    frozen=json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    by_name={r['name']:r for r in frozen['cases']}
    added=json.loads((HERE/'a6_affine_inputs.json').read_text())
    by_name.update({r['name']:r for r in added['cases']})
    for row in data['cases']:
        assert row['pd']==by_name[row['name']]['pd']
        assert len(combinatorial_data(row['snappea_triangulation']))==row['tetrahedra']
    results=[replay(row) for row in data['cases']]
    corrupt=deepcopy(data['cases'][0]);corrupt['original_relators']=[]
    try:replay(corrupt)
    except AssertionError:pass
    else:raise AssertionError('missing deletion-witness relators accepted')
    topology=snappy_replay(data['cases']) if args.snappy else None
    sources={'sources/snappy-3.3.2-fundamental_group.c':
             '8b5e62b703d92672100feef7a5d2dcc3406b01929c4e2faf04428263e3b29191',
             'sources/snappy-3.3.2-fundamental_group.pyx':
             '2537075b91b1b3d606fd7574fce4fa6c63d899cc2b4badb1827dcfd95ea8b1f5',
             'sources/snappy-3.3.2-init.py':
             'ec32a0e700665926426bc3a178680463f40c5a50695776b653caa99c2c112277'}
    for name,digest in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    files=['RANK_THREE_KNOT_GROUP_UPPERS.tex','verify_rank_three_knot_groups.py',
           'rank_three_knot_group_inputs.json','seven_sheet_cover_inputs.json',
           'a6_affine_inputs.json',
           'audit_diagram_identity.py']
    receipt={'status':'PASS','cases':results,'corrupted_cell_relators_rejected':True,
             'free_group_Nielsen_controls':controls,
             'primary_software_source_sha256':sources,
             'snappy_reconstructed':args.snappy,'topology_replay':topology,
             'scope':'All retained algebra moves independently checked as Nielsen transformations, generator definitions or deletions with a literal unit relator. Original cell presentation and PD-to-exterior construction remain trusted exact SnapPy combinatorial algorithms; optional --snappy independently reruns those bindings and named-census planar-map checks, with no numerical geometry acceptance.',
             'input_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}}
    (HERE/'rank_three_knot_group_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','knots':len(results),
                      'snappy_reconstructed':args.snappy,
                      'unit_deletions':[r['move_counts']['unit_deletions'] for r in results]},indent=2))


if __name__=='__main__':main()
