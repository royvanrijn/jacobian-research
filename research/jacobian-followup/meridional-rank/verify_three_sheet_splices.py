#!/usr/bin/env python3
"""Replay sharp-four certificates from independent planar arc traces and amalgams."""
from collections import Counter
import hashlib
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram  # source-index transport only
from verify_integer_lifts import replay
from verify_three_sheet_structure import evaluate,normalized,invert,rotations,canonical

HERE=Path(__file__).resolve().parent


def power(word,n):return word*n if n>=0 else invert(word)*(-n)


def planar_arcs(pd):
    ends={}
    for c,crossing in enumerate(pd):
        for h,edge in enumerate(crossing):ends.setdefault(edge,[]).append((c,h))
    assert all(len(v)==2 for v in ends.values());over={};incoming={};paths=[]
    for arc in range(len(pd)):
        current=(arc,2);path=[]
        for _ in range(2*len(pd)):
            pair=ends[pd[current[0]][current[1]]];other=pair[1] if pair[0]==current else pair[0]
            c,h=other;path.append([c,h])
            if h%2==0:
                assert h==0 and c not in incoming;incoming[c]=arc;break
            assert c not in over;over[c]=arc;current=(c,(h+2)%4)
        else:raise AssertionError('arc tracing did not terminate')
        paths.append(path)
    assert len(over)==len(incoming)==len(pd)
    return [(over[c],incoming[c],c) for c in range(len(pd))],paths


def branch_words(crossings,labels,new_to_old,tree):
    roots=[(0,1),(0,2),(1,2)];result=[]
    for o,a,b in crossings:
        for start in range(3):
            sheet=start;word=[]
            for arc in [o,a,o,b]:
                p,q=roots[labels[arc]];name=new_to_old[arc]+1
                if sheet==p:
                    if name not in tree:word.append(name)
                    sheet=q
                elif sheet==q:
                    if name not in tree:word.append(-name)
                    sheet=p
            assert sheet==start;result.append(word)
    return result


def check_splice_trace(row):
    words=normalized(row['initial_relators']);images={1:[1],2:[2]}
    for step in row['steps']:
        before=sum(map(len,words))
        if step['operation']=='nielsen':
            forward={int(k):v for k,v in step['forward'].items()};backward={int(k):v for k,v in step['inverse'].items()}
            assert set(forward)==set(backward)=={1,2}
            for g in [1,2]:assert evaluate(forward[g],backward)==evaluate(backward[g],forward)==[g]
            words=normalized([evaluate(w,forward) for w in words]);images={g:evaluate(w,forward) for g,w in images.items()}
        else:
            assert step['operation']=='relator_shortening'
            r=step['used_relator'];target=step['target_relator'];signed=step['used_signed_rotation'];w=step['target_rotation'];k=step['pattern_length']
            assert r in words and target in words and r!=target
            assert signed in rotations(r)+rotations(invert(r)) and w in rotations(target)
            assert k>len(r)//2 and w[:k]==signed[:k]
            new=words[:];new.remove(target);new.append(invert(signed[k:])+w[k:]);words=normalized(new)
        assert step['before_measure']==before and step['after_measure']==sum(map(len,words))
    assert words==row['final_relators']==row['target_relators']
    return images


def amalgam_check(kind,e,target):
    # Generators a,b,d,e are numbered1,2,3,4. Edge relations solve d,e.
    dword=[-2]*3+[1] if kind=='A' else [2]*2+[1]
    eword=power([1],3*e)+[2]
    substitution={1:[1],2:[2],3:dword,4:eword}
    first_vertex=[3,3]+[-1]*3
    if kind=='A':
        second_vertex=[4,4]+[-2]*3
        edges=[[3,-1]+[2]*3,[1]*3+power([4,-2],-e)]
        edge_matrix=[[0,e],[-1,0]]
    else:
        second_vertex=[2,2]+[-4]*3
        edges=[[3,-1]+[-2]*2,[1]*3+power([2,-4],e)]
        edge_matrix=[[0,-e],[1,0]]
    assert all(evaluate(w,substitution)==[] for w in edges)
    assert normalized([evaluate(first_vertex,substitution),evaluate(second_vertex,substitution)])==target
    assert abs(edge_matrix[0][0]*edge_matrix[1][1]-edge_matrix[0][1]*edge_matrix[1][0])==1
    return {'first_vertex':'<d,a | d^2=a^3>','second_vertex':'<e,b | e^2=b^3>' if kind=='A' else '<b,e | b^2=e^3>',
            'eliminated_vertex_generator_words':{3:dword,4:eword},'edge_basis_matrix':edge_matrix,
            'edge_injectivity':'Each (meridian,fibre) subgroup is Z^2: quotient by fibre is C2*C3 and meridian powers have nonempty reduced form; fibre has infinite order by abelianization. Amalgam normal form embeds each nonabelian vertex group.',
            'cover_group_rank':2}


def torus_normal(word):
    fibre=0;syllables=[]
    for a in word:
        g=abs(a);mod=2 if g==1 else 3;value=1 if a>0 else -1
        if syllables and syllables[-1][0]==g:value+=syllables.pop()[1]
        carry,remainder=divmod(value,mod);fibre+=carry
        if remainder:syllables.append((g,remainder))
    return fibre,tuple(syllables)


def abelian_cover_check(row):
    N=row['abelianization_order'];r=row['generator_shifts'][0];rows=[]
    source=next(x for x in json.loads((HERE/'three_sheet_structural_results.json').read_text())['results'] if x['name']==row['name'])
    for w in source['final_relators']:
        for start in range(N):
            vertex=start;full=[0]*(2*N)
            for letter in w:
                g=abs(letter);step=r if g==1 else 1
                destination=(vertex+(step if letter>0 else -step))%N
                initial=vertex if letter>0 else destination
                full[(g-1)*N+initial]+=1 if letter>0 else -1;vertex=destination
            assert vertex==start
            rows.append(full[:N]+[full[-1]])  # kill y edges0..N-2 of its spanning tree
    assert rows==row['cyclic_cover_matrix']
    diagonal=replay(rows,row['smith']);assert diagonal==[1]*(N+1)
    assert row['H1_free_rank']==0 and not row['H1_torsion']
    return {'cyclic_degree':N,'universal_abelian_cover_H1':'0','commutator_subgroup':'perfect','all_solvable_quotients':'cyclic'}


def main():
    report=json.loads((HERE/'three_sheet_splice_results.json').read_text());assert report['status']=='COMPLETE'
    structure=json.loads((HERE/'three_sheet_structure_verification.json').read_text());assert structure['status']=='PASS'
    provenance=json.loads((HERE/'splice_knot_provenance.json').read_text());assert provenance['status']=='PASS exact named PD provenance'
    metabelian=json.loads((HERE/'three_sheet_metabelian_results.json').read_text());assert metabelian['status']=='COMPLETE'
    for receipt in [report,structure,provenance,metabelian]:
        for p,h in receipt['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    original={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    upper={r['name']:r for r in json.loads((HERE/'small_symmetric_full_panel_results.json').read_text())['results']}
    # Uniform normal-form proof is in the note; these finite checks catch signs.
    assert torus_normal([1,1])==torus_normal([2,2,2])==(1,())
    assert torus_normal([1,2,-1,-2])[1]
    mu=[1,-2];fibre=[1,1];assert torus_normal(mu+fibre+invert(mu)+invert(fibre))==(0,())
    for n in range(-8,9):
        for m in range(-3,4):
            normal=torus_normal(power(mu,n)+power(fibre,m))
            assert (normal==(0,()))==(n==m==0)
    results=[]
    for row in report['results']:
        name=row['name'];old=original[name];pd=panel[name]['pd'];cr,paths=planar_arcs(pd)
        source_crossings=oriented_diagram(pd);new_to_old=[c[2] for c in source_crossings];assert sorted(new_to_old)==list(range(13))
        old_to_new={a:i for i,a in enumerate(new_to_old)};labels=[old['arc_labels'][a] for a in new_to_old]
        assert len(set(labels))>=2  # distinct transpositions generate S3
        words=branch_words(cr,labels,new_to_old,old['tree_generators'])
        assert Counter(tuple(canonical(w)) for w in words)==Counter(tuple(canonical(w)) for w in old['relators'])
        seeds=[old_to_new[s] for s in upper[name]['seeds']];assert len(set(seeds))==4
        seen=set(seeds);trace=[]
        while True:
            advance=next(((o,a,b) if a in seen else (o,b,a) for o,a,b in cr if o in seen and ((a in seen)^(b in seen))),None)
            if advance is None:break
            o,a,b=advance;seen.add(b);trace.append([o,a,b])
        assert len(seen)==13
        images=check_splice_trace(row);amalgam=amalgam_check(row['model'],row['epsilon'],row['target_relators'])
        homology=abelian_cover_check(next(r for r in metabelian['results'] if r['name']==name))
        assert next(r for r in provenance['results'] if r['name']==name)['literal_PD_matches_after_adding_one']
        results.append({'name':name,'independent_planar_arc_paths':paths,'new_to_source_arc_indices':new_to_old,
                        'S3_arc_labels':labels,'four_bridge_seed_arcs':seeds,'four_bridge_seed_trace':trace,
                        'branched_word_multiset_matches':True,'splice_basis_transport':images,'model':row['model'],'epsilon':row['epsilon'],
                        'amalgam':amalgam,'abelian_cover':homology,'cover_group_rank':2,
                        'meridional_lower':4,'bridge_upper':4,'mr_equals_bridge':4})
    assert [r['name'] for r in results]==['13a_3145','13a_3147','13a_3148','13a_3149']
    paths=['verify_three_sheet_splices.py','three_sheet_splice_plan.json','three_sheet_splice_results.json',
           'three_sheet_structure_verification.json','three_sheet_structural_results.json','three_sheet_nonabelian_results.json',
           'three_sheet_metabelian_results.json','splice_knot_provenance.json','joint_a5_panel.json',
           'small_symmetric_full_panel_results.json','parabolic_lifts.py','verify_three_sheet_structure.py','verify_integer_lifts.py']
    receipt={'status':'PASS exact sharp-four evidence','sharp_named_knots':4,'results':results,
             'total_three_sheet_covers':26,'proved_cyclic_covers':22,'proved_rank2_covers':4,'unresolved_cover_groups':0,
             'proof_assurance':'Written uniform peripheral-Z2 injection and free-product-with-amalgamation normal-form arguments, unformalized and externally unreviewed. Exact presentation, named PD, upper seeds, branch words, peripheral elimination and homology replay pass. General meridional branched-cover rank theorem and published Wirtinger-number/bridge equality are explicit dependencies.',
             'priority_scope':'Four named sharp certificates newly retained here; targets appear in the fetched external open list and pinned[3,4] data. Worldwide previously-unresolved status and known-family/priority exclusion are not established. No full MRC solution.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'three_sheet_splice_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'sharp_knots':[r['name'] for r in results],
                      'cover_group_ranks':[r['cover_group_rank'] for r in results],'unresolved_cover_groups':0},indent=2))


if __name__=='__main__':main()
