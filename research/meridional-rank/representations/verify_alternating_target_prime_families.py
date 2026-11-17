#!/usr/bin/env python3
"""Exact controls for the written alternating-target prime-family proof.

Default replay uses only the Python standard library.  It neither enumerates
large alternating groups nor replaces the quantified mathematical proof by
sampling.  Optional --topology independently reconstructs the literal PD.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
from itertools import combinations, product
import json
from math import factorial, gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'alternating_target_prime_families_inputs.json'
RECEIPT = HERE / 'alternating_target_prime_families_verification.json'
TEX = HERE / 'ALTERNATING_TARGET_PRIME_FAMILIES.tex'

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))

def inverse(a):
    out = [0] * len(a)
    for i, j in enumerate(a): out[j] = i
    return tuple(out)

def cyc(N, *cycles):
    p = list(range(N))
    for cycle in cycles:
        for a, b in zip(cycle, cycle[1:] + cycle[:1]): p[a] = b
    return tuple(p)

def cycle_lengths(p, include_fixed=False):
    unused = set(range(len(p))); ans = []
    while unused:
        a = min(unused); q = a; count = 0
        while q in unused:
            unused.remove(q); q = p[q]; count += 1
        if include_fixed or count > 1: ans.append(count)
    return sorted(ans)

def conjugate(a, b):
    return compose(compose(a, b), inverse(a))

def evaluate(word, generators):
    p = tuple(range(len(generators[0])))
    for i in word: p = compose(p, generators[i])
    return p

def hurwitz(word, initial):
    cur = list(initial)
    for k in word:
        i = abs(k) - 1
        require(0 <= i < len(cur)-1, 'braid generator range')
        u, v = cur[i:i+2]
        cur[i:i+2] = [v, conjugate(v,u)] if k > 0 else [conjugate(u,v),u]
    return cur

def family_word(beta, b, n):
    return (beta + [j for j in range(4,b) for _ in range(3)]) * n

def family_seed(quartet, b):
    N = 2*b
    return [tuple(x)+tuple(range(8,N)) for x in quartet] + [cyc(N,(1,2*j),(6,2*j+1)) for j in range(4,b)]

def strand_permutation(word, b):
    tags = list(range(b))
    for k in word:
        i = abs(k)-1; tags[i], tags[i+1] = tags[i+1], tags[i]
    return tuple(tags)

def seesaws(word, b):
    out = []
    for i in range(2,b):
        indices = [abs(k) for k in word if abs(k) in (i-1,i)]
        out.append(sum(j == 0 or indices[j] != indices[j-1] for j in range(len(indices))))
    return out

def validate_domain(b, n, four=False):
    require(isinstance(b,int) and isinstance(n,int), 'integer parameters')
    require((b == 4 and n >= 1) if four else (b >= 5 and n >= 2), 'prime-family parameter domain')
    require(gcd(n,b) == 1, 'coprime parameters required for a knot')

def components(adjacency):
    unused = set(adjacency); out = []
    while unused:
        block = {min(unused)}; todo = list(block)
        for x in todo:
            for y in adjacency[x] - block: block.add(y); todo.append(y)
        unused -= block; out.append(sorted(block))
    return sorted(out)

def verify_literal(data, beta, quartet):
    b,n = data['b'],data['n']; validate_domain(b,n)
    word = family_word(beta,b,n)
    require(data['literal_braid_word'] == word, 'literal word binding')
    pd = data['pd']; m = len(pd); N = 2*b
    require(m == len(word) == 40, 'literal crossing count')
    require(all(len(row) == 4 for row in pd), 'PD row arity')
    counts = Counter(e for row in pd for e in row)
    require(set(counts) == set(range(2*m)) and set(counts.values()) == {2}, 'PD labels occur twice')
    require(data['signed_pd'] == [[1 if k>0 else -1]+row for k,row in zip(word,pd)], 'signed PD binding')
    occurrences = defaultdict(list)
    for j,row in enumerate(pd):
        for k,e in enumerate(row): occurrences[e].append(4*j+k)
    pairing = {}
    for positions in occurrences.values():
        a,bp = positions; pairing[a]=bp; pairing[bp]=a
    face_perm = tuple(4*(pairing[h]//4)+(pairing[h]%4+1)%4 for h in range(4*m))
    faces = len(cycle_lengths(face_perm,include_fixed=True))
    require(m-2*m+faces == 2, 'spherical PD rotation system')
    link_adj = {h:{pairing[h],4*(h//4)+(h%4+2)%4} for h in range(4*m)}
    require(len(components(link_adj)) == 1, 'one component of literal PD')
    adjacency = {e:set() for e in counts}
    for a,o,c,d in pd: adjacency[o].add(d); adjacency[d].add(o)
    blocks = components(adjacency)
    require(blocks == data['arc_edge_components'] and len(blocks)==m, 'arc component convention')
    arc = {e:i for i,block in enumerate(blocks) for e in block}
    images = [tuple(p) for p in data['arc_permutations']]
    require(len(images)==m, 'one permutation per arc')
    require(all(sorted(p)==list(range(N)) and cycle_lengths(p)==[2,2] for p in images), 'double-transposition arc images')
    crossings = [[arc[o],arc[a],arc[c],1 if k>0 else -1] for k,(a,o,c,d) in zip(word,pd)]
    require(crossings == data['signed_wirtinger_crossings'], 'signed Wirtinger binding')
    for o,a,c,sign in crossings:
        over = images[o] if sign>0 else inverse(images[o])
        require(conjugate(over,images[a]) == images[c], 'signed literal Wirtinger relation')
    initial = family_seed(quartet,data['b'])
    require([list(p) for p in initial] == data['generator_permutations'], 'seed formula binding')
    cur=initial[:]; tags=list(range(b)); seed_edges=[None]*b
    wire=[None]*b; topwire=[None]*b
    for k,row in zip(word,pd):
        i=abs(k)-1; u,v=cur[i:i+2]
        if k>0:
            ports=[u,v,conjugate(v,u),v]; incoming=[0,3]; outgoing=[1,2]; new=[v,conjugate(v,u)]
        else:
            ports=[v,u,conjugate(u,v),u]; incoming=[1,0]; outgoing=[2,3]; new=[conjugate(u,v),u]
        for column,port in [(i,incoming[0]),(i+1,incoming[1])]:
            if wire[column] is None: wire[column]=row[port]; topwire[column]=row[port]
            require(wire[column]==row[port], 'literal incoming braid port wiring')
            if tags[column] is not None: seed_edges[tags[column]]=row[port]
        for edge,color in zip(row,ports): require(images[arc[edge]]==color, 'literal braid arc-color binding')
        wire[i:i+2]=[row[j] for j in outgoing]; cur[i:i+2]=new; tags[i]=tags[i+1]=None
    require(cur==initial and wire==topwire, 'color and edge closure wiring')
    require(all(e is not None for e in seed_edges), 'all seed strands touched')
    seed_arcs=[arc[e] for e in seed_edges]
    require(seed_arcs==data['seed_arcs'], 'literal five seed arcs')
    require([images[j] for j in seed_arcs]==initial, 'literal generating tuple present')
    known=set(seed_arcs); relations={tuple(row[:3]) for row in crossings}
    for o,a,c in data['seed_trace']:
        require((o,a,c) in relations or (o,c,a) in relations, 'trace uses literal crossing')
        require(o in known and a in known and c not in known, 'valid Wirtinger coloring move')
        known.add(c)
    require(known==set(range(m)) and len(data['seed_trace'])==m-b, 'five seeds color every arc')
    return {'b':b,'n':n,'crossings':m,'pd_faces':faces,'pd_euler_characteristic':2,
            'components':1,'signed_wirtinger_relations_checked':m,'arc_permutations_checked':m,
            'seed_arcs':seed_arcs,'seed_moves_checked':len(data['seed_trace']),
            'full_alternating_target_order_from_written_generation_theorem':factorial(N)//2,
            'target_order_not_obtained_by_large_group_enumeration':True}

def connected_tree(edges, N=6):
    parent=list(range(N))
    def find(x):
        while parent[x]!=x: x=parent[x]
        return x
    for a,b in edges:
        a,b=find(a),find(b)
        if a==b: return False
        parent[a]=b
    return len({find(i) for i in range(N)})==1

def canonical_ordered_tree(edges):
    # Any conjugacy must send the first edge to itself. Anchor each of its
    # two orientations, then discover vertices deterministically using edge
    # order. A connected tree makes this labeling exhaustive and unique.
    encodings=[]
    for a,b in (edges[0],edges[0][::-1]):
        labels={a:0,b:1}
        while len(labels)<6:
            before=len(labels)
            for u,v in edges:
                if u in labels and v not in labels: labels[v]=len(labels)
                if v in labels and u not in labels: labels[u]=len(labels)
            require(len(labels)>before, 'connected tree discovery')
        encodings.append(tuple(tuple(sorted((labels[u],labels[v]))) for u,v in edges))
    return min(encodings)

def symmetric_control(data, beta):
    edges=list(combinations(range(6),2)); reps=Counter(); normalized=0
    for tail in product(range(15),repeat=4):
        ordered=[edges[0]]+[edges[i] for i in tail]
        if connected_tree(ordered):
            normalized+=1; reps[canonical_ordered_tree(ordered)]+=1
    require(normalized==10368 and len(reps)==data['expected_full_conjugacy_orbits']==216, 'complete S6 tree tuple orbit count')
    require(set(reps.values())=={48}, 'normalized tree orbit multiplicities')
    word=family_word(beta,5,2); returning=[]
    for ordered in sorted(reps):
        seed=[cyc(6,edge) for edge in ordered]
        if hurwitz(word,seed)==seed: returning.append([list(e) for e in ordered])
    require(len(returning)==data['expected_returning_conjugacy_orbits']==5, 'S6 returning orbit count')
    witness=data['witness_transpositions']; seed=[cyc(6,tuple(edge)) for edge in witness]
    require(connected_tree(witness) and hurwitz(word,seed)==seed, 'literal full S6 positive witness')
    return {'scope':'all transitive ordered five-transposition tuples on six points modulo simultaneous conjugacy',
            'normalized_first_transposition_tuples':normalized,'full_conjugacy_orbits':len(reps),
            'returning_conjugacy_orbits':len(returning),'returning_representatives':returning,
            'positive_witness':witness,'arbitrary_Coxeter_absence_claim':False}

def expect_rejection(name, operation):
    try: operation()
    except (AssertionError,ValueError,IndexError,KeyError): return {'name':name,'rejected':True}
    raise AssertionError('negative control unexpectedly accepted: '+name)

def run(data):
    quartet=[tuple(p) for p in data['seed_quartet']]; beta=data['beta']
    expected=[cyc(8,(0,1),(2,3)),cyc(8,(1,4),(2,3)),cyc(8,(0,5),(1,6)),cyc(8,(1,2),(6,7))]
    require(quartet==expected, 'base quartet literal definition')
    require(len(beta)==17 and beta==[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2], 'base signed word')
    for key,target in [('ordinary_a',cyc(8,(0,1,2,3,4,5,6))),('ordinary_b',cyc(8,(0,1,7)))]:
        entry=data['base_short_words'][key]
        require(tuple(entry['permutation'])==target and evaluate(entry['word_in_quartet'],quartet)==target, 'short base generation identity: '+key)
    require(hurwitz(beta,quartet)==quartet, 'seventeen-letter base return')
    require(cycle_lengths(strand_permutation(beta,4))==[4], 'base four-cycle permutation')
    # Exact identities supporting the symbolic induction in the manuscript.
    induction=[]
    for m in (8,10,12,14):
        N=m+2; u,v=m,m+1; z=cyc(N,(1,u),(6,v))
        zp=conjugate(cyc(N,(1,0,2)),z)
        zpp=conjugate(cyc(N,(6,0,2)),z)
        require(compose(z,zp)==cyc(N,(1,u,0)), 'first induction three-cycle identity')
        require(compose(z,zpp)==cyc(N,(6,v,0)), 'second induction three-cycle identity')
        induction.append({'old_degree':m,'new_degree':N,'two_three_cycle_identities':True})
    # Independent elementary identities used in the ordinary-rank proof.
    require(conjugate(cyc(6,(0,1,2)),cyc(6,(0,2,3)))==cyc(6,(0,3,1)), 'edge completion identity')
    require(compose(cyc(6,(0,1,2)),cyc(6,(0,2,3)))==cyc(6,(1,2,3)), 'three-cycle factorization identity')
    controls=[]
    for b,n in [(4,1),(4,3),(5,2),(5,3),(6,5),(7,2),(8,3),(9,2)]:
        validate_domain(b,n,four=b==4); seed=family_seed(quartet,b); word=family_word(beta,b,n)
        require(all(cycle_lengths(p)==[2,2] for p in seed), 'control generator class')
        require(hurwitz(word,seed)==seed, 'control full braid tuple return')
        require(cycle_lengths(strand_permutation(family_word(beta,b,1),b))==[b], 'control b-cycle')
        require(len(cycle_lengths(strand_permutation(word,b),include_fixed=True))==gcd(b,n), 'component formula control')
        for j in range(4,b): require(cycle_lengths(compose(seed[j-1],seed[j]))==[3,3], 'consecutive extra-generator product')
        expected_runs=[4*n,6*n+1]+[2*n]*(b-4)
        require(seesaws(word,b)==expected_runs and min(expected_runs)>=4, 'exact linear seesaws')
        require(len(word)==n*(3*b+5), 'crossing count control')
        genus_numerator=len(word)-b+1; require(genus_numerator%2==0, 'integral genus control')
        controls.append({'b':b,'n':n,'word_length':len(word),'linear_seesaws':expected_runs,
                         'genus_from_written_homogeneous_theorem':genus_numerator//2,
                         'meridional_rank_from_written_proof':b})
    literal=verify_literal(data['literal_b5_n2'],beta,quartet)
    cox=symmetric_control(data['finite_symmetric_control'],beta)
    bad_pd=deepcopy(data['literal_b5_n2']); bad_pd['pd'][0][0]=bad_pd['pd'][0][1]
    bad_color=deepcopy(data['literal_b5_n2']); old=bad_color['arc_permutations'][0]
    bad_color['arc_permutations'][0]=next(p for p in bad_color['generator_permutations'] if p!=old)
    bad_word=data['base_short_words']['ordinary_a']['word_in_quartet'][:-1]
    bad_seed=family_seed(quartet,5); bad_seed[4]=cyc(10,(1,8),(2,9))
    negatives=[
        expect_rejection('corrupted short base generation word',lambda: require(evaluate(bad_word,quartet)==cyc(8,(0,1,2,3,4,5,6)), 'bad base word')),
        expect_rejection('corrupted literal PD edge multiplicity',lambda: verify_literal(bad_pd,beta,quartet)),
        expect_rejection('single corrupted literal arc permutation',lambda: verify_literal(bad_color,beta,quartet)),
        expect_rejection('extension generator violating cubic return',lambda: require(hurwitz([4,4,4],bad_seed)==bad_seed,'bad extension return')),
        expect_rejection('extra unmatched signed base crossing',lambda: require(hurwitz(beta+[1],quartet)==quartet,'bad braid return')),
        expect_rejection('n=1 outside b>=5 prime-family domain',lambda: validate_domain(5,1)),
        expect_rejection('noncoprime power is not a knot',lambda: validate_domain(5,5)),
        expect_rejection('even four-braid power is not a knot',lambda: validate_domain(4,2,four=True))]
    sources={}
    for name,entry in data['sources'].items():
        p=HERE/entry['path']; digest=hashlib.sha256(p.read_bytes()).hexdigest()
        require(digest==entry['sha256'], 'primary source hash: '+name)
        sources[name]={**entry,'bytes':p.stat().st_size}
    return {'base_short_word_multiplications_checked':19,'base_braid_letters_checked':17,
            'base_A8_generation':'two explicit words yield the ordinary generators in the written lemma',
            'induction_identity_controls':induction,'finite_parameter_controls':controls,
            'uniform_assertions_justification':'written group induction, braid insertion, seesaw formulas, and imported source theorems; finite controls are supporting checks',
            'literal_b5_n2':literal,'finite_S6_control':cox,'negative_controls':negatives,'primary_sources':sources}

def topology_replay(data):
    import importlib.metadata
    from spherogram import Link
    literal=data['literal_b5_n2']; link=Link(braid_closure=literal['literal_braid_word'])
    require([list(row) for row in link.PD_code()]==literal['pd'], 'fresh topology literal PD')
    require(json.loads(json.dumps(link.DT_code()))==literal['DT'], 'fresh topology DT')
    require(len(link.link_components)==1, 'fresh topology one component')
    return {'status':'PASS','backend':'Spherogram','version':importlib.metadata.version('spherogram'),
            'fresh_literal_PD_and_DT_match':True,'components':1,'named_census_lookup_used':False}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--topology',action='store_true',help='also reconstruct the literal braid PD with Spherogram')
    parser.add_argument('--check-only',action='store_true',help='do not replace the canonical receipt')
    parser.add_argument('--pdf',type=Path,help='record a separately built manuscript PDF hash')
    args=parser.parse_args(); data=json.loads(INPUT.read_text()); result=run(data)
    result['schema']='alternating-target-prime-families-verification-v1'
    result['status']='PASS'; result['checked_at_utc']=datetime.now(timezone.utc).isoformat()
    result['review_boundaries']={'independent_replay':False,'formal_verification':False,'external_review':False}
    result['topology_replay']=topology_replay(data) if args.topology else {'status':'not requested','required_for_standard_library_checks':False}
    result['artifacts']={p.name:{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in (INPUT,Path(__file__),TEX)}
    if args.pdf:
        pdf=args.pdf.resolve(); result['pdf_artifact']={'path':str(pdf),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'bytes':pdf.stat().st_size}
    if not args.check_only:
        if RECEIPT.exists():
            previous=json.loads(RECEIPT.read_text())
            for key in ('internal_pdf_reviews','review_history'):
                if key in previous: result[key]=previous[key]
            if not args.pdf and previous.get('artifacts')==result['artifacts'] and 'pdf_artifact' in previous:
                result['pdf_artifact']=previous['pdf_artifact']
        text=json.dumps(result,indent=2)+'\n'; temporary=RECEIPT.with_suffix('.json.tmp'); temporary.write_text(text); temporary.replace(RECEIPT)
    print('PASS: A8 short words; exact induction controls; literal b=5,n=2 PD/coloring; eight parameter controls; 216 S6 tree orbits; eight rejecting controls.')
    if args.topology: print('PASS: fresh Spherogram literal PD and DT reconstruction.')
    if not args.check_only: print(RECEIPT)

if __name__=='__main__': main()
