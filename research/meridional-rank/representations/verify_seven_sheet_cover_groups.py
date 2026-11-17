#!/usr/bin/env python3
"""Replay exact group simplifications for the retained58 seven-sheet covers.

This checker builds an involution-folded group presentation. The separate
homology checker uses full directed edges. No historical GL3 enumeration
or symbolic group simplifier is imported by this verifier.
"""
import hashlib
import json
from copy import deepcopy
from pathlib import Path

HERE = Path(__file__).resolve().parent


def inverse_word(word):
    return [-x for x in reversed(word)]


def free_reduce(word):
    result = []
    for x in word:
        if result and result[-1] == -x:
            result.pop()
        else:
            result.append(x)
    return result


def normal_relation(word):
    word = free_reduce(word)
    while len(word)>1 and word[0] == -word[-1]:
        word = word[1:-1]
    if not word:
        return ()
    # Cyclic conjugation and inversion preserve the normal closure.
    return min(tuple(v[i:]+v[:i]) for v in (word,inverse_word(word))
               for i in range(len(word)))


def normalize_relations(words):
    return sorted(set(map(normal_relation,words)) - {()})


def folded_cover(case, action):
    pd = case['pd']
    labels = sorted(x for row in pd for x in row)
    assert labels == sorted(list(range(1,27))*2)
    components = [{v} for v in range(1,27)]
    for row in pd:
        a = next(c for c in components if row[1] in c)
        b = next(c for c in components if row[3] in c)
        if a is not b:
            components.remove(a)
            components.remove(b)
            components.append(a|b)
    components = sorted(sorted(c) for c in components)
    assert len(components) == 13
    arc = {v:i for i,c in enumerate(components) for v in c}
    images = []
    for matrix in action['meridian_matrices_F2']:
        image = tuple(sum(((r&v).bit_count()%2)<<j for j,r in
                          enumerate(matrix))-1 for v in range(1,8))
        assert sorted(image) == list(range(7))
        assert all(image[image[v]] == v for v in range(7))
        assert sum(image[v]!=v for v in range(7)) == 4
        images.append(image)
    assert len(images) == 13
    # In the filled cover, each fixed-sheet edge is identity and the two
    # directed edges of a moving pair are inverses. Keep one pair edge.
    edges, edge = [], {}
    for i,image in enumerate(images):
        for v in range(7):
            if v < image[v]:
                edges.append((i,v,image[v]))
                edge[i,v] = len(edges)
                edge[i,image[v]] = -len(edges)
    assert len(edges) == 26
    connected, tree = [{v} for v in range(7)], []
    for e,(_,v,w) in enumerate(edges,1):
        a = next(c for c in connected if v in c)
        b = next(c for c in connected if w in c)
        if a is not b:
            connected.remove(a)
            connected.remove(b)
            connected.append(a|b)
            tree.append(e)
    assert len(tree) == 6 and connected == [set(range(7))]
    words = []
    for row in pd:
        over,under1,under2 = arc[row[1]],arc[row[0]],arc[row[2]]
        for start in range(7):
            v,word = start,[]
            # Involutory meridians let either crossing sign use o a o b.
            for a in (over,under1,over,under2):
                e = edge.get((a,v))
                if e is not None and abs(e) not in tree:
                    word.append(e)
                v = images[a][v]
            assert v == start, 'crossing permutation relation fails'
            words.append(word)
    generators = set(range(1,27)) - set(tree)
    assert len(generators) == 20
    return normalize_relations(words),generators,tree


def replay(case,action,certificate):
    words,generators,tree = folded_cover(case,action)
    assert tree == certificate['tree_edges']
    for step in certificate['eliminations']:
        x = step['generator']
        defining = tuple(step['defining_relator'])
        assert x in generators, 'generator already absent'
        assert defining in words, 'defining relation absent'
        occurrences = [i for i,a in enumerate(defining) if abs(a)==x]
        assert len(occurrences)==1, 'generator is not single-occurrence'
        position = occurrences[0]
        tail = list(defining[position+1:]+defining[:position])
        replacement = inverse_word(tail) if defining[position]>0 else tail
        assert all(abs(a)!=x for a in replacement)
        words.remove(defining)
        new_words = []
        for word in words:
            substituted = []
            for a in word:
                if a == x:
                    substituted.extend(replacement)
                elif a == -x:
                    substituted.extend(inverse_word(replacement))
                else:
                    substituted.append(a)
            new_words.append(substituted)
        words = normalize_relations(new_words)
        generators.remove(x)
        assert all(abs(a) in generators for word in words for a in word)
    assert sorted(generators) == certificate['final_generators'], 'final generators differ'
    assert [list(w) for w in words] == certificate['final_relators'], 'final relations differ'
    expected = action['expected_H1_F3_dimension']
    if expected:
        assert len(generators)==1 and len(words)==1
        generator = next(iter(generators))
        assert words == [(-generator,)*3]
        assert len(certificate['eliminations']) == 19
        group = 'C3'
    else:
        assert not generators and not words
        assert len(certificate['eliminations']) == 20
        group = 'trivial'
    return {'name':case['name'],'module':action['module'],
            'cover_group':group,'group_rank':int(group=='C3'),
            'eliminations':len(certificate['eliminations'])}


def main():
    inputs = json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    certificates = json.loads((HERE/'seven_sheet_cover_group_certificates.json').read_text())
    cases = {(c['name'],a['module']):(c,a) for c in inputs['cases']
             for a in c['actions']}
    assert len(cases) == len(certificates['cases']) == 58
    results = []
    for certificate in certificates['cases']:
        case,action = cases.pop((certificate['name'],certificate['module']))
        results.append(replay(case,action,certificate))
    assert not cases
    assert sum(r['cover_group']=='C3' for r in results)==48
    assert sum(r['cover_group']=='trivial' for r in results)==10
    # A changed defining word cannot silently preserve an accepted proof.
    bad = deepcopy(certificates['cases'][0])
    bad['eliminations'][0]['defining_relator'].append(99)
    case = inputs['cases'][0]
    try:
        replay(case,case['actions'][0],bad)
    except AssertionError as error:
        assert 'defining relation absent' in str(error)
    else:
        raise AssertionError('corrupt defining word accepted')
    bad = deepcopy(certificates['cases'][0])
    bad['final_relators'][0].pop()
    try:
        replay(case,case['actions'][0],bad)
    except AssertionError as error:
        assert 'final relations differ' in str(error)
    else:
        raise AssertionError('corrupt final relation accepted')
    paths = ['verify_seven_sheet_cover_groups.py','seven_sheet_cover_inputs.json',
             'seven_sheet_cover_group_certificates.json',
             'GENERAL_BRANCHED_COVER_RANK_BOUNDS.tex']
    receipt = {'status':'PASS','actions':58,'cyclic_order3_cover_groups':48,
               'trivial_cover_groups':10,'results':results,
               'negative_controls':['corrupt defining word rejected',
                                    'corrupt final relation rejected'],
               'scope':'Exact folded branched-cover presentations and one-occurrence Tietze certificate replay on58 fixed witnesses. Planarity/source identity and four-seed bridge upper bounds belong to the separate diagram/homology checker. No exhaustive quotient search or novelty claim.',
               'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                               for p in paths}}
    (HERE/'seven_sheet_cover_group_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','actions':58,'C3':48,'trivial':10,
                      'corruptions_rejected':2},indent=2))


if __name__ == '__main__':
    main()
