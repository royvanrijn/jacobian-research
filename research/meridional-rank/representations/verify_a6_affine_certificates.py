#!/usr/bin/env python3
"""Exact positive A6 affine knot certificates; standard library by default.

--cover-homology adds the independent full directed-edge six-sheet control.
--census rereads the named records from the pinned optional census packages.
Neither control is a premise of the literal matrix lower-bound proof.
"""
import argparse
from collections import defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'a6_affine_inputs.json'
NAMES = ['13a_' + str(i) for i in [304,362,623,629,3232,3259,3374,3375]] + [
    '13n_' + str(i) for i in [834,838,839,991,992,993,1367,1370,1371,1374,
                            1375,1376,3814,3840,3911,3914]]
SOURCE_PINS = {'sources/wirtinger-systems-published-2020.pdf':
               '884b32bd4c6f9a8fe250e6b2d6342cbf0dadcb3d8e3a416c76cb7efe7e5e40b0'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rank(rows, columns):
    a = [[v % 3 for v in row] for row in rows]
    assert all(len(row) == columns for row in a)
    r = 0
    for j in range(columns):
        i = next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None:
            continue
        a[r],a[i] = a[i],a[r]
        factor = pow(a[r][j],-1,3)
        a[r] = [factor*v % 3 for v in a[r]]
        for i in range(r+1,len(a)):
            if a[i][j]:
                factor = a[i][j]
                a[i] = [(v-factor*w) % 3 for v,w in zip(a[i],a[r])]
        r += 1
        if r == len(a):
            break
    return r

def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) % 3 for col in zip(*b)] for row in a]

def difference(a):
    return [[(v-int(i==j)) % 3 for j,v in enumerate(row)] for i,row in enumerate(a)]

def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def projection_group(gens):
    identity = tuple(range(6))
    found = {identity}
    queue = [identity]
    for p in queue:
        for q in gens:
            r = compose(q,p)
            if r not in found:
                found.add(r)
                queue.append(r)
    return found

def diagram(pd):
    ends = defaultdict(list)
    over = defaultdict(set)
    for i,row in enumerate(pd):
        assert len(row) == len(set(row)) == 4
        for h,e in enumerate(row):
            ends[e].append(4*i+h)
        over[row[1]].add(row[3])
        over[row[3]].add(row[1])
    n = len(pd)
    assert len(ends) == 2*n and all(len(v)==2 for v in ends.values())
    twin = {a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
    visited = set()
    v = 0
    for _ in range(2*n):
        assert v not in visited
        visited.add(v)
        v = twin[4*(v//4)+(v+2)%4]
    assert v == 0 and len({v//4 for v in visited}) == n
    unseen = set(twin)
    faces = 0
    while unseen:
        start = v = min(unseen)
        while True:
            unseen.remove(v)
            w = twin[v]
            v = 4*(w//4)+(w+1)%4
            if v == start:
                break
        faces += 1
    assert faces == n+2, 'literal diagram is not spherical'
    unseen = set(ends)
    arcs = []
    while unseen:
        start = min(unseen)
        component = {start}
        queue = [start]
        for e in queue:
            for f in over[e]-component:
                component.add(f)
                queue.append(f)
        unseen -= component
        arcs.append(sorted(component))
    arcs.sort()
    assert len(arcs) == n
    where = {e:i for i,arc in enumerate(arcs) for e in arc}
    return [(where[b],where[a],where[c]) for a,b,c,d in pd], twin

def named_binding(case):
    _,a = diagram(case['pd'])
    _,b = diagram(case['census_pd'])
    w = case['census_map_witness']
    f = w['dart_images']
    orientation,swap = w['orientation'],w['over_under_swap']
    # These certificates preserve both plane orientation and crossing data.
    assert orientation == 1 and swap == 0
    assert sorted(f) == list(range(4*len(case['pd'])))
    for x,y in enumerate(f):
        assert f[a[x]] == b[y]
        assert x % 2 == y % 2
        assert f[4*(x//4)+(x+1)%4] == 4*(y//4)+(y+1)%4

def essential_matrix(m):
    # Quotient by c=(1,1,1,1,1,1,0,0); use e0,...,e4,e6,e7.
    columns = [0,1,2,3,4,6,7]
    return [[(m[i][j]-m[5][j]) % 3 for j in columns] for i in range(5)] + [
            [m[i][j] for j in columns] for i in [6,7]]

def check_case(case):
    assert case['name'] in NAMES
    assert len(case['pd']) == 13
    named_binding(case)
    crossings,_ = diagram(case['pd'])
    seeds = case['seed_arcs']
    assert len(seeds) == len(set(seeds)) == 4
    assert all(0 <= i < 13 for i in seeds)
    known = set(seeds)
    for o,a,b in case['seed_trace']:
        assert o in known and a in known and b not in known
        assert (o,a,b) in crossings or (o,b,a) in crossings
        known.add(b)
    assert known == set(range(13))
    perms = case['permutations']
    columns = case['cocycle_columns_F3']
    matrices = case['block_matrices_F3']
    assert len(perms) == len(columns) == len(matrices) == 13
    rebuilt = []
    for p,z,m in zip(perms,columns,matrices):
        assert sorted(p) == list(range(6))
        assert all(p[p[i]] == i for i in range(6))
        assert sum(p[i] != i for i in range(6)) == 4
        assert len(z) == 6 and all(len(row)==2 for row in z)
        assert all(v in [0,1,2] for row in z for v in row)
        a = [[int(p[j]==i) for j in range(6)]+z[i] for i in range(6)] + [
            [0]*6+[1,0], [0]*6+[0,1]]
        assert a == m, 'literal block matrix changed'
        rebuilt.append(a)
    identity = [[int(i==j) for j in range(8)] for i in range(8)]
    assert all(multiply(m,m)==identity for m in matrices)
    assert all(rank(difference(m),8)==2 for m in matrices)
    for o,a,b in crossings:
        assert multiply(multiply(matrices[o],matrices[a]),matrices[o]) == matrices[b], 'Wirtinger relation fails'
    stacked_rank = rank([r for a in seeds for r in difference(matrices[a])],8)
    assert stacked_rank == 7, 'common fixed space is larger than constant line'
    common = [1]*6+[0,0]
    assert all(sum(x*y for x,y in zip(row,common)) % 3 == 0
               for m in matrices for row in difference(m))
    group = projection_group([tuple(perms[a]) for a in seeds])
    assert len(group) == 360
    # Double transpositions are even, so this order identifies the full A6.
    assert all(sum(p[i]>p[j] for i in range(6) for j in range(i+1,6)) % 2 == 0 for p in group)
    essential = [essential_matrix(m) for m in matrices]
    assert all(rank(difference(m),7)==2 for m in essential)
    assert rank([r for a in seeds for r in difference(essential[a])],7)==7
    return {'name':case['name'],'projection_order':360,'ambient_dimension':8,
            'ambient_common_fixed_dimension':1,'essential_dimension':7,
            'essential_common_fixed_dimension':0,'meridian_displacement_rank':2,
            'meridian_class_generating_rank':4,'meridional_lower':4,
            'bridge_upper':4,'mr_equals_bridge':4}

def cover_homology(case):
    """All13 lifted generators, crossing relators and all cycle fillings.

    A backward edge differs from a forward edge by its meridian-cycle
    filling. Thus swapping either crossing sign gives the same filled row
    span, even at fixed sheets, and the uniform relator convention suffices.
    """
    perms = case['permutations']
    crossings,_ = diagram(case['pd'])
    n = 6
    edges = [(g,v) for g in range(13) for v in range(n)]
    parent = list(range(n))
    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v
    tree = set()
    for j,(g,v) in enumerate(edges):
        w = perms[g][v]
        if root(v) != root(w):
            parent[root(v)] = root(w)
            tree.add(j)
    assert len(tree) == 5
    chords = [j for j in range(78) if j not in tree]
    rows = []
    for o,a,b in crossings:
        for start in range(n):
            v = start
            row = [0]*78
            # Inverse edges are traversed backwards at the new endpoint.
            for g,sign in [(o,1),(a,1),(o,-1),(b,-1)]:
                if sign == 1:
                    row[g*n+v] += 1
                    v = perms[g][v]
                else:
                    v = perms[g][v]
                    row[g*n+v] -= 1
            assert v == start
            rows.append([row[j] for j in chords])
    for g,p in enumerate(perms):
        unseen = set(range(n))
        while unseen:
            start = v = min(unseen)
            row = [0]*78
            while True:
                unseen.remove(v)
                row[g*n+v] += 1
                v = p[v]
                if v == start:
                    break
            rows.append([row[j] for j in chords])
    assert len(rows) == 78+52 and len(chords) == 73
    r = rank(rows,73)
    assert r == 71
    return {'columns':73,'lifted_relator_rows':78,'meridian_cycle_fill_rows':52,
            'rank_F3':71,'H1_F3_dimension':2,'tree_edges':sorted(tree)}

def census_reread(data):
    import importlib.metadata
    from spherogram import Link
    import snappy_15_knots
    provenance = data['census_provenance']
    assert importlib.metadata.version('snappy_15_knots') == provenance['version']
    assert importlib.metadata.version('spherogram') == provenance['spherogram_version']
    package = Path(snappy_15_knots.__file__).parent
    for pin in provenance['source_files']:
        source = package / pin['path']
        assert source.stat().st_size == pin['bytes'] and digest(source) == pin['sha256']
    for case in data['cases']:
        link = Link('K'+case['name'].replace('_',''))
        assert [list(row) for row in link.PD_code()] == case['census_pd']
        assert json.loads(json.dumps(link.DT_code())) == case['census_DT']

def rejection(call):
    try:
        call()
    except (AssertionError,IndexError,KeyError,ValueError):
        return True
    raise AssertionError('corrupted certificate accepted')

# Exact Schreier kernel and two-generator controls.
def basis_mod3(rows):
 basis={}
 for row in rows:
  v=[int(x)%3 for x in row]
  for k in sorted(basis):
   if v[k]:v=[(a-v[k]*b)%3 for a,b in zip(v,basis[k])]
  k=next((k for k,x in enumerate(v) if x),None)
  if k is not None:
   factor=pow(v[k],-1,3);basis[k]=[factor*x%3 for x in v]
 return [basis[k] for k in sorted(basis)]

ID6=tuple(range(6))
ONE_AFFINE=(ID6,(0,)*12)

def affine_multiply(a,b):
 perm=tuple(a[0][b[0][i]] for i in range(6));v=list(a[1])
 for i in range(6):
  for j in range(2):v[2*a[0][i]+j]=(v[2*a[0][i]+j]+b[1][2*i+j])%3
 return perm,tuple(v)

def affine_inverse(a):
 perm=tuple(a[0].index(i) for i in range(6));v=[0]*12
 for i in range(6):
  for j in range(2):v[2*perm[i]+j]=-a[1][2*i+j]%3
 return perm,tuple(v)

def affine_word(word,generators):
 result=ONE_AFFINE
 for i in word:result=affine_multiply(result,generators[i])
 return result

def quotient_translation(v):
 return tuple((v[2*i+j]-v[10+j])%3 for i in range(5) for j in range(2))

def schreier_structure(generators):
 sections={ID6:ONE_AFFINE};queue=[ID6]
 for p in queue:
  for s in generators:
   v=affine_multiply(sections[p],s)
   if v[0] not in sections:sections[v[0]]=v;queue.append(v[0])
 rows=[]
 for p in queue:
  for s in generators:
   v=affine_multiply(sections[p],s)
   t=affine_multiply(v,affine_inverse(sections[v[0]]))
   assert t[0]==ID6
   rows.append(t[1])
 basis=basis_mod3(rows)
 quotient_basis=basis_mod3([quotient_translation(v) for v in basis])
 return {'projection_order':len(queue),'kernel_dimension_8D':len(basis),'image_order_8D':len(queue)*3**len(basis),'kernel_dimension_7D':len(quotient_basis),'image_order_7D':len(queue)*3**len(quotient_basis),'schreier_loops_checked':len(rows)},basis,quotient_basis

def quotient_matrix(perm,z):
 # Quotient coordinates are (w0-w5,...,w4-w5,s0,s1).
 # The section takes w5=0; this formula is D M E, with D E=I.
 out=[[0]*7 for _ in range(7)]
 for i in range(5):
  for j in range(5):out[i][j]=(int(perm[j]==i)-int(perm[j]==5))%3
  for j in range(2):out[i][5+j]=(z[i][j]-z[5][j])%3
 out[5][5]=out[6][6]=1
 return out

def matmul(a,b):
 return [[sum(x*y for x,y in zip(row,col))%3 for col in zip(*b)] for row in a]

def group_structure(case):
 permutations=case['permutations'];columns=case['cocycle_columns_F3'];matrices=case['block_matrices_F3']
 arc_elements=[(tuple(p),tuple(int(v)%3 for row in z for v in row)) for p,z in zip(permutations,columns)]
 generators=[arc_elements[i] for i in case['seed_arcs']]
 full,basis,qbasis=schreier_structure(generators)
 assert full['projection_order']==360
 assert full['kernel_dimension_8D']==6 and full['kernel_dimension_7D']==4
 pair_words=([[0,1,2],[0,2,3]] if case['name'] in ('13a_629','13a_3375') else [[0,2],[1,3]])
 assert pair_words == case['ordinary_pair_words_on_seed_indices']
 pair,_,_=schreier_structure([affine_word(w,generators) for w in pair_words])
 for key in ('projection_order','kernel_dimension_8D','image_order_8D','kernel_dimension_7D','image_order_7D'):
  assert pair[key]==full[key],(case['name'],key)
 qmatrices=[quotient_matrix(p,z) for p,z in zip(permutations,columns)]
 D=[[int(i==j)-int(j==5) if i<5 else int(j==i+1) for j in range(8)] for i in range(7)]
 D=[[x%3 for x in row] for row in D]
 I7=[[int(i==j) for j in range(7)] for i in range(7)]
 displacement=[]
 for m,q in zip(matrices,qmatrices):
  assert matmul(D,m)==matmul(q,D)
  assert matmul(q,q)==I7
  delta=[[(a-b)%3 for a,b in zip(r,s)] for r,s in zip(q,I7)]
  assert len(basis_mod3(delta))==2
  displacement.extend(delta)
 assert len(basis_mod3(displacement))==7
 return {**full,'ordinary_pair_words_zero_based_seed_positions':pair_words,'ordinary_generating_rank_8D':2,'ordinary_generating_rank_7D':2,'meridian_class_generating_rank_8D':4,'meridian_class_generating_rank_7D':4,'kernel_basis_8D':basis,'kernel_basis_7D':qbasis,'quotient_common_fixed_dimension':0,'quotient_meridian_displacement_rank':2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cover-homology',action='store_true')
    parser.add_argument('--census',action='store_true')
    parser.add_argument('--output',type=Path,default=HERE/'a6_affine_verification.json')
    args = parser.parse_args()
    data = json.loads(INPUT.read_text())
    assert data['field'] == 3 and data['ambient_dimension'] == 8
    assert [c['name'] for c in data['cases']] == NAMES
    assert data['source_sha256'] == SOURCE_PINS
    for name,expected in SOURCE_PINS.items():
        assert digest(HERE/name) == expected
    results = [check_case(c) for c in data['cases']]
    for case,result in zip(data['cases'],results):
        result['finite_image_structure'] = group_structure(case)
    if args.cover_homology:
        for case,result in zip(data['cases'],results):
            result['independent_cover_homology'] = cover_homology(case)
    if args.census:
        census_reread(data)
    bad = deepcopy(data['cases'][0])
    bad['block_matrices_F3'][0][0][0] = 2
    controls = {'changed_matrix':rejection(lambda:check_case(bad))}
    bad = deepcopy(data['cases'][0])
    bad['seed_trace'][0][2] = bad['seed_arcs'][0]
    controls['invalid_seed_move'] = rejection(lambda:check_case(bad))
    bad = deepcopy(data['cases'][0])
    bad['census_map_witness']['dart_images'][0:2] = [1,0]
    controls['invalid_census_map'] = rejection(lambda:check_case(bad))
    bad = deepcopy(data['cases'][0])
    for z,m in zip(bad['cocycle_columns_F3'],bad['block_matrices_F3']):
        for i in range(6):
            z[i] = [0,0]
            m[i][6:] = [0,0]
    controls['principal_zero_lift'] = rejection(lambda:check_case(bad))
    report = {'status':'PASS','python_version':platform.python_version(),
              'named_cases':24,'results':results,'negative_controls':controls,
              'optional_cover_homology_replayed':args.cover_homology,
              'optional_named_census_reread':args.census,
              'input_sha256':{'a6_affine_inputs.json':digest(INPUT),
                              Path(__file__).name:digest(Path(__file__)),**SOURCE_PINS},
              'scope':'Direct literal matrix, exact four-seed, planar named-map and meridian-class rank proofs. Optional controls are reported separately; no complete negative search or worldwide priority claim.',
              'assurance':'Internal exact replay, with published bridge=Wirtinger theorem and named-census provenance imported. No formal verification or external specialist review.'}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','named_cases':24,'mr_equals_bridge':4,
                      'cover_homology_replayed':args.cover_homology,
                      'named_census_reread':args.census,'image_orders':[262440,29160],
                      'ordinary_generating_ranks':[2,2],'negative_controls':controls},indent=2))

if __name__ == '__main__':
    main()
