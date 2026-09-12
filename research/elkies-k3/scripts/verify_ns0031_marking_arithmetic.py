#!/usr/bin/env python3
"""Replay NS0031's finite modular obstruction from a compact input projection.

Standard library only. No producer imports, searches, catalogue rebuilding or
file writes. The arithmetic K3-period implication and X0(37)(Q) classification
remain external theorem inputs. This is not an independent proof of that map.
"""

if not __debug__:
    raise RuntimeError('NS0031 replay requires assertions; remove Python -O/PYTHONOPTIMIZE')

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / 'artifacts/generated-results/elkies-k3-ns0031-replay-inputs-v1.json'
CERTIFICATE = ROOT / 'artifacts/generated-results/elkies-k3-ns0031-qq-marking-obstruction-v1.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def field(row, path):
    for key in path.split('/'):
        row = row[key]
    return row


def check_source_projection(packet, root):
    """Optional provenance audit; absent originals fail without reconstruction."""
    for source in packet['sources']:
        relative = Path(source['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path = root / relative
        assert path.is_file(), f'missing original input: {relative}'
        assert digest(path) == source['sha256'], f'changed original input: {relative}'
        rows = json.loads(path.read_text())[source['collection']]
        selected = [r for r in rows if r[source['key']] == source['value']]
        assert len(selected) == 1, f'ambiguous source selector: {relative}'
        for key, value in source['fields'].items():
            assert field(selected[0], key) == value, (relative, key)


def matmul(a, b, modulus=None):
    value = (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
             a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])
    return value if modulus is None else tuple(x % modulus for x in value)


def determinant(a):
    return a[0]*a[3]-a[1]*a[2]


def split_order(gram, certificate):
    """Reduce Clifford words directly using ei*ej+ej*ei = Gram[i,j]."""
    @lru_cache(None)
    def reduce(word):
        for i in range(len(word)-1):
            a, b = word[i:i+2]
            if a < b:
                continue
            rest = word[:i]+word[i+2:]
            result = {w: F(gram[a][a], 2)*v for w, v in reduce(rest).items()} if a == b else {
                w: F(gram[a][b])*v for w, v in reduce(rest).items()}
            if a != b:
                swapped = word[:i]+(b,a)+word[i+2:]
                for w, v in reduce(swapped).items():
                    result[w] = result.get(w, F(0))-v
            return {w:v for w,v in result.items() if v}
        return {word:F(1)}

    words = [(), (0,1), (0,2), (1,2)]
    matrices = [(F(1),F(0),F(0),F(1)), (F(0),F(1),F(0),F(0)),
                (F(4),F(0),F(0),F(0)), (F(1),F(1,4),F(148),F(0))]
    lookup = dict(zip(words, matrices))
    for a, b in product(words, repeat=2):
        expansion = reduce(a+b)
        image = tuple(sum(v*lookup[w][k] for w,v in expansion.items()) for k in range(4))
        assert image == matmul(lookup[a],lookup[b]), 'Clifford multiplication mismatch'
    recorded = certificate['even_clifford_order']
    raw = [tuple(F(x) for row in m for x in row) for m in recorded['split_basis_before_conjugation']]
    assert matrices == raw
    integral = [(m[0],4*m[1],m[2]/4,m[3]) for m in matrices]
    expected = [(1,0,0,1),(0,4,0,0),(4,0,0,0),(1,1,37,0)]
    assert integral == expected
    assert [list(m[:2])+list(m[2:]) for m in integral] == [
        [x for row in m for x in row] for m in recorded['integral_basis_after_diag_4_1_conjugation']]
    # D=x, C=37w, B=4y+w, A=x+4z+w. Solving these four equations
    # gives exactly the three recorded integral congruences; index = 4*4*37.
    assert recorded['reduced_discriminant'] == 4*4*37 == 592


def cycle_lengths(permutation):
    assert sorted(permutation) == list(range(len(permutation)))
    seen, lengths = set(), []
    for start in range(len(permutation)):
        if start in seen:
            continue
        point, size = start, 0
        while point not in seen:
            seen.add(point)
            size += 1
            point = permutation[point]
        lengths.append(size)
    return sorted(lengths)


def modular_signature(certificate):
    # Work in (Z/4)[theta]/(theta^2-theta-1); its reduction is the field F4.
    assert all((x*x-x-1) % 2 for x in range(2))
    units = [(a,b) for a,b in product(range(4),repeat=2) if (a*a+a*b-b*b) % 2]
    trace_det = sorted({((2*a+b) % 4,(a*a+a*b-b*b) % 4) for a,b in units})
    subgroup = {((a+b)%4,b,b,a) for a,b in units if (a*a+a*b-b*b) % 4 == 1}
    group = {g for g in product(range(4),repeat=4) if determinant(g) % 4 == 1}
    assert len(units)==12 and len(subgroup)==6 and len(group)==48
    cosets, labels = [], {}
    for g in sorted(group):
        if g in labels:
            continue
        coset = {matmul(h,g,4) for h in subgroup}
        assert len(coset)==6 and not (coset & labels.keys())
        labels.update({h:len(cosets) for h in coset})
        cosets.append(g)
    lines = [(0,1)]+[(1,x) for x in range(37)]
    points = [(i,l) for i in range(len(cosets)) for l in lines]
    position = {x:i for i,x in enumerate(points)}
    signatures=[]
    for operator in [(0,-1,1,0),(0,-1,1,1),(1,1,0,1)]:
        a,b,c,d=operator
        permutation=[]
        for i,(x,y) in points:
            u,v=(x*a+y*c)%37,(x*b+y*d)%37
            line=(1,v*pow(u,-1,37)%37) if u else (0,1)
            j=labels[matmul(cosets[i],operator,4)]
            permutation.append(position[(j,line)])
        signatures.append(cycle_lengths(permutation))
    e2,e3=signatures[0].count(1),signatures[1].count(1)
    cusps=signatures[2]
    index=len(points)
    genus=1+F(index,12)-F(e2,4)-F(e3,3)-F(len(cusps),2)
    recorded=certificate['norm_one_modular_curve']
    assert (index,e2,e3,cusps,genus)==(304,0,4,[4,4,148,148],23)
    for key,value in [('index_in_PSL2Z',index),('elliptic_orbits_order_2',e2),
                      ('elliptic_orbits_order_3',e3),('cusp_widths',cusps),('genus',genus),
                      ('mod_4_full_nonsplit_cartan_order',len(units)),('mod_4_norm_one_order',len(subgroup))]:
        assert recorded[key]==value,key
    assert list(map(list,trace_det))==certificate['x0_37_rational_point_gate']['nonsplit_cartan_trace_determinants_mod_4']
    assert (2,3) not in trace_det
    return set(trace_det)


def curve_invariants(ainvs):
    a1,a2,a3,a4,a6=ainvs
    b2,b4,b6=a1*a1+4*a2,2*a4+a1*a3,a3*a3+4*a6
    b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    c4=b2*b2-24*b4
    disc=-b2*b2*b8-8*b4**3-27*b6**2+9*b2*b4*b6
    return F(c4**3,disc),disc


def verify(packet, certificate):
    assert packet['schema']=='elkies-k3.ns0031-replay-inputs.v1'
    assert certificate['schema']=='elkies-k3.ns0031-qq-marking-obstruction.v1'
    assert certificate['surface_id']=='K3-d1b1381f87d69f1c'
    assert (certificate['ns_id'],certificate['ns_rank'],certificate['ns_determinant'])==('NS0031',19,1184)
    assert [(s['collection'],s['key'],s['value']) for s in packet['sources']]==[
        ('surfaces','surface_id','K3-d1b1381f87d69f1c'),
        ('surfaces','surface_id','K3-d1b1381f87d69f1c'),
        ('sources','source_id','NS0031-S001')]
    assert {s['path']:s['sha256'] for s in packet['sources']}==certificate['input_hashes']
    assert len(packet['sources'])==3
    catalogue, arithmetic, source = [s['fields'] for s in packet['sources']]
    gram=catalogue['surface_key/transcendental_gram']
    assert gram==arithmetic['literal_transcendental_gram']==certificate['transcendental_lattice']['gram']
    assert gram==[[0,0,4],[0,74,1],[4,1,-2]]
    assert 4*(0*1-74*4)==certificate['transcendental_lattice']['determinant']==-1184
    assert arithmetic['rational_isotropy/primitive_isotropic_vector']==[1,0,0]
    assert arithmetic['rational_isotropy/integral_u_split/isotropic_divisibility']==4
    assert arithmetic['clifford/quaternion_discriminant']==1
    assert arithmetic['clifford/integral_even_clifford_order/reduced_discriminant']==592
    assert source=={'determinant':1184,'source/root_type':'A1+2A7',
                    'source/root_lattice_primitive':True,'source/torsion':1,
                    'source/mw_height_gram':[['2','1'],['1','41/8']]}
    frame=certificate['source_frame']
    assert (frame['source_id'],frame['root_type'],frame['mw_rank'])==('NS0031-S001','A1+2A7',2)
    assert frame['mw_height_gram']==source['source/mw_height_gram']
    assert frame['torsion']==1 and frame['root_lattice_primitive'] is True
    split_order(gram,certificate)
    allowed=modular_signature(certificate)
    gate=certificate['x0_37_rational_point_gate']
    assert gate['noncuspidal_rational_j_values']==[-7*11**3,-7*137**3*2083**3]
    assert len(gate['excluded_lifts'])==2
    for expected, row in zip(gate['noncuspidal_rational_j_values'],gate['excluded_lifts']):
        ainvs=row['minimal_ainvariants']
        j,disc=curve_invariants(ainvs)
        assert j==expected==row['j'] and disc % 19
        a1,a2,a3,a4,a6=ainvs
        count=1+sum((y*y+a1*x*y+a3*y-x*x*x-a2*x*x-a4*x-a6)%19==0
                    for x,y in product(range(19),repeat=2))
        trace=20-count
        assert trace==row['frobenius_trace']==-6 and row['good_prime']==19
        assert [trace%4,19%4]==row['trace_determinant_mod_4']==[2,3]
        assert (trace%4,3) not in allowed and ((-trace)%4,3) not in allowed
        assert row['in_nonsplit_cartan_trace_determinant_set'] is False
    return 'PASS NS0031 finite arithmetic: Clifford order, 304 cosets, genus23 and both p19 exclusions; external moduli/classification inputs are not proved here.'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-source-projection',action='store_true',
                        help='also read and compare the three original catalogues; never rebuild them')
    args=parser.parse_args()
    packet=json.loads(PACKET.read_text())
    assert digest(CERTIFICATE)==packet['certificate_sha256'], 'primary certificate changed'
    certificate=json.loads(CERTIFICATE.read_text())
    if args.check_source_projection:
        check_source_projection(packet,ROOT)
    print(verify(packet,certificate))


if __name__=='__main__':
    main()
