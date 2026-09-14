#!/usr/bin/env sage -python
"""Exact arithmetic companion to the common-quartic branch-stratum proof.

This is one checker, not an independent replay or a full-chart elimination.
Old point-count and branch-comparison proofs are retained inputs.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import GF, QQ, ZZ, PolynomialRing, matrix
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-common-quartic-branch-strata-v1'
SOURCE = 'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
FROB = 'artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_v1.json'
FROB_REPLAY = 'artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_verification_v1.json'
INPUTS = [SOURCE, FROB, FROB_REPLAY,
    'artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1/independent-replay.json',
    'artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1/independent-replay-final.json',
    'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/result.json',
    'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/independent-replay.json',
    'elkies-k3/Q80_ALL_SMOOTH_GENUS_ONE_BISECTIONS_2026-09-13.md',
    'elkies-k3/R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md',
    'elliptic-curves/rank-jump/THE_LAST_GLOBAL_CLASS_IS_ABSENT.md',
    'elliptic-curves/notes/DET1092_SURFACE_BRAUER_TRIVIALITY_2026-09-09.md',
    'elliptic-curves/notes/DET1092_UNRAMIFIED_KUMMER_OBSTRUCTION_2026-09-09.md']
NAMES = ['published-r17', 'alternate-q80', 'curve302-parent', 'x1092-class1']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write_new(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')


def freeze(out):
    write_new(out / 'input.json', {
        'schema': 'common-quartic-branch-strata-input-v1',
        'bindings': {p: digest(ROOT / p) for p in INPUTS},
        'checker_sha256': digest(Path(__file__)),
        'parents': NAMES, 'good_patch_values': list(range(1,10)), 'good_prime': 1009,
        'cpu_seconds': 40, 'address_space_bytes': 4*1024**3,
        'sage_version': version,
        'scope': 'Universal factor identities and degrees; arithmetic hypotheses for a written Q80 same-branch genus-one exclusion. No complete coefficient elimination or point census.',
    })
    print('Frozen branch-stratum proof inputs', flush=True)


def universal_checks():
    R = PolynomialRing(QQ, names='g,d,u,v,a,b,z,c,l')
    g,d,u,v,a,b,z,c,l = R.gens()
    h, D = c*g*u*v, c*l*g*d
    r1, r2 = (u*a+v*b)/2, (v*b-u*a)/2
    A = l*d*a*b - 3*z*z - 3*z*h - h*h
    B = D*r2*r2 + 2*z**3 + 3*z*z*h + z*h*h-l*z*d*a*b
    x1, x2 = z+h, z
    f = lambda x: x**3+A*x+B
    assert f(x1) == D*r1*r1 and f(x2) == D*r2*r2
    assert h*(x1*x1+x1*x2+x2*x2+A) == D*(r1*r1-r2*r2)
    assert f(x1)*f(x2) == (D*r1*r2)**2
    W = h*(2*x1+x2)*(x1+2*x2)
    assert (4*A**3+27*B**2+W**2).subs({d:0}) == 0
    # The parity triple at a branch point is the complement of its chosen root.
    patterns = [[int(i != j) for i in range(3)] for j in range(3)]
    for i in range(3):
        for j in range(3):
            assert all((x+y) % 2 == 0 for x,y in zip(patterns[i],patterns[j])) == (i == j)
    strata = []
    for k in range(5):
        for j in range(5-k):
            degrees = [k,4-k,j,4-k-j,4-j,k+j,4,0,0]
            weighted = lambda p: max(sum(e*w for e,w in zip(ex, degrees)) for ex in p.exponents())
            assert weighted(D) == 4 and weighted(h) == 4
            assert weighted(r1) == weighted(r2) == 4
            assert weighted(A) <= 8 and weighted(B) <= 12
            I = k+2*j
            cross = 4-I
            hm, hp = 8+2*I, 24-2*I
            assert hm+hp == 32 and 64-cross**2 >= 48
            # 21 parameters, including the one-dimensional a,b rescaling.
            assert k+(4-k)+j+(4-k-j)+(5-j)+(k+j+1)+5+2 == 21
            strata.append({'k':k, 'j':j, 'degrees':degrees, 'intersection':I,
                           'height_matrix':[[8,cross],[cross,8]],
                           'difference_height':hm, 'sum_height':hp})
    assert len(strata) == 15
    last = strata[-1]
    assert (last['k'],last['j']) == (4,0)
    assert last['sum_height'] == last['difference_height'] == 16
    assert (16-16)//2 == 0 and 1+(16-16)//4 == 1
    # Degree-six Galois cover with 24 transposition inertia elements.
    assert (-2*6+24*3+2)//2 == 31
    return strata


def patch_checks(packet):
    S = PolynomialRing(GF(packet['good_prime']), 't'); t=S.gen()
    parents = json.loads((ROOT/SOURCE).read_text())['parents']
    assert [r['name'] for r in parents] == packet['parents'] == NAMES
    places = packet['good_patch_values']
    assert places == list(range(1,10)) and len(set(places)) > 4+4
    out = []
    for row in parents:
        A,B = [S([QQ(v) for v in row[key]]) for key in ('A','B')]
        Delta = 4*A**3+27*B**2
        assert Delta.degree() == 24 and Delta.gcd(Delta.derivative()) == 1
        residues = [int(Delta(s)) for s in places]
        assert all(residues)
        # Coefficient-level weighted transformation includes lower-degree terms.
        for s in places:
            transform = lambda f,n: sum(f[i]*(s*t+1)**i*t**(n-i) for i in range(n+1))
            assert transform(Delta,24) == 4*transform(A,8)**3+27*transform(B,12)**2
            assert transform(Delta,24)[24] == Delta(s)
        out.append({'name':row['name'], 'prime':packet['good_prime'],
                    'good_patch_values':places, 'discriminant_residues':residues})
    return out


def retained_arithmetic():
    old = json.loads((ROOT/FROB).read_text())
    replay = json.loads((ROOT/FROB_REPLAY).read_text())
    assert replay['bindings'][FROB] == digest(ROOT/FROB)
    # Check the immediate producer and independent-replay bindings as frozen.
    for data in (old,replay):
        for p, expected in data['bindings'].items():
            assert digest(ROOT/p) == expected, p
        assert data['status'] == 'PASS'
    assert old['prime'] == 131 and old['finite_field_generic_rank'] == 17
    assert abs(matrix(ZZ,old['gram']).det()) == old['gram_determinant'] == 948
    assert old['kummer_character_rank'] == replay['character_rank'] == 17
    assert old['two_saturation'] == 'PROVED'
    assert old['normalized_Artin_Tate_factor'] == replay['Artin_Tate_factor'] == '948/131'
    p = ZZ(131)
    T1, T2 = [ZZ(x)-17*p**i for i,x in enumerate(old['Frobenius_traces'],1)]
    assert [T1,T2] == [-343,27783]
    assert -p-212 == T1 and p*p+212**2-2*p*p == T2
    assert (4*p+2*212)/948 == 1
    return {'prime':131, 'arithmetic_picard_rank':19,
            'arithmetic_picard_determinant':948, 'finite_field_brauer_order':1,
            'new_global_deduction':'Br(X)[2] consists of constants on the determinant948 surface',
            'counts_and_saturation':'inherited source-bound arithmetic proofs'}


def old_control_patches():
    data=json.loads((ROOT/'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/result.json').read_text())
    R=PolynomialRing(QQ,'t'); t=R.gen()
    out=[]
    for row in data['controls']:
        # Both old controls have h=-1. Move a good finite point to infinity.
        s=1
        transform=lambda f,n: sum(f[i]*(s*t+1)**i*t**(n-i) for i in range(n+1))
        D0=R(list(map(QQ,row['D'])))
        s0=R(list(map(QQ,row['s'])))
        D=transform(D0,4)
        x1,x2=R(0),t**4
        r1,r2=transform(s0,4),transform(s0+1,4)
        h=x1-x2; c=h.leading_coefficient()
        g=D.gcd(h).monic(); d=R(D/g).monic()
        h0=R(h/c/g); u=h0.gcd(r1-r2).monic(); v=R(h0/u)
        a=R((r1-r2)/u); b=R((r1+r2)/v)
        lam=D.leading_coefficient()/c
        z=x2
        A=lam*d*a*b-3*z*z-3*z*h-h*h
        B=D*(v*b-u*a)**2/4+2*z**3+3*z*z*h+z*h*h-lam*z*d*a*b
        assert A == transform(R(list(map(QQ,row['A']))),8)
        assert B == transform(R(list(map(QQ,row['B']))),12)
        assert h.degree() == D.degree() == 4
        assert (g.degree(),u.degree()) == (0,4)
        assert R(v).gcd(R(a)) == 1 and d.gcd(R(u*v)) == 1
        # The scalar, rather than only proportionality, must be preserved.
        assert D == c*lam*g*d
        out.append({'name':row['name'], 'patch':s, 'k':0, 'j':4,
                    'genus':row['normalization_genus'], 'literal_cover_scalar':str(c*lam)})
    assert len(out)==2
    return out


def check(out, record):
    packet=json.loads((out/'input.json').read_text())
    assert packet['checker_sha256'] == digest(Path(__file__))
    assert set(packet['bindings']) == set(INPUTS)
    for p,expected in packet['bindings'].items():
        assert digest(ROOT/p) == expected, p
    resource.setrlimit(resource.RLIMIT_CPU,(40,45))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time()
    result={'schema':'common-quartic-branch-strata-result-v1','status':'PASS',
            'input_sha256':digest(out/'input.json'),
            'strata':universal_checks(), 'patches':patch_checks(packet),
            'arithmetic_input':retained_arithmetic(), 'control_patches':old_control_patches(),
            'q80_k4_genus_one':'EXCLUDED by the written descent and smooth-image proof',
            'fixed_parent_full_family':'UNKNOWN', 'positive_mw17_endpoint':False,
            'independent_replay':False, 'formal_verification':False}
    result['cpu_seconds']=round(time.process_time()-start,6)
    if record:
        write_new(out/'result.json',result)
    else:
        saved=json.loads((out/'result.json').read_text())
        assert {k:v for k,v in saved.items() if k!='cpu_seconds'} == {k:v for k,v in result.items() if k!='cpu_seconds'}
    print(json.dumps({'status':'PASS','factor_strata':15,'good_patches_per_parent':9,
                      'q80_exclusion':'k=4, genus1 only','full_target':'OPEN',
                      'cpu_seconds':result['cpu_seconds']}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze',action='store_true')
    parser.add_argument('--record',action='store_true')
    parser.add_argument('--output',type=Path,default=DEFAULT)
    args=parser.parse_args()
    if args.freeze:
        assert not args.record
        freeze(args.output)
    else:
        check(args.output,args.record)


if __name__=='__main__':
    main()
