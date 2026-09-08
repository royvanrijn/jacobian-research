#!/usr/bin/env sage-python
"""Generic-only fixed-place compatibility codes; first seed is evaluation only.

Nine unchanged addresses and the equation-defined bad primes of302.
One case per call, 25 seconds, no new rational points or full Selmer group.
"""
import argparse
import hashlib
import json
import signal
import sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,pari,prod

CAS = Path(__file__).resolve().parent
sys.path.insert(0,str(CAS))
from research_runtime.local_kummer import LocalSquareclasses
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_seed_local_code_v3'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(p,d):
    text = json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)


def signs_at_roots(f,x):
    seq = [f,f.derivative()]
    while seq[-1].degree() > 0:
        seq.append(-(seq[-2] % seq[-1]))
    def variations(signs):
        signs = [s for s in signs if s]
        return sum(a != b for a,b in zip(signs,signs[1:]))
    minus = variations([int(g.leading_coefficient().sign())*(-1)**g.degree() for g in seq])
    plus = variations([int(g.leading_coefficient().sign()) for g in seq])
    at = variations([int(g(x).sign()) for g in seq])
    n,left = minus-plus,minus-at
    assert n in [1,3] and f(x)
    return [0]*left+[1]*(n-left)


def build(index):
    paths = [ART/'curve302_recovered_mw17_parent_v1.json',
             ART/'rank_jump_curve302_strict_constructor_arithmetic_inputs_v1.json',
             ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent,arith,roster = [json.loads(p.read_text()) for p in paths]
    primes = [int(p) for p,e in arith['discriminant_factors']]
    protocol = {'classification':'generic-only fixed-place compatibility diagnostic',
        'places':primes+['infinity'],
        'cases':[{'index':i,'parameter':r['parameter']} for i,r in enumerate(roster['cases'])],
        'rule':'Compute each generic local image and its fixed character coordinates. Freeze the complete linear relation space on their direct product before reading the historical seed.',
        'limits':{'wall_seconds_per_case':25,'cases':9,'finite_places':20,
                  'new_addresses':0,'point_searches':0,'full_Selmer_runs':0,
                  'class_groups':0,'global_unit_groups':0,'V3_inputs':0,'pilot_changes':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__),CAS/'research_runtime/local_kummer.py']}}
    retain(DIR/'protocol.json',protocol)
    case = protocol['cases'][index]
    tau = QQ(case['parameter'])
    T = PolynomialRing(QQ,'t')
    def ev(r):
        return T(r['numerator'])(tau)/T(r['denominator'])(tau)
    E = EllipticCurve(QQ,[ev(r) for r in parent['a_invariants']])
    literal = [E([ev(r) for r in P]) for P in parent['basis_weierstrass_coordinates']]
    R = PolynomialRing(QQ,'x')
    x = R.gen()
    den = tau.denominator()
    f = x**3+5*den**4*x*x+(16*E.a4()+8)*den**8*x+(64*E.a6()+16)*den**12
    assert all(a in ZZ for a in f) and f.discriminant()
    coordinates = [(den**4*4*P[0],den**6*(8*P[1]+4*P[0]+4)) for P in literal]
    assert all(f(X)==Y*Y and Y for X,Y in coordinates)
    input_data = {'classification':'generic-only arithmetic input', 'case':case,
        'polynomial':list(map(str,f.list())),
        'generic_cubic_points':[[str(X),str(Y)] for X,Y in coordinates],
        'places':primes,'scale_denominator':str(den),
        'limits':protocol['limits'], 'inputs':protocol['inputs']}
    retain(DIR/f'case-{index:02d}-input.json',input_data)
    pari.allocatemem(64000000,268435456,silent=True)
    nf = pari.nfinit([pari(f),primes])  # maximal only at the declared primes
    basis = [R(a) for a in nf.nf_get_zk()]
    theta = pari(x).Mod(pari(f))
    betas = [pari(X)-theta for X,Y in coordinates]
    local = []
    engines = []
    matrices = []

    def encode_element(v):
        return list(map(str,R(v.lift()).list()))

    def square_test_input(value,engine,expected):
        vals = [int(pari.idealval(nf,value,P)) for P in engine.primes]
        factors = pari.matrix(len(vals),2,
            [entry for P,v in zip(engine.primes,vals) for entry in [P,v//2]])
        scale = pari.nfbasistoalg(nf,pari.idealappr(nf,factors))
        reduced = value/(scale*scale)
        assert [int(pari.idealval(nf,reduced,P)) for P in engine.primes] == [v%2 for v in vals]
        return {'scale':encode_element(scale),'expected_square':bool(expected)}

    for p in primes:
        engine = LocalSquareclasses(nf,p)
        engines.append(engine)
        signatures = [list(engine.signature(beta)) for beta in betas]
        M = matrix(GF(2),signatures).transpose()
        pivots = list(M.pivots())
        columns = M.matrix_from_columns(pivots)
        coeffs = [list(map(int,columns.solve_right(M.column(j)))) for j in range(17)]
        compressed = matrix(GF(2),coeffs).transpose()
        assert compressed.rank() == len(pivots)
        relations = []
        for j,beta in enumerate(betas):
            quotient = beta/prod(betas[k]**a for k,a in zip(pivots,coeffs[j]))
            relations.append(square_test_input(quotient,engine,True))
        nonzero = []
        for mask in range(1,1<<len(pivots)):
            value = prod(betas[k] for j,k in enumerate(pivots) if mask>>j&1)
            nonzero.append({'mask':mask,**square_test_input(value,engine,False)})
        row = {'place':p,'generic_rank':len(pivots),
            'full_local_point_dimension':int(engine.point_kummer_dimension),
            'prime_count':len(engine.primes),'basis_indices':pivots,
            'generic_coordinates':coeffs,'relation_square_tests':relations,
            'basis_nonsquare_tests':nonzero}
        local.append(row)
        matrices.append(compressed)
        retain(DIR/f'case-{index:02d}-place-{p}.json',row)
    real = [signs_at_roots(f,X) for X,Y in coordinates]
    M = matrix(GF(2),real).transpose()
    pivots = list(M.pivots())
    B = M.matrix_from_columns(pivots)
    real_coords = [list(map(int,B.solve_right(M.column(j)))) for j in range(17)]
    local.append({'place':'infinity','generic_rank':len(pivots),
        'full_local_point_dimension':int(len(real[0])==3),
        'basis_indices':pivots,'generic_coordinates':real_coords,'root_signs':real})
    matrices.append(matrix(GF(2),real_coords).transpose())
    H = matrices[0]
    for m in matrices[1:]:
        H = H.stack(m)
    checks = H.left_kernel().basis_matrix()
    generic = {'classification':'generic-only frozen compatibility code',
        'case':case,'polynomial':list(map(str,f.list())),
        'local_order_basis':[list(map(str,b.list())) for b in basis],
        'local':local,'generic_joint_rank':int(H.rank()),
        'product_of_generic_local_images_dimension':H.nrows(),
        'all_generic_local_images_full':all(r['generic_rank']==r['full_local_point_dimension'] for r in local),
        'compatibility_checks':[list(map(int,row)) for row in checks.rows()],
        'generic_matrix':[list(map(int,row)) for row in H.rows()],
        'meaning':'These are relations among finite local images, not Selmer dimensions or proof of rational points in a complementary local class.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'protocol.json',DIR/f'case-{index:02d}-input.json',Path(__file__)]}}
    retain(DIR/f'case-{index:02d}-generic.json',generic)
    # Evaluation boundary: only after the generic code is immutable.
    if index == 8:
        source = ART/'det1092_first_centre_rr_net_replay_v1.json'
        first = json.loads(source.read_text())
        xp,yp = map(QQ,first['reconstructed_point_literal302'])
        X,Y = 4*xp,8*yp+4*xp+4
        assert f(X)==Y*Y
        beta = pari(X)-theta
        coords = []
        tests = []
        for row,engine in zip(local[:-1],engines):
            signatures = matrix(GF(2),[list(engine.signature(b)) for b in betas]).transpose().matrix_from_columns(row['basis_indices'])
            cc = list(map(int,signatures.solve_right(matrix(GF(2),list(engine.signature(beta))).transpose())[:,0]))
            quotient = beta/prod(betas[j]**c for j,c in zip(row['basis_indices'],cc))
            tests.append({'place':row['place'],'coordinates':cc,**square_test_input(quotient,engine,True)})
            coords.extend(cc)
        sr = signs_at_roots(f,X)
        signatures = matrix(GF(2),[real[j] for j in local[-1]['basis_indices']]).transpose()
        cc = list(map(int,signatures.solve_right(matrix(GF(2),sr).transpose())[:,0]))
        coords.extend(cc)
        col = matrix(GF(2),coords).transpose()
        syndrome = list(map(int,(checks*col)[:,0]))
        assert any(syndrome)
        # At most32 generic compatibility checks exist; no subgroup census.
        assert checks.nrows() <= 5
        separators = []
        offsets = []
        start = 0
        for row in local:
            offsets.append((row['place'],start,start+row['generic_rank']))
            start += row['generic_rank']
        for mask in range(1,1<<checks.nrows()):
            word = sum((checks[j] for j in range(checks.nrows()) if mask>>j&1),checks[0]*0)
            if (word*col)[0]:
                support = [p for p,lo,hi in offsets if any(word[lo:hi])]
                separators.append({'mask':mask,'word':list(map(int,word)),
                                   'places':support,'place_count':len(support)})
        separators.sort(key=lambda r:(r['place_count'],r['mask']))
        result = {'classification':'retrospective evaluation of a generic-only frozen code',
            'status':'PASS_FIRST_SEED_LOCAL_COMPATIBILITY_DEFECT',
            'point_cubic':[str(X),str(Y)],'coordinates':coords,'syndrome':syndrome,
            'local_square_tests':tests,'real_signs':sr,
            'all_separating_checks':separators,
            'minimal_place_count':separators[0]['place_count'],
            'minimal_separator':separators[0],
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [source,DIR/f'case-{index:02d}-generic.json']}}
        retain(DIR/'first-seed-evaluation.json',result)
    print('PASS_CONSTRUCTED_GENERIC_LOCAL_CODE',index,H.nrows(),H.rank(),checks.nrows(),generic['all_generic_local_images_full'],flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    DIR.mkdir(parents=True,exist_ok=True)
    parser = argparse.ArgumentParser()
    parser.add_argument('--case',type=int,choices=range(9),required=True)
    build(parser.parse_args().case)
