#!/usr/bin/env sage-python
"""Frozen rank14 triangle dictionary, exact pencils and bounded302 inverses.

prepare: reconstruct289 words, the123 admitted edge configurations, and
         the ten rank14 pencils before new target tests.
probe:   one frozen pencil/prime, at most64 samples,300 seconds.
run:     reuse the two certified controls; test remaining rows sequentially
         at1013,1021,1009, stopping each after an exact exclusion;1200 seconds.
report:  collect exact exclusions and leave every failure/survivor UNKNOWN.
replay:  rebuild all QQ models and check every saved modular polynomial.
         Full fibre-normalization replay is probe --check per saved row.
"""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import runpy
import signal
import subprocess
import time
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, gcd, lcm, matrix, pari, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
SOURCE = ART/'elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
RANK_HELPER = ROOT/'elkies-k3/scripts/certify_curve302_triangle_mw14.sage'
FIBRE_HELPER = ROOT/'elkies-k3/scripts/probe_curve302_triangle_mw14.sage'
PROTOCOL = ART/'elkies-k3-curve302-shortword-triangles-protocol-v1.json'
OUT = ART/'elkies-k3-curve302-shortword-triangles-v1.json'
LOCAL = ROOT/'artifacts/local/elkies-k3/curve302-shortword-triangles-v1'
PUBLIC = ROOT/'elliptic-curves/cas/icarm_curve302.py'
PRIMES = [1013,1021,1009]


def digest(p): return sha256(p.read_bytes()).hexdigest()
def rows(m): return [[int(c) for c in row] for row in m.rows()]
def rat(f): return {'n':list(map(str,f.numerator().list())), 'd':list(map(str,f.denominator().list()))}
def decode(field,R,f): return R.fraction_field()(R([field(QQ(c)) for c in f['n']])/R([field(QQ(c)) for c in f['d']]))


def make_pencil(field,record):
    source = json.loads(SOURCE.read_text()); R = PolynomialRing(field,'u'); u = R.gen(); K = R.fraction_field()
    def load(f):
        return K(R([field(QQ(c)) for c in f['numerator_coefficients_low_to_high']])/R([field(QQ(c)) for c in f['denominator_coefficients_low_to_high']]))
    A = R([field(QQ(c)) for c in source['weierstrass_model']['A_coefficients_low_to_high']])
    B = R([field(QQ(c)) for c in source['weierstrass_model']['B_coefficients_low_to_high']])
    delta = 4*A**3+27*B**2
    assert delta.degree() == 24 and delta.gcd(delta.derivative()) == 1
    E = EllipticCurve(K,[A,B]); points = {}
    def point(word):
        answer = E(0)
        for i,a in enumerate(word):
            if a:
                if i not in points:
                    data = source['sections']['records'][i]
                    points[i] = E(load(data['X']),load(data['Y']))
                answer += ZZ(a)*points[i]
        return answer
    P,Q,Z = [point(record[key]) for key in ['P_word','Q_word','Z_word']]
    def meet(P):
        h = P[0].denominator(); h = h//h.gcd(h.derivative())
        assert h.degree() == 1, 'Outside finite simple-pole compiler'
        return -h[0]/h[1]
    rp,rq,lam = meet(P),meet(Q),meet(P-Q)
    assert len({rp,rq,lam}) == 3, 'Collision outside this compiler'
    assert P[0](lam) == Q[0](lam) and P[1](lam) == Q[1](lam)
    def principal(P,r):
        a = (u-r)*(-P[1]/P[0]); k = a(r)/(u-r)+a.derivative()(r)
        assert P[0].valuation(u-r) == -2 and P[1].valuation(u-r) == -3
        assert (-P[1]/P[0]-k).valuation(u-r) >= 1
        assert P[0].numerator().degree()-P[0].denominator().degree() <= 4
        assert P[1].numerator().degree()-P[1].denominator().degree() <= 6
        return k
    kp,kq = principal(P,rp),principal(Q,rq)
    c1 = (lam-rp)/((u-rp)*(u-lam)); c2 = -(lam-rq)/((u-rq)*(u-lam))
    c0 = (kp(lam)-kq(lam))/(u-lam)-c1*kp-c2*kq
    assert c1+c2
    z = c1*(Z[1]+P[1])/(Z[0]-P[0])+c2*(Z[1]+Q[1])/(Z[0]-Q[0])+c0
    assert max(z.numerator().degree(),z.denominator().degree()) == 1
    S = PolynomialRing(K,'s'); s = S.gen(); X = PolynomialRing(S,'x'); x = X.gen()
    den = (x-P[0])*(x-Q[0]); ny = c1*(x-Q[0])+c2*(x-P[0])
    n0 = c1*P[1]*(x-Q[0])+c2*Q[1]*(x-P[0])+c0*den
    f,rem = ((s*den-n0)**2-(x**3+A*x+B)*ny**2).quo_rem(den)
    assert not rem and f.degree() == 3
    return K,f,(c1,c2,c0),z,(rp,rq,lam)


def make_protocol():
    source = json.loads(SOURCE.read_text()); G = matrix(ZZ,source['sections']['height_gram'])
    N = block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G); basis = matrix.identity(ZZ,17)
    def sec(w): return vector(ZZ,[1,(w*G*w)//2]+list(w))
    words = list(basis.rows())+[basis[i]+s*basis[j] for i,j in combinations(range(17),2) for s in [-1,1]]
    six = [w for w in words if w*G*w == 6]
    zeros = sorted([w for w in words if w*G*w == 4],key=lambda w:tuple(w))
    O = sec(vector(ZZ,17)); seen = set(); roster = []; models = []
    exact_roots = runpy.run_path(str(RANK_HELPER))['exact_roots']
    for v,w in combinations(six,2):
        if abs(v*G*w) != 3: continue
        w = w if v*G*w == 3 else -w
        edges = []
        for a in [v,w,v-w]:
            if next(c for c in a if c)<0: a = -a
            edges.append(tuple(a))
        key = tuple(sorted(edges))
        if key in seen: continue
        seen.add(key)
        D = O+sec(v)+sec(w)
        choices = [z for z0 in zeros for z in [z0,-z0] if D*N*sec(z) == 1]
        if not choices: continue
        z = choices[0]; Z = sec(z)
        assert D*N*D == 0 and all(a*N*b == 1 for a,b in [(O,sec(v)),(O,sec(w)),(sec(v),sec(w))])
        C = matrix(ZZ,[D*N,Z*N]).right_kernel_matrix(); H = -C*N*C.transpose()
        U = matrix(ZZ,pari(H).qflllgram()); assert abs(U.det()) == 1
        C = U.transpose()*C; H = -C*N*C.transpose(); assert H.det() == 948
        result = pari(H).qfminim(2,10000,2); roots = matrix(ZZ,result[2]).transpose()
        assert int(result[0]) == 2*roots.nrows() and all(r*H*r == 2 for r in roots.rows())
        rank = 17-roots.rank()
        record = {'index':len(roster),'P_word':list(map(int,v)),'Q_word':list(map(int,w)),
                  'Z_word':list(map(int,z)),'rank':int(rank),'root_count':int(result[0])}
        roster.append(record)
        if rank != 14: continue
        exact,nodes = exact_roots(H)
        assert exact == {tuple(a) for r in roots.rows() for a in [r,-r]}
        assert roots.nrows() == 4
        # Primitive root span gives trivial torsion and a full abstract quotient.
        root_basis = matrix(ZZ,[r for r in roots.rows()]).row_module().basis_matrix()
        elementary = list(root_basis.smith_form()[0].diagonal())
        assert elementary == [1,1,1]
        K,f,c,zmap,meet = make_pencil(QQ,record)
        model = dict(record,frame_basis=rows(C),frame_gram=rows(H),root_vectors=rows(roots),
                     exact_root_nodes=nodes,fibre_class=list(map(int,D)),zero_class=list(map(int,Z)),
                     raw_cubic_coefficients=[[rat(a) for a in f[i].list()] for i in range(4)],
                     pencil_coefficients=list(map(rat,c)),zero_parameter=rat(zmap),intersection_parameters=list(map(str,meet)))
        # The old two certified pencils are exact transport/compilation controls.
        if record['index'] in [0,12]:
            name = 'elkies-k3-curve302-triangle'+('' if record['index'] == 0 else '-6-8')+'-mw14-qq-v1.json'
            control = json.loads((ART/name).read_text())
            def legacy(a): return {'n':a['numerator_coefficients_low_to_high'],'d':a['denominator_coefficients_low_to_high']}
            assert model['pencil_coefficients'] == list(map(legacy,control['pencil_coefficients']))
            control_record = dict(record,Z_word=[int(i == 6) for i in range(17)])
            control_zero = make_pencil(QQ,control_record)[3]
            assert rat(control_zero) == legacy(control['zero_parameter'])
            model['control_source'] = str((ART/name).relative_to(ROOT))
        models.append(model)
        print('PREPARED',record['index'],flush=True)
    assert (len(words),len(six),len(zeros),len(roster),len(models)) == (289,65,67,123,10)
    assert [m['index'] for m in models] == [0,12,17,30,33,56,84,98,108,120]
    return {'schema':'curve302.shortword-triangles.protocol.v1','status':'TEN_EXACT_MW14_PENCILS_PREPARED',
            'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),SOURCE,RANK_HELPER,FIBRE_HELPER,PUBLIC]},
            'word_count':len(words),'height_six_words':len(six),'height_four_zero_words':len(zeros),
            'roster':roster,'models':models,'primes':PRIMES,
            'limits':{'sample_parameters':64,'required_good_samples':52,'seconds_per_prime':300,'run_seconds':1200,'workers':1},
            'boundary':'Finite dictionary, deduplicated only by its three unoriented edge vectors. Not a complete section ball or fibration classification. Branch equality is not an equivalence proof; all-rational inverse exclusions are per compiled pencil only.'}


def model(index):
    protocol = json.loads(PROTOCOL.read_text())
    for p,h in protocol['input_sha256'].items(): assert digest(ROOT/p) == h,p
    return next(m for m in protocol['models'] if m['index'] == index)


def probe(index,prime,check=False):
    m = model(index); field = GF(prime); K,f,c,z,meet = make_pencil(field,m)
    R = K.ring()
    assert list(c) == [decode(field,R,a) for a in m['pencil_coefficients']]
    assert z == decode(field,R,m['zero_parameter'])
    assert list(meet) == list(map(lambda a:field(QQ(a)),m['intersection_parameters']))
    assert all(list(f[i]) == [decode(field,R,a) for a in m['raw_cubic_coefficients'][i]] for i in range(4))
    fibre_j = runpy.run_path(str(FIBRE_HELPER))['fibre_j']
    checkpoint = LOCAL/f'row{index}-p{prime}-samples.json'
    samples = []; skipped = []; start = 0
    if not check and checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        assert saved['protocol_sha256'] == digest(PROTOCOL)
        samples = saved['samples']; skipped = saved['skipped']; start = saved['last_parameter']+1
    for value in range(start,64):
        if len(samples) == 52: break
        j = fibre_j(field,K,f,value)
        if j is None: skipped.append(value)
        else: samples.append([value,j])
        if not check:
            checkpoint.write_text(json.dumps({'protocol_sha256':digest(PROTOCOL),'samples':samples,'skipped':skipped,'last_parameter':value},sort_keys=True)+'\n')
        print('SAMPLE',index,prime,value,j,'good',len(samples),flush=True)
    assert len(samples) == 52, 'Insufficient samples: UNKNOWN'
    S = PolynomialRing(field,'s'); s = S.gen()
    M = matrix(field,[[field(a)**i for i in range(25)]+[-field(j)*field(a)**i for i in range(25)] for a,j in samples[:49]])
    kernel = M.right_kernel(); v = next(v for v in kernel.basis() if any(v[25:]))
    n,d = S(list(v[:25])),S(list(v[25:])); common = n.gcd(d); n,d = n//common,d//common
    assert max(n.degree(),d.degree()) == 24, 'Degree loss: UNKNOWN'
    assert all(n(a) == field(j)*d(a) and d(a) for a,j in samples)
    public = runpy.run_path(str(PUBLIC)); target = EllipticCurve(QQ,list(map(QQ,public['GENERAL_WEIERSTRASS_COEFFICIENTS'])))
    assert field(target.discriminant()) != 0
    comparison = n-field(target.j_invariant())*d
    roots = [int(a) for a in field if comparison(a) == 0]; infinity = comparison[24] == 0
    return {'schema':'curve302.shortword-triangle.inverse.v1','index':index,'prime':prime,
            'status':'EXCLUDED_ALL_RATIONAL_PARAMETERS' if not roots and not infinity else 'MODULAR_SURVIVOR',
            'protocol_sha256':digest(PROTOCOL),'samples':samples,'skipped':skipped,
            'j_numerator':list(map(int,n.list())),'j_denominator':list(map(int,d.list())),
            'comparison':list(map(int,comparison.list())),'finite_roots':roots,'infinity_root':bool(infinity),
            'QQ_pencil_reduction_checked':True,'degree':24}


def result_path(index,prime): return LOCAL/f'row{index}-p{prime}.json'


def run():
    start = time.monotonic(); results = []
    for m in json.loads(PROTOCOL.read_text())['models']:
        index = m['index']
        if index in [0,12]: continue
        for prime in PRIMES:
            if time.monotonic()-start >= 1200: break
            out = result_path(index,prime)
            if out.exists():
                d = json.loads(out.read_text()); assert d['protocol_sha256'] == digest(PROTOCOL)
                if d['status'] == 'EXCLUDED_ALL_RATIONAL_PARAMETERS': break
                continue
            log = LOCAL/f'row{index}-p{prime}.log'
            print('START',index,prime,flush=True)
            with log.open('a') as stream:
                try:
                    p = subprocess.run(['sage','-python',str(Path(__file__)),'probe','--index',str(index),'--prime',str(prime)],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=min(310,1200-(time.monotonic()-start)))
                    code = p.returncode
                except subprocess.TimeoutExpired:
                    code = 'TIMEOUT'
            results.append({'index':index,'prime':prime,'exit_code':code,'log':str(log.relative_to(ROOT))})
            (LOCAL/'supervision.json').write_text(json.dumps({'protocol_sha256':digest(PROTOCOL),'attempts':results,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
            if out.exists():
                result = json.loads(out.read_text()); print('RESULT',index,prime,result['status'],result['finite_roots'],flush=True)
                if result['status'] == 'EXCLUDED_ALL_RATIONAL_PARAMETERS': break
            else: print('UNRESOLVED',index,prime,code,flush=True)


def report():
    records = []
    for m in json.loads(PROTOCOL.read_text())['models']:
        index = m['index']
        if index in [0,12]:
            suffix = '' if index == 0 else '-6-8'
            p = ART/f'elkies-k3-curve302-triangle{suffix}-mw14-mod1013-v1.json'
            old = json.loads(p.read_text())
            assert not old['finite_target_roots'] and not old['infinity_possible'] and old['j_degree'] == 24
            records.append({'index':index,'status':'EXCLUDED_ALL_RATIONAL_PARAMETERS','control_source':str(p.relative_to(ROOT)),'control_sha256':digest(p)})
            continue
        results = [json.loads(result_path(index,p).read_text()) for p in PRIMES if result_path(index,p).exists()]
        assert all(r['protocol_sha256'] == digest(PROTOCOL) for r in results)
        excluded = any(r['status'] == 'EXCLUDED_ALL_RATIONAL_PARAMETERS' for r in results)
        records.append({'index':index,'status':'EXCLUDED_ALL_RATIONAL_PARAMETERS' if excluded else 'UNKNOWN','results':results})
    return {'schema':'curve302.shortword-triangles.v1','protocol_sha256':digest(PROTOCOL),'records':records,
            'excluded':sum(r['status'] == 'EXCLUDED_ALL_RATIONAL_PARAMETERS' for r in records),
            'unknown':sum(r['status'] == 'UNKNOWN' for r in records),
            'boundary':'Only these ten exact rank14 triangle pencils. Any failed or missing computation remains unknown. No global parent exclusion or new302 parent.'}


def replay():
    assert make_protocol() == json.loads(PROTOCOL.read_text())
    saved = json.loads(OUT.read_text()); assert report() == saved
    public = runpy.run_path(str(PUBLIC)); target = EllipticCurve(QQ,list(map(QQ,public['GENERAL_WEIERSTRASS_COEFFICIENTS'])))
    for row in saved['records']:
        for r in row.get('results',[]):
            field = GF(r['prime']); R = PolynomialRing(field,'s')
            n,d = R(r['j_numerator']),R(r['j_denominator'])
            assert n.gcd(d) == 1 and max(n.degree(),d.degree()) == 24
            assert len(r['samples']) == 52 and len({a for a,j in r['samples']}) == 52
            assert all(n(a) == field(j)*d(a) and d(a) for a,j in r['samples'])
            equation = n-field(target.j_invariant())*d
            assert list(map(int,equation.list())) == r['comparison']
            assert [int(a) for a in field if equation(a) == 0] == r['finite_roots']
            assert bool(equation[24] == 0) == r['infinity_root']
    print('PASS_REPLAY_QQ_AND_SAVED_MODULAR_POLYNOMIALS',saved['excluded'],saved['unknown'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command',choices=['prepare','probe','run','report','replay'])
    parser.add_argument('--index',type=int); parser.add_argument('--prime',type=int,choices=PRIMES); parser.add_argument('--check',action='store_true'); args = parser.parse_args()
    LOCAL.mkdir(parents=True,exist_ok=True)
    signal.alarm(300 if args.command == 'probe' else 1250 if args.command == 'run' else 120)
    if args.command == 'prepare':
        assert not PROTOCOL.exists(); PROTOCOL.write_text(json.dumps(make_protocol(),sort_keys=True,indent=2)+'\n')
    elif args.command == 'probe':
        result = probe(args.index,args.prime,args.check); path = result_path(args.index,args.prime)
        if args.check: assert result == json.loads(path.read_text())
        else:
            assert not path.exists(); path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        print(result['status'],flush=True)
    elif args.command == 'run': run()
    elif args.command == 'report':
        assert not OUT.exists(); result = report(); OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n'); print(result['excluded'],result['unknown'],flush=True)
    else: replay()
