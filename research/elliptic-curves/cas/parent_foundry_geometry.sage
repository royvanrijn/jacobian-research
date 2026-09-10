#!/usr/bin/env sage-python
"""Target-free construction of marked A1 neighbours on the rational X948 K3.

The proposals are genuine new pencils D=O+P_w, never deleted MW17 bases.
Exact section identities and a positive Shioda Gram certify generic rank 16;
the I2+22I1 profile and the source Picard-19 theorem certify its upper bound.
All runtime inputs contain only generic geometry. No record fibre is read.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction

from sage.all import (QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector,
                      pari, gcd, lcm, prime_range)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
sys.path.insert(0, str(CAS))
from research_runtime.store import checkpoint
import audit_r17_constant_scaling as scaling

ATLAS = ROOT/'artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json'
TABLE = ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv'


def load(name, path):
    return SourceFileLoader(name, str(path)).load_module()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(value):
    return list(map(str, value.list())) or ['0']


def encode(value):
    return {'numerator_coefficients_low_to_high': record(value.numerator()),
            'denominator_coefficients_low_to_high': record(value.denominator())}


def decode(value, ring):
    return ring.fraction_field()(ring(value['numerator_coefficients_low_to_high'])) / ring(value['denominator_coefficients_low_to_high'])


def source_class(v, gram):
    return vector(QQ, [(v*gram*v-2)/2, 1] + list(v))


def pairing(a, b, gram):
    return a[0]*b[1]+a[1]*b[0]-vector(QQ,a[2:])*gram*vector(QQ,b[2:])


def section_candidates(gram, trace):
    """Bounded floating shell proposals, each tested by the exact norm identity.

    This is deliberately not used to assert completeness of a shell. Finding
    sixteen sections with a positive exact Gram is the acceptance condition.
    Adapted from the historical degree-one shell enumerator (2026-09-04).
    """
    n = gram.nrows()
    change = gram.LLL_gram().transpose()
    reduced = change*gram*change.transpose()
    parity = tuple(int(x)%2 for x in trace*change.inverse())
    lower = matrix(QQ,n,n); diagonal = []
    for i in range(n):
        lower[i,i] = 1
        diagonal.append(reduced[i,i]-sum(lower[i,k]**2*diagonal[k] for k in range(i)))
        for j in range(i+1,n):
            lower[j,i] = (reduced[j,i]-sum(lower[j,k]*lower[i,k]*diagonal[k] for k in range(i)))/diagonal[i]
    lf = [[float(x) for x in row] for row in lower.rows()]
    df = list(map(float,diagonal)); coords = [0]*n; answer = []
    def visit(i, cost):
        if i<0:
            z = vector(ZZ,coords)*change
            if z*gram*z == 12:
                v = vector(ZZ,(trace-z)/2)
                if v*gram*v-trace*gram*v != 1:
                    raise ArithmeticError('degree-one identity failed')
                answer.append(v)
            return
        center = sum(lf[j][i]*coords[j] for j in range(i+1,n))
        radius = math.sqrt(max(0,(12.25-cost)/df[i]))
        lo, hi = math.ceil(-center-radius-1e-9), math.floor(-center+radius+1e-9)
        for value in range(lo+(parity[i]-lo)%2,hi+1,2):
            extra = df[i]*(value+center)**2
            if cost+extra<=12.25+1e-9:
                coords[i] = value; visit(i-1,cost+extra)
    visit(n-1,0)
    return sorted(answer,key=lambda v:(sum(abs(x) for x in v),tuple(v)))


def shioda(vectors, zero, trace, gram):
    """Intersection formula with the unique A1 correction; zero meets P_w."""
    if zero*gram*zero != 4:
        raise ArithmeticError('zero must meet the identity component')
    z = source_class(zero,gram)
    o = vector(QQ,[-1,1]+[0]*17)
    data = [(source_class(v,gram),pairing(source_class(v,gram),z,gram),
             pairing(source_class(v,gram),o,gram)) for v in vectors]
    if any(c not in (0,1) for _,_,c in data):
        raise ArithmeticError('section component index is not 0 or 1')
    return matrix(QQ,len(vectors),len(vectors),
        lambda i,j: (4+2*data[i][1]-data[i][2]/2) if i==j else
        2+data[i][1]+data[j][1]-pairing(data[i][0],data[j][0],gram)-data[i][2]*data[j][2]/2)


def affine_invariant(a, b):
    """Exact separating invariant of j's poles, for a unique I2 at infinity.

    Every fibration equivalence must preserve that unique double pole, hence
    acts affinely here. Depress the monic degree-22 finite discriminant.
    Ratios of equal weights are invariant under every affine base change.
    Different keys prove inequivalence even over Qbar. Equal keys are only
    a possible collision and are never counted as a new fibration.
    """
    d = a.parent()(4*a**3+27*b**2)
    if d.degree()!=22 or d.gcd(d.derivative()).degree()!=0:
        raise ArithmeticError('unique-I2 invariant requires 22 simple finite poles')
    t = d.parent().gen(); shift = -d[21]/(22*d[22])
    p = d(t+shift)/d[22]
    if not p[20] or not p[19]:
        raise ArithmeticError('separating chart unresolved: weight 2 or 3 vanishes')
    values = []
    for i in range(21):
        weight = 22-i
        values.append(str(p[i]/p[20]**(weight//2) if weight%2==0 else
                          p[i]*p[19]/p[20]**((weight+3)//2)))
    return {'key':hashlib.sha256(json.dumps(values).encode()).hexdigest(),
            'values':values, 'method':'unique-I2 depressed discriminant weight ratios v1'}


def compact(a, b, hints):
    """Bounded auxiliary minimization; exact weighted identities decide validity."""
    values = lambda f:[Fraction(str(x)) for x in f.list()]
    u0,_ = scaling.scale_for(values(a),values(b)); u0=QQ(str(u0))
    aa,bb = a/u0**4,b/u0**6
    def primitive(f):
        g=f*lcm(q.denominator() for q in f)
        return g/gcd(ZZ(q) for q in g)
    p,q=primitive(aa),primitive(bb)
    common=abs(gcd([ZZ(p.discriminant()),ZZ(q.discriminant()),ZZ(p.resultant(q))]))
    remainder=common; primes=[]
    for prime in prime_range(10000):
        if remainder%prime==0:
            primes.append(prime); remainder//=prime**remainder.valuation(prime)
    exponent,root=1,remainder
    if remainder>1:
        for k in range(2,133):
            candidate,exact=remainder.nth_root(k,truncate_mode=True)
            if exact: exponent,root=k,candidate
    reducer=load('parent_foundry_compact_helpers',CAS/'compact_mw16_base_v2.sage')
    factors,pieces=reducer.bounded_factor(root,{
        'maximum_power_exponent':132,'prime_root_bit_gate':512,
        'composite_factor_bit_gate':192},hints)
    primes=sorted(set(primes+[v for v,e in factors]))
    if root**exponent!=remainder or math.prod(v**e for v,e in factors)!=root:
        raise ArithmeticError('factor hint reconstruction failed')
    pari.allocatemem(256000000,silent=True); pari.set_real_precision_bits(8192)
    fn=pari('(P,L)->{my(m,n);my(Q=hyperellminimalmodel(P,&m,L));my(S=hyperellred(Q,&n));[m,n]}')
    m,n=fn(pari(p),pari(primes))
    maps=[matrix(QQ,2,2,[v[1][i,j] for i,j in ((0,0),(0,1),(1,0),(1,1))]) for v in (m,n)]
    transform=maps[0]*maps[1]
    helper=load('parent_foundry_homogeneous',CAS/'reduce_r17_family_base.sage')
    x,y,z,w=transform.list()
    A,B=helper.homogeneous(aa,8,x,y,z,w),helper.homogeneous(bb,12,x,y,z,w)
    extra,_=scaling.scale_for(values(A),values(B)); extra=QQ(str(extra))
    A,B,u=A/extra**4,B/extra**6,u0*extra
    if not transform.det() or helper.homogeneous(a,8,x,y,z,w)!=u**4*A or helper.homogeneous(b,12,x,y,z,w)!=u**6*B:
        raise ArithmeticError('exact compactification identity failed')
    return A,B,transform,u,{'before_bits':scaling.bits(values(a)+values(b)),
        'after_bits':scaling.bits(values(A)+values(B)),'primes':list(map(str,primes)),
        'factor_pieces':pieces,'base_matrix_a_b_c_d':list(map(str,transform.list())),
        'scale_u':str(u)}


def construct(priority, output, source=None, trace_coordinates=None):
    atlas=json.loads(ATLAS.read_text())
    native=source or next(f for f in atlas['families'] if f['family']=='11952')
    if trace_coordinates is None:
        rows=list(csv.DictReader(TABLE.open(),delimiter='\t'))
        row=rows[priority-1]
        if int(row['priority_rank'])!=priority or int(row['minimal_unoriented_count'])!=1:
            raise ArithmeticError('construction must select the proved A1 stratum')
        trace_coordinates=list(map(int,row['section_basis_w'].split()))
    trace_word=vector(ZZ,trace_coordinates)
    gram=matrix(QQ,native['generic_height_gram'])
    determinant=gram.det()
    if trace_word*gram*trace_word!=8 or determinant not in (948,1092):
        raise ArithmeticError('marked source lattice changed')
    if source is not None:
        import numpy as np
        from visibility_lattice_fast import IntegerExactParity
        change=gram.LLL_gram().transpose()
        solver=IntegerExactParity((2*change*gram*change.transpose()).rows())
        residue=tuple(int(x)%2 for x in trace_word*change.inverse())
        starts,_=solver.babai(np.asarray([residue],dtype=np.int64))
        proof=solver.solve(residue,tuple(map(int,starts[0])),2000000)
        if proof['norm']!=16 or len(proof['minima'])!=2:
            output.mkdir(parents=True,exist_ok=True)
            checkpoint(output/'not-a1.json',{'status':'NOT_A1_STRATUM',
                'source_family':native['family'],'priority':priority,'trace_word':trace_coordinates,
                'scaled_minimum':proof['norm'],'signed_minima':len(proof['minima']),
                'boundary':'Rejects this norm-eight A1 construction proposal, not any specialized elliptic rank.'})
            return None
    T=PolynomialRing(QQ,'t'); K=T.fraction_field()
    a,b=(T(native[k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    E=EllipticCurve(K,[a,b])
    basis=[E(decode(s['X'],T),decode(s['Y'],T)) for s in native['sections']]
    def combination(word):
        return sum((n*p for n,p in zip(word,basis) if n),E(0))
    trace=combination(trace_word)
    builder=load('parent_foundry_chord_builder',CAS/'build_a1_mw16_target_free_parameter_candidates.sage')
    chord=load('parent_foundry_chord',builder.CHORD)
    h,nx,m0,quartic,child_a,child_b=builder.child_geometry(trace,a,b,T,chord)
    delta=child_a.parent()(-16*(4*child_a**3+27*child_b**2))
    if [child_a.degree(),child_b.degree(),delta.degree()]!=[8,12,22] or delta.gcd(delta.derivative()).degree()!=0:
        raise ArithmeticError('not an I2+22I1 pencil')
    identity=affine_invariant(child_a,child_b)
    output.mkdir(parents=True,exist_ok=True)
    checkpoint(output/'equation.json',{'priority':priority,'trace_word':list(map(int,trace_word)),
        'A':record(child_a),'B':record(child_b),'novelty_invariant':identity})
    print('PARENT_EQUATION',priority,flush=True)
    vectors=section_candidates(gram,trace_word)
    zero=next(v for v in vectors if v*gram*v==4)
    pool=[v for v in vectors if v!=zero]
    full=shioda(pool,zero,trace_word,gram)
    chosen=list(full.pivots())[:16]
    if len(chosen)!=16:
        raise ArithmeticError('insufficient explicit section proposals')
    # Prefer a saturated displayed basis; keep the exact index if none found.
    for shift in range(min(32,len(pool))):
        order=list(range(shift,len(pool)))+list(range(shift))
        pick=[]
        for i in order:
            trial=pick+[i]
            if full.matrix_from_rows_and_columns(trial,trial).rank()==len(trial):pick=trial
            if len(pick)==16:break
        g=full.matrix_from_rows_and_columns(pick,pick)
        if g.det()==determinant/2: chosen=pick; break
    selected=[pool[i] for i in chosen]
    g=shioda(selected,zero,trace_word,gram)
    if not g.is_positive_definite() or not (g.det()/(determinant/2)).is_square():
        raise ArithmeticError('generic independence or index identity failed')
    print('PARENT_LATTICE',priority,str(g.det()),len(vectors),flush=True)
    A,B,M,u,compactification=compact(child_a,child_b,{
        'source':native,'trace_x':encode(trace[0]),'trace_y':encode(trace[1]),
        'h':record(h),'nx':record(nx),'m0':record(m0)})
    checkpoint(output/'compactification.json',compactification)
    R=PolynomialRing(QQ,'z'); z=R.gen(); L=R.fraction_field()
    aa,bb,cc,dd=map(QQ,M.list()); lam=L(aa*z+bb)/(cc*z+dd); den=L(cc*z+dd)
    S=PolynomialRing(L,'v'); v=S.gen()
    q=S([L(c(lam)) for c in quartic.list()])
    sx,rem=((S(m0)+lam*S(h)**2)**2-S(nx)).quo_rem(S(h)**2)
    if rem:raise ArithmeticError('chord sum identity failed')
    points=[]; maps=[]
    for i,word in enumerate([zero]+selected):
        p=combination(word)
        base=K((((p[1]+trace[1])/(p[0]-trace[0]))*h-m0)/h**2)
        if max(base.numerator().degree(),base.denominator().degree())!=1:
            raise ArithmeticError('old curve is not a degree-one section')
        n,d=base.numerator(),base.denominator()
        old_t=(n[0]-lam*d[0])/(lam*d[1]-n[1])
        if n(old_t)/d(old_t)!=lam:raise ArithmeticError('base inversion failed')
        w=(2*p[0](old_t)-sx(old_t))/h(old_t)
        if w**2!=q(old_t):raise ArithmeticError('quartic section identity failed')
        points.append((old_t,w)); maps.append(encode(base))
        print('PARENT_QUARTIC_SECTION',priority,i,flush=True)
    t0,w0=points[0]
    ee,d,c,b,a=[L(q(v+t0)[i]) for i in range(5)]
    if not w0 or ee!=w0**2:raise ArithmeticError('pointed zero failed')
    a1=d/w0; a2=c-d**2/(4*w0**2); a3=2*w0*b; a4=-4*w0**2*a; a6=a2*a4
    b2=a1**2+4*a2; b4=a1*a3+2*a4; b6=a3**2+4*a6
    if -81*(b2**2-24*b4)/48!=child_a(lam) or -729*(-b2**3+36*b2*b4-216*b6)/864!=child_b(lam):
        raise ArithmeticError('pointed Jacobian normalization failed')
    sections=[]
    for i,(old_t,w) in enumerate(points[1:]):
        delta=old_t-t0
        xg=(2*w0*(w+w0)+d*delta)/delta**2
        yg=(4*w0**2*(w+w0)+2*w0*d*delta+(2*w0*c-d**2/(2*w0))*delta**2)/delta**3
        x=9*(xg+b2/12)*den**4/u**2
        y=27*(yg+(a1*xg+a3)/2)*den**6/u**3
        if y**2!=x**3+R(A)*x+R(B):raise ArithmeticError('compact section identity failed')
        sections.append({'basis_index':i,'X':encode(x),'Y':encode(y)})
        print('PARENT_SECTION',priority,i,flush=True)
    result={'schema':'parent-foundry.marked-a1.v1','status':'PASS_EXACT_CONSTRUCTED_MW16',
        'family':f'x{determinant}-{native["family"]}-a1-{priority:05d}',
        'source_surface':f'X{determinant}; new fibration on the declared existing K3',
        'source_family':native['family'],'priority':priority,
        'source_sha256':native.get('source_certificate_sha256',digest(ATLAS)),
        'source_input_sha256':hashlib.sha256(json.dumps(native,sort_keys=True).encode()).hexdigest(),
        'source_generic_gram':native['generic_height_gram'],'source_determinant':str(determinant),
        'source_input':native,
        'table_sha256':digest(ROOT/'artifacts/generated-results/elliptic-curves/curve302_parent_degree2_multisection_orbits_v1.tsv')
            if native['family']=='det1092' else digest(TABLE),
        'trace_word':list(map(int,trace_word)),
        'zero_word':list(map(int,zero)),'section_words':[list(map(int,w)) for w in selected],
        'base_maps':maps,'generic_height_gram':[[str(x) for x in row] for row in g.rows()],
        'generic_height_gram_determinant':str(g.det()),'generic_rank_lower_bound':16,
        'generic_rank_upper_bound':16,'rank_completion':'THEOREM_BLOCKED_EXACT_MW16',
        'upper_bound_reason':f'Picard(X{determinant})=19 and A1 root rank 1: 19-2-1=16.',
        'fibre_configuration':'I2 + 22 I1','sections':sections,
        'A_coefficients_low_to_high':record(A),'B_coefficients_low_to_high':record(B),
        'raw_A':record(child_a),'raw_B':record(child_b),'novelty_invariant':identity,
        'compactification':compactification,
        'claim_boundary':'Exact new equation and generic sections; novelty is gated separately against known/accepted parents. No new K3, no rank-17 completion, no specialized rank claim.'}
    checkpoint(output/'parent.json',result)
    verify(output/'parent.json')
    return result


def verify(path):
    p=json.loads(path.read_text()); R=PolynomialRing(QQ,'z')
    A,B=(R(p[k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    for s in p['sections']:
        x,y=decode(s['X'],R),decode(s['Y'],R)
        if y*y!=x*x*x+A*x+B:raise ArithmeticError('section replay failed')
    source=p.get('source_input') or next(f for f in json.loads(ATLAS.read_text())['families'] if f['family']=='11952')
    if p.get('source_input_sha256') and hashlib.sha256(json.dumps(source,sort_keys=True).encode()).hexdigest()!=p['source_input_sha256']:
        raise ArithmeticError('source input binding changed')
    g=shioda([vector(ZZ,w) for w in p['section_words']],vector(ZZ,p['zero_word']),
        vector(ZZ,p['trace_word']),matrix(QQ,source['generic_height_gram']))
    if g!=matrix(QQ,p['generic_height_gram']) or not g.is_positive_definite():
        raise ArithmeticError('intersection replay failed')
    rawA,rawB=R(p['raw_A']),R(p['raw_B'])
    if affine_invariant(rawA,rawB)!=p['novelty_invariant']:
        raise ArithmeticError('novelty invariant replay failed')
    helper=load('parent_verify_homogeneous',CAS/'reduce_r17_family_base.sage')
    a,b,c,d=map(QQ,p['compactification']['base_matrix_a_b_c_d']);u=QQ(p['compactification']['scale_u'])
    if a*d==b*c or helper.homogeneous(rawA,8,a,b,c,d)!=u**4*A or helper.homogeneous(rawB,12,a,b,c,d)!=u**6*B:
        raise ArithmeticError('coefficient transport replay failed')
    checkpoint(path.parent/'verified.json',{'status':'PASS_EXACT_PARENT_REPLAY',
        'parent_sha256':digest(path),'generic_rank_lower_bound':16,
        'section_identities':16,'scope':'Section identities, marked intersections, coefficient transport and separating invariant; reconstruction provenance retained.'})
    print('PARENT_VERIFIED',p['family'],p['compactification']['after_bits'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--priority',type=int)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    parser.add_argument('--request',type=Path)
    args=parser.parse_args()
    if args.verify:verify(args.verify)
    elif args.request:
        request=json.loads(args.request.read_text())
        construct(request['priority'],args.output,request['source'],request['trace_word'])
    else:construct(args.priority,args.output)
