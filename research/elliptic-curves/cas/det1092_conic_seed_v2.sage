#!/usr/bin/env sage-python
"""Produce a proved rank18 seed from an integer n or recognize an input t.

Use --n for the unconditional progression; --parameter returns UNKNOWN
outside its image. Only algebraic maps are evaluated: no rational-point
search, exceptional point, catalogue rank, or amplifier input is used.
"""
import argparse,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_conic_seed_progression_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def solve_degree_two(g):
    assert g and g.degree()<=2
    if g.degree()==0:return []
    if g.degree()==1:return [-g[0]/g[1]]
    delta=g[1]**2-4*g[2]*g[0]
    if delta<0 or not delta.is_square():return []
    root=delta.sqrt()
    return sorted(set([(-g[1]+root)/(2*g[2]),(-g[1]-root)/(2*g[2])]))
def produce(args):
    parent_path=ART/'curve302_recovered_mw17_parent_v1.json'
    cover_path=ART/'det1092_orbit8044_rank18_base_change_v2.json'
    proof_path=DIR/'replay.json'
    parent,cover,proof=map(read,[parent_path,cover_path,proof_path])
    provenance(cover);provenance(proof)
    # Bind current generic maps transitively to the already verified proof.
    provenance(read(DIR/'protocol.json'))
    assert proof['status']=='PASS_INDEPENDENT_ORACLE_FREE_RANK18_PROGRESSION'
    B=ZZ(proof['B'])
    R=PolynomialRing(QQ,'t');K=R.fraction_field();U=PolynomialRing(QQ,'u');F=U.fraction_field()
    def dec(d,R,K):return K(R(d['numerator']))/R(d['denominator'])
    T=dec(cover['parametrization']['t_of_u'],U,F)
    if args.n is not None:
        n=ZZ(args.n);u=B*n;tau=QQ(T(u))
    else:
        tau=QQ(args.parameter)
        equation=T.numerator()-tau*T.denominator()
        roots=solve_degree_two(equation)
        accepted=[v for v in roots if T.denominator()(v) and v/B in ZZ]
        if not accepted:
            return {'status':'UNKNOWN_OUTSIDE_PROVED_PROGRESSION_IMAGE',
                'classification':'exact noncoverage of this construction; not a seed-nonexistence claim',
                'parameter':str(tau),'rational_preimages':[str(v) for v in roots],
                'condition':'Some rational root u of num(t(u))-t*den(t(u)) must have u/B an integer.',
                'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [parent_path,cover_path,proof_path,Path(__file__)]}}
        u=accepted[0];n=ZZ(u/B);assert T(u)==tau
    def at_t(d):return QQ(dec(d,R,K)(tau))
    E=EllipticCurve(QQ,[at_t(d) for d in parent['a_invariants']])
    assert E.discriminant()
    generic=[E([at_t(d) for d in P]) for P in parent['basis_weierstrass_coordinates']]
    seed=E([dec(cover['lift'][key],U,F)(u) for key in ['x_of_u','y_of_u']])
    S=dec(cover['parametrization']['s_of_u'],U,F)
    q=R(cover['curve_over_Q']['q_coefficients'])
    assert S(u)**2==q(tau)
    # The proof applies to every integer n, without another point/rank search.
    assert n in ZZ and u==B*n
    return {'status':'PROVED_RANK_AT_LEAST_18_WITH_NON_GENERIC_SEED',
        'classification':'constructive verified application of the uniform progression theorem',
        'parameter':str(tau),'n':str(n),'conic_parameter_u':str(u),'B':str(B),
        'elliptic_a_invariants':list(map(str,E.a_invariants())),
        'generic_points':[list(map(str,P)) for P in generic],
        'seed':list(map(str,seed)),
        'cover_point':[str(tau),str(S(u))],
        'certified_subgroup_rank_lower_bound':18,
        'independence':'All18 reductions agree with the independent rank18 finite-group certificate, since u=B*n. The reduction at23 has odd order, excluding rational2-torsion.',
        'scope':'No full-rank upper bound, new-record claim, catalogue comparison or amplification.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [parent_path,cover_path,proof_path,Path(__file__)]}}
if __name__=='__main__':
    signal.alarm(25)
    parser=argparse.ArgumentParser();choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--n');choice.add_argument('--parameter');parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=produce(args);payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        if args.output.exists():assert args.output.read_text()==payload
        else:args.output.write_text(payload)
        print(result['status'],str(args.output),flush=True)
    else:print(payload,flush=True)
