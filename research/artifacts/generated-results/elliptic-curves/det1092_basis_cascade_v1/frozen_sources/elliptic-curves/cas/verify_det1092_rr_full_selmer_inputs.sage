#!/usr/bin/env sage-python
"""Replay integral maps, exact model discriminants and worker input isolation.

No PARI model reduction, point search or Selmer calculation is invoked.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,gcd,lcm,prime_range
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_rr_full_selmer_inputs_v1'
OUT=DIR/'replay.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    pp=DIR/'protocol.json';protocol=json.loads(pp.read_text())
    for p,h in protocol['inputs'].items():assert sha(ROOT/p)==h
    paths=[ART/'det1092_rr_generic_point_controls_v2'/('case-%02d.json'%i) for i in range(9)]+[
        ART/'det1092_rr_residual_jacobian_class_v1.json']
    R=PolynomialRing(QQ,'x');x=R.gen();primes=list(prime_range(2,98));rows=[];inputs=[pp]
    for i,path in enumerate(paths):
        cp=DIR/('case-%02d.json'%i);d=json.loads(cp.read_text());raw=json.loads(path.read_text())
        curve=raw['curve'] if i==9 else raw;q=R(curve['q']);scale=QQ(curve['scale'])
        assert d['source_path']==str(path.relative_to(ROOT)) and d['source_sha256']==sha(path)
        for p,h in d['inputs'].items():assert sha(ROOT/p)==h
        assert d['case_index']==i and d['status']=='PASS_EXACT_INTEGRAL_MODEL_AND_EQUATION_ONLY_INPUT'
        assert d['arm']==('RETROSPECTIVE_CALIBRATION' if i==9 else 'GENERIC_SOURCE_ONLY_COHORT')
        g=R(d['primitive_sextic']);lam=QQ(d['scalar']);k=QQ(d['s_equals_y_scale_times_old_y'])
        assert g.degree()==6 and all(a.denominator()==1 for a in g)
        assert gcd(ZZ(a) for a in g)==1 and g.leading_coefficient()>0
        assert lam*g==scale*q
        for n,trial in [(abs(lam.numerator()),d['square_strip_numerator']),
                        (lam.denominator(),d['square_strip_denominator'])]:
            assert [row['p'] for row in trial]==primes
            assert all(n.valuation(row['p'])==row['valuation'] for row in trial)
        F=R(d['pre_reduction_polynomial']);assert k*k*F==scale*q and k>0
        assert all(a.denominator()==1 for a in F)
        P,Q=R(d['reduced_P']),R(d['reduced_Q'])
        assert P.degree()<=6 and Q.degree()<=3 and all(a.denominator()==1 for a in list(P)+list(Q))
        a,b,c,e=map(QQ,d['new_to_integral_model']['mobius']);H=R(d['new_to_integral_model']['H'])
        assert a*e-b*c==1 and all(z.denominator()==1 for z in [a,b,c,e])
        assert 2*H==Q
        # Verify the complete pullback relation, not a sample of points.
        pullback=sum(scale*q[j]*(a*x+b)**j*(c*x+e)**(6-j) for j in range(7))
        assert pullback==k*k*(P+H*H)
        # Direct Sylvester determinant independently checks the saved
        # discriminant, which the constructor obtained via discriminant().
        disc=-g.sylvester_matrix(g.derivative()).determinant()/g.leading_coefficient()
        assert disc and ZZ(disc)==ZZ(d['primitive_polynomial_discriminant'])
        trial=d['disc_trial_division'];assert [row['p'] for row in trial]==primes
        rem=ZZ(d['unfactored_disc_cofactor']);assert rem>0
        for row in trial:
            p=ZZ(row['p']);v=ZZ(row['valuation']);assert v>=0 and rem%p!=0
            rem*=p**v
        assert rem==abs(disc)
        assert ZZ(d['unfactored_disc_cofactor']).nbits()==d['disc_cofactor_bits']
        worker_path=DIR/d['worker_input'];worker=json.loads(worker_path.read_text())
        assert set(worker)=={'schema','P','Q'}
        assert worker['schema']=='det1092.equation-only-hyperelliptic-input.v1'
        assert worker['P']==d['reduced_P'] and worker['Q']==d['reduced_Q']
        opaque=hashlib.sha256(json.dumps(worker,sort_keys=True).encode()).hexdigest()[:20]
        assert d['worker_input']=='input-'+opaque+'.json'
        assert all(d['backend'][key] is None for key in ['full_Selmer_group','class_group_completeness','other_Sha_classes'])
        assert max(int(abs(ZZ(a)).nbits()) for a in F)==d['coefficient_bits_before']
        assert max(int(abs(ZZ(a)).nbits()) for a in list(P)+list(Q))==d['coefficient_bits_after']
        rows.append({'case_index':i,'arm':d['arm'],'worker_input':d['worker_input'],
                     'map_identity':True,'polynomial_discriminant_verified':True,
                     'coefficient_bits_after':d['coefficient_bits_after'],'unfactored_cofactor_bits':d['disc_cofactor_bits'],
                     'full_Selmer_group':'NOT_COMPUTED'})
        inputs.extend([path,cp,worker_path])
    return {'classification':'verified application: descent input preparation only',
            'status':'PASS_TEN_EXACT_MODELS_AND_ISOLATED_INPUTS_NOT_A_SELMER_COMPUTATION',
            'sage_version':sage_version,'cases':rows,
            'oracle_boundary':'The nine-member cohort is fixed by the earlier generic-section rule. The tenth first-unlock member is explicitly retrospective. Equation-only input isolation is verified, not a claim that historical curve selection was blind.',
            'completion_gate':'Full Selmer groups require a general genus-two descent backend and a separately agreed resource budget. No global class/unit computation is authorized or run by these scripts.',
            'claim_boundary':'Coefficient heights and polynomial discriminant cofactors are properties of the displayed inputs, not lower bounds on intrinsic field discriminants or proof that all descent algorithms are infeasible.',
            'limits':protocol['limits'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);d=verify();payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(d['status'],flush=True)
