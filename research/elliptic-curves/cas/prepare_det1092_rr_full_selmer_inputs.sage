#!/usr/bin/env sage-python
"""Prepare one bounded, equation-only full-Selmer input. No descent is run."""
import argparse,hashlib,json,signal,shutil
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,gcd,lcm,prime_range,pari
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_full_selmer_inputs_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def source_paths():
    return [ART/'det1092_rr_generic_point_controls_v2'/('case-%02d.json'%i) for i in range(9)]+[
        ART/'det1092_rr_residual_jacobian_class_v1.json']
def prepare(i):
    paths=source_paths();protocol={'classification':'frozen input preparation, not a completed Selmer experiment',
        'cohort':'Nine source-only r0(t) family members; the historical first-unlock RR member is a separate retrospective calibration arm.',
        'comparison_boundary':'Masking computation inputs does not make the historically chosen first-unlock RR member selection-blind.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'limits':{'wall_seconds_per_case':25,'cases':10,'trial_primes_through':97,
                  'model_reductions_per_case':1,'point_searches':0,'Selmer_runs':0,'class_group_runs':0,
                  'full_integer_factorizations':0,'pilot_changes':0},
        'execution_gate':'Do not run full descent until a suitable backend and explicit class-group/resource cap are available.',
        'script_sha256':sha(Path(__file__))}
    OUT.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    raw=json.loads(paths[i].read_text())
    curve=raw['curve'] if i==9 else raw
    R=PolynomialRing(QQ,'x');x=R.gen();q=R(curve['q']);scale=QQ(curve['scale'])
    assert q.degree()==6 and q.gcd(q.derivative())==1
    # Primitive field equation and exact scalar, without factoring an integer.
    den=lcm(a.denominator() for a in q);a=R(q*den);content=gcd(ZZ(v) for v in a)
    g=a/content
    if g.leading_coefficient()<0:g=-g;content=-content
    lam=scale*content/den;assert scale*q==lam*g
    small=list(prime_range(2,98))
    def strip(n):
        n=ZZ(n);square=ZZ(1);rows=[]
        for p in small:
            e=n.valuation(p);n//=p**e;square*=p**(e//2)
            rows.append({'p':int(p),'valuation':int(e)})
        # Retain the odd valuations; only extract a cofactor if it is exactly
        # a square. No probable-prime or unfactored-cofactor assumption.
        odd=ZZ(1)
        for row in rows:odd*=ZZ(row['p'])**(row['valuation']%2)
        root=n.sqrt() if n.is_square() else ZZ(1)
        residual=odd*(n//(root*root));square*=root
        return square,residual,rows
    ns,N,ntrial=strip(abs(lam.numerator()));ds,D,dtrial=strip(lam.denominator())
    yscale=QQ(ns)/(ds*D);F=R(lam.sign()*N*D*g)
    assert all(a.denominator()==1 for a in F) and yscale*yscale*F==scale*q
    prefix={'classification':'exact model preparation; full Selmer status remains NOT_COMPUTED',
            'case_index':i,'arm':'RETROSPECTIVE_CALIBRATION' if i==9 else 'GENERIC_SOURCE_ONLY_COHORT',
            'source_path':str(paths[i].relative_to(ROOT)),'source_sha256':sha(paths[i]),
            'primitive_sextic':list(map(str,g.list())),'scalar':str(lam),
            'square_strip_numerator':ntrial,'square_strip_denominator':dtrial,
            'pre_reduction_polynomial':list(map(str,F.list())),'s_equals_y_scale_times_old_y':str(yscale),
            'software':{'sage':sage_version,'pari':str(pari.version())}}
    retain(OUT/('case-%02d-integral.json'%i),prefix)
    # Cremona--Stoll reduction changes coordinates, not the curve or its
    # discriminant. It does not enumerate points, fields, or classes.
    reduction=pari('(f)->{my(m);my(r=hyperellred(f,&m));[r,m]}')(pari(F))
    PQ,m=reduction[0],reduction[1];P,Q=R(PQ[0]),R(PQ[1])
    e=QQ(m[0]);M=m[1];aa,bb,cc,dd=[QQ(M[j,k]) for j,k in [(0,0),(0,1),(1,0),(1,1)]]
    H=R(m[2]);assert e==1 and aa*dd-bb*cc==1
    transform=sum(F[j]*(aa*x+bb)**j*(cc*x+dd)**(6-j) for j in range(7))
    assert 2*H==Q and transform==P+H*H
    primitive_disc=ZZ(g.discriminant());cofactor=abs(primitive_disc);trial=[]
    for p in small:
        v=cofactor.valuation(p);cofactor//=p**v
        trial.append({'p':int(p),'valuation':int(v)})
    worker={'schema':'det1092.equation-only-hyperelliptic-input.v1',
            'P':list(map(str,P.list())),'Q':list(map(str,Q.list()))}
    opaque=hashlib.sha256(json.dumps(worker,sort_keys=True).encode()).hexdigest()[:20]
    retain(OUT/('input-'+opaque+'.json'),worker)
    result={**prefix,'status':'PASS_EXACT_INTEGRAL_MODEL_AND_EQUATION_ONLY_INPUT',
            'worker_input':'input-'+opaque+'.json','reduced_P':worker['P'],'reduced_Q':worker['Q'],
            'new_to_integral_model':{'mobius':list(map(str,[aa,bb,cc,dd])),'H':list(map(str,H.list()))},
            'new_to_original':'T=(a*x+b)/(c*x+d); s=yscale*(y+H(x))/(c*x+d)^3',
            'primitive_polynomial_discriminant':str(primitive_disc),
            'disc_trial_division':trial,'unfactored_disc_cofactor':str(cofactor),
            'disc_cofactor_bits':int(cofactor.nbits()),
            'coefficient_bits_before':max(int(abs(ZZ(a)).nbits()) for a in F),
            'coefficient_bits_after':max(int(abs(ZZ(a)).nbits()) for a in list(P)+list(Q)),
            'arithmetic_boundary':'Polynomial discriminant is not a certified field discriminant or complete bad-prime support. The retained cofactor is not assumed prime or squarefree.',
            'backend':{'magma_executable':shutil.which('magma'),'full_Selmer_group':None,
                       'class_group_completeness':None,'other_Sha_classes':None},
            'limits':protocol['limits'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [paths[i],OUT/'protocol.json',Path(__file__)]}}
    retain(OUT/('case-%02d.json'%i),result)
    print('PASS case',i,'coefficient bits',result['coefficient_bits_before'],result['coefficient_bits_after'],
          'unfactored discriminant cofactor bits',result['disc_cofactor_bits'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);prepare(args.case)
