#!/usr/bin/env sage-python
"""Frozen equation-pair transport gate: ten pairs, old64 primes,25 seconds.

No exceptional point coordinate is read. The first carrier equation was
historically calibrated; all nine controls are the unchanged old-S0 rule.
No isogeny, number-field, Selmer or point search is invoked.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_genus1_transport_gate_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==s,('IMMUTABLE_MISMATCH',str(p))
    else:p.write_text(s)
def residue(v,p):return int(GF(p)(QQ(v)))
def main():
    parentpath=ART/'curve302_recovered_mw17_parent_v1.json'
    oldprotocol=ART/'det1092_genus1_picard_controls_v1/protocol.json'
    selectionpath=ART/'det1092_genus1_picard_controls_v1/selection.json'
    firstpath=ART/'det1092_genus1_picard_image_v1/case-01-generic.json'
    controls=[ART/f'det1092_genus1_picard_controls_v1/case-{i:02d}-generic.json' for i in range(9)]
    paths=[parentpath,oldprotocol,selectionpath,firstpath,*controls,Path(__file__)]
    parent,prior,selection,first=map(read,paths[:4])
    protocol={'classification':'equation-only transport-obstruction test on already selected members',
      'rule':'Use the fixed first carrier equation at original t=0 and all nine unchanged source-only carriers at their frozen addresses. Use every prime in the old64-prime pool. Record all outcomes; report first witnesses for irreducibility, transpositions, unequal quadratic characters, and distinct ordinary Frobenius fields. No replacement pair, prime, point or member.',
      'limits':{'wall_seconds':25,'pairs':10,'primes':64,'new_addresses':0,
        'exceptional_coordinate_inputs':0,'point_searches':0,'number_fields':0,
        'class_groups':0,'Selmer_groups':0,'V3_inputs':0,'pilot_changes':0},
      'primes':prior['primes'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(OUT/'protocol.json',protocol)
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']))
    assert ai[:3]==[1,1,1]
    cases=[('302-first-carrier',QQ(0),firstpath)]+[(r['label'],QQ(r['tau']),p) for r,p in zip(selection['cases'],controls)]
    results=[]
    for i,(label,tau,path) in enumerate(cases):
        data=read(path);q=R(data['quartic']);t0,s=map(QQ,data['basepoint'])
        q0,q1,q2,q3,q4=q(t+t0).list();assert q0==s*s and s and q.gcd(q.derivative())==1
        J=[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        assert J==list(map(QQ,data['Jacobian_cubic_coefficients'] if i==0 else data['cubic']))
        E=[QQ(5),16*ai[3](tau)+8,64*ai[4](tau)+16]
        models=[J,E]
        fixed={'label':label,'tau':str(tau),'source':str(path.relative_to(ROOT)),
          'carrier_abc':list(map(str,J)),'original_abc':list(map(str,E))}
        retain(OUT/f'case-{i:02d}-equations.json',fixed)
        witnesses={'carrier_irreducible':None,'carrier_transposition':None,
          'original_irreducible':None,'original_transposition':None,
          'distinct_quadratic_fields':None,'geometric_nonisogeny':None}
        trials=[]
        for p in protocol['primes']:
            trial={'p':p};reductions=[]
            if any(v.denominator()%p==0 for model in models for v in model):
                trial['status']='SKIP_MODEL_DENOMINATOR';trials.append(trial);continue
            for model in models:
                abc=list(map(lambda x:residue(x,p),model));a,b,c=abc
                reductions.append(abc)
            discs=[(a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c)%p for a,b,c in reductions]
            trial['abc']=reductions;trial['discriminants_mod_p']=discs
            if 0 in discs:
                trial['status']='SKIP_BAD_REDUCTION';trials.append(trial);continue
            trial['status']='COMPLETE_GOOD_PAIR';types=[];traces=[];counts=[]
            P=PolynomialRing(GF(p),'x');x=P.gen()
            for abc in reductions:
                a,b,c=abc;f=x**3+a*x*x+b*x+c
                types.append(sorted(int(g.degree()) for g,e in f.factor() for _ in range(e)))
                n=int(EllipticCurve(GF(p),[0,a,0,b,c]).cardinality())
                counts.append(n);traces.append(p+1-n)
            trial['factor_types']=types;trial['point_counts']=counts;trial['traces']=traces
            for prefix,typ in zip(['carrier','original'],types):
                key=prefix+('_irreducible' if typ==[3] else '_transposition' if typ==[1,2] else '_other')
                if key in witnesses and witnesses[key] is None:witnesses[key]=p
            signs=[int(GF(p)(d).is_square()) for d in discs]
            trial['discriminant_is_square']=signs
            if signs[0]!=signs[1] and witnesses['distinct_quadratic_fields'] is None:
                witnesses['distinct_quadratic_fields']=p
            ds=[a*a-4*p for a in traces];ordinary=all(a%p!=0 for a in traces)
            trial['ordinary_pair']=ordinary;trial['Frobenius_discriminants']=ds
            trial['distinct_Frobenius_fields']=ordinary and not ZZ(ds[0]*ds[1]).is_square()
            if trial['distinct_Frobenius_fields'] and witnesses['geometric_nonisogeny'] is None:
                witnesses['geometric_nonisogeny']=p
            trials.append(trial)
        passed=all(v is not None for v in witnesses.values())
        result={'label':label,'tau':str(tau),'status':'CERTIFIED_TRANSPORT_OBSTRUCTION' if passed else 'UNKNOWN_IN_FIXED_FOOTPRINT',
          'witnesses':witnesses,'trials':trials,
          'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/f'case-{i:02d}-equations.json']}}
        retain(OUT/f'case-{i:02d}.json',result);results.append({'label':label,'status':result['status'],'witnesses':witnesses})
        print('CHECKPOINT',i,label,result['status'],witnesses,flush=True)
    retain(OUT/'panel.json',{'classification':'verified application of ordinary-reduction and joint2-torsion obstructions',
      'cases':results,'all_certified':all(r['status']=='CERTIFIED_TRANSPORT_OBSTRUCTION' for r in results),
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',*[OUT/f'case-{i:02d}.json' for i in range(10)]]}})
if __name__=='__main__':
    signal.alarm(25);OUT.mkdir(parents=True,exist_ok=True);main()
