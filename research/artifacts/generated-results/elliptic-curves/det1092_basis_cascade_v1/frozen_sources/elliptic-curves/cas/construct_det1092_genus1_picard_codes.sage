#!/usr/bin/env sage-python
"""Freeze full inherited Picard-image finite codes before marked evaluation.

All18 nontrivial candidate generators (F and17 old sections), both fixed
curves, all64 old primes. Incomplete reductions are skipped as whole places,
never by dropping an inherited generator.25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2]
DIR=ROOT/'artifacts/generated-results/elliptic-curves/det1092_genus1_picard_image_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def main():
    protocol=read(DIR/'protocol.json');all_cases=[]
    for index in range(2):
        path=DIR/f'case-{index:02d}-generic.json';data=read(path)
        abc=list(map(QQ,data['Jacobian_cubic_coefficients']))
        J=EllipticCurve(QQ,[0,abc[0],0,abc[1],abc[2]])
        HF=list(map(QQ,data['old_F_restriction']['point']));J(HF)
        R=PolynomialRing(QQ,'t')
        divs=[]
        for d in data['old_section_restrictions']:
            divs.append(None if d['status'].startswith('ZERO_') else
                        (R(d['g']),R(d['Jacobian_Xmod'])))
        rows=[];owners=[];trials=[]
        for p in protocol['primes']:
            trial={'p':p}
            if any(v.denominator()%p==0 for v in abc):
                trial['status']='SKIP_JACOBIAN_DENOMINATOR';trials.append(trial);continue
            if J.discriminant().valuation(p)!=0:
                trial['status']='SKIP_BAD_JACOBIAN_REDUCTION';trials.append(trial);continue
            a,b,c=[residue(v,p) for v in abc]
            roots=[x for x in range(p) if (x**3+a*x*x+b*x+c)%p==0]
            trial.update(cubic=[c,b,a,1],rational_roots=roots)
            if not roots:
                trial.update(status='COMPLETE_ZERO_LOCAL_QUOTIENT',rows=[],norms=[])
                trials.append(trial);continue
            if any(v.denominator()%p==0 for pair in divs if pair for f in pair for v in f):
                trial['status']='SKIP_FULL_DIVISOR_DENOMINATOR';trials.append(trial);continue
            S=PolynomialRing(GF(p),'t');raw=[];norms=[];failed=False
            for root in roots:
                values=[]
                if HF[0].valuation(p)<0:values.append(1)
                else:
                    x=residue(HF[0],p);value=(x-root)%p
                    if not value:value=(3*root*root+2*a*root+b)%p
                    assert value;values.append(value)
                for pair in divs:
                    if pair is None:values.append(1);continue
                    gp,Xp=map(S,pair)
                    norm=int(gp.resultant(Xp-root));values.append(norm)
                    if not norm:failed=True
                norms.append(values)
                if failed:break
                raw.append([int(pow(v,(p-1)//2,p)==p-1) for v in values])
            if failed:
                trial.update(status='SKIP_NONUNIT_DIVISOR_CLASS',partial_norms=norms)
                trials.append(trial);continue
            assert len(raw)==len(roots)
            if len(roots)==3:assert all(sum(r[j] for r in raw)%2==0 for j in range(18))
            trial.update(status='COMPLETE_FULL_PICARD_REDUCTION',norms=norms,rows=raw)
            rows+=raw;owners += [p]*len(raw);trials.append(trial)
        M=matrix(GF(2),len(rows),18,[v for r in rows for v in r]);rank=int(M.rank())
        assert rank<=data['rank_upper_bound']
        result={'classification':'full inherited finite image, not a Selmer group',
          'index':index,'trials':trials,'matrix_rows':rows,'row_primes':owners,
          'rank':rank,'rank_upper_bound':data['rank_upper_bound'],
          'status':'FULL_RANK_PICARD_IMAGE' if rank==data['rank_upper_bound'] else 'UNRESOLVED_IMAGE_RANK',
          'marked_point_evaluated':False,
          'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [path,DIR/'protocol.json',Path(__file__)]}}
        retain(DIR/f'case-{index:02d}-codes.json',result)
        all_cases.append({'index':index,'rank':rank,'upper':data['rank_upper_bound'],
                          'complete_primes':sum(r['status'].startswith('COMPLETE_') for r in trials)})
        print('CHECKPOINT_FROZEN_FULL_PICARD_CODE',all_cases[-1],flush=True)
    retain(DIR/'codes-frozen.json',{'status':'GENERIC_PICARD_CODES_FROZEN_BEFORE_MARKED_POINT',
       'cases':all_cases,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
         [DIR/'case-00-codes.json',DIR/'case-01-codes.json',DIR/'protocol.json',Path(__file__)]}})
if __name__=='__main__':
    signal.alarm(25);main()
