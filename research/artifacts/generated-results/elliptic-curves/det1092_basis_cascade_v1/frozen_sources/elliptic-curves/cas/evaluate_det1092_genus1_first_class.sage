#!/usr/bin/env sage-python
"""Evaluate the reconstructed first class only after full generic-code freeze.

No saved point coordinates, later points or point search. Calibration is in
the previously fixed curve equation and is not relabelled prospective.25s.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
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
    paths=[DIR/'codes-frozen.json',DIR/'case-01-generic.json',DIR/'case-01-codes.json',Path(__file__)]
    frozen,data,codes=map(read,paths[:3])
    for d in [frozen,data,codes]:
        for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest
    R=PolynomialRing(QQ,'t');t=R.gen();q=R(data['quartic'])
    t0,s=map(QQ,data['basepoint']);sq=q(0).sqrt();u=-t0
    q0,q1,q2,q3,q4=q(t+t0).list();assert q0==s*s and u
    X=(2*s*(sq+s)+q1*u)/(u*u)
    Y=((X*X-4*s*s*q4)*u-q1*X-2*s*s*q3)/(2*s)
    a,b,c=map(QQ,data['Jacobian_cubic_coefficients'])
    J=EllipticCurve(QQ,[0,a,0,b,c]);P=J([X,Y])
    column=[];evaluations=[]
    for trial in codes['trials']:
        if not trial['status'].startswith('COMPLETE_'):continue
        p=trial['p'];cc,bb,aa,one=trial['cubic'];assert one==1
        bits=[]
        for root in trial['rational_roots']:
            if X.valuation(p)<0:value=1
            else:
                value=(residue(X,p)-root)%p
                if not value:value=(3*root*root+2*aa*root+bb)%p
            assert value;bits.append(int(pow(value,(p-1)//2,p)==p-1))
        column+=bits;evaluations.append({'p':p,'bits':bits})
    M=matrix(GF(2),codes['matrix_rows']);Pcol=vector(GF(2),column)
    A=M.augment(matrix(GF(2),len(column),1,column));rank=int(A.rank())
    separator=None
    if rank>M.rank():
        separator=next(v for v in M.left_kernel().basis() if v*Pcol)
        assert separator*M==0 and separator*Pcol==1
    result={'classification':'retrospective marked rational class on the elliptic carrier',
      'status':'FIRST_CLASS_OUTSIDE_FULL_PICARD_IMAGE' if separator is not None and M.rank()==data['rank_upper_bound'] else 'UNRESOLVED_FIRST_CLASS',
      'quartic_point':['0',str(sq)],'Jacobian_point':list(map(str,P[:2])),
      'generic_rank':int(M.rank()),'generic_rank_upper_bound':data['rank_upper_bound'],
      'augmented_rank':rank,'prime_evaluations':evaluations,
      'separating_row':list(map(int,separator)) if separator is not None else None,
      'boundary':'Nonzero separation is against all inherited Picard generators, not just one base point. Rank equality alone would not prove global membership. This calibrated-carrier class is not yet an oracle-free selector or a Selmer/Sha computation.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'first-class.json',result)
    print(result['status'],'full image',M.rank(),'augmented',rank,flush=True)
if __name__=='__main__':
    signal.alarm(25);main()
