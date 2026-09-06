#!/usr/bin/env python3
"""Portable arithmetic replay of the labelled native-carrier Sha obstruction."""
import argparse
from fractions import Fraction as F
from math import gcd, isqrt, prod
from pathlib import Path
import re
import retrospective as r
from cover_experiment import evaluate, mul, sub, trim, sqrtq
from carrier_ct_local_witness import hilbert, valuation_unit

HERE=Path(__file__).resolve().parent
LABEL=r.OUT/'rank_jump_labelled_carrier_ct_v2.json'
LABEL_INPUT=r.OUT/'rank_jump_labelled_carrier_ct_inputs_v2.json'
LOCAL=r.OUT/'rank_jump_carrier_ct_local_witness_v1.json'
PRIMES=r.OUT/'rank_jump_carrier_ct_primality_v1.json'
OUTPUT=r.OUT/'rank_jump_labelled_carrier_ct_verification_v1.json'


def capture_primes():
    from sage.all import ZZ
    proof={}
    def certify(p):
        if str(p) in proof:return
        if p==2:proof['2']={'prime':2};return
        factors=[(int(q),int(e)) for q,e in ZZ(p-1).factor(proof=True)]
        for q,e in factors:certify(q)
        a=next(a for a in range(2,1001) if pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q,e in factors))
        proof[str(p)]={'prime':p,'predecessor_factorization':factors,'witness':a}
    for row in r.read(LOCAL)['local_evaluations']:certify(row['prime'])
    r.write_new(PRIMES,{'schema':'rank-jump.carrier-ct-primality.v1','proofs':proof,
                        'local_result_sha256':r.digest(LOCAL.read_bytes())})


def gp_polynomial(text):
    text=text.replace(' ','');assert re.fullmatch(r'[0-9x^*/+\-]+',text)
    coeff={}
    for term in text.replace('-','+-').split('+'):
        if not term:continue
        if 'x' not in term:k,c=0,F(term)
        else:
            left,right=term.split('x');k=int(right[1:]) if right else 1
            left=left.rstrip('*');c=F(-1 if left=='-' else 1 if left=='' else left)
        coeff[k]=coeff.get(k,F(0))+c
    return trim([coeff.get(k,F(0)) for k in range(max(coeff)+1)])


def determinant(rows):
    a=[list(map(F,row)) for row in rows];ans=F(1)
    for k in range(len(a)):
        i=next((i for i in range(k,len(a)) if a[i][k]),None)
        if i is None:return F(0)
        if i!=k:a[i],a[k]=a[k],a[i];ans=-ans
        pivot=a[k][k];ans*=pivot
        for i in range(k+1,len(a)):
            c=a[i][k]/pivot
            a[i]=[x-c*y for x,y in zip(a[i],a[k])]
    return ans


def resultant(f,g):
    m,n=len(f)-1,len(g)-1
    if n==0:return g[0]**m
    rows=[]
    for i in range(n):rows.append([0]*i+list(reversed(f))+[0]*(n-1-i))
    for i in range(m):rows.append([0]*i+list(reversed(g))+[0]*(m-1-i))
    return determinant(rows)


def remainder(f,T):
    f=trim(f);assert T[-1]==1
    while len(f)>=len(T):
        shift=len(f)-len(T);c=f[-1]
        f=sub(f,[0]*shift+[c*x for x in T])
    return f+[F(0)]*(3-len(f))


def productK(a,b,T):return remainder(mul(a,b),T)


def inverseK(a,T):
    columns=[productK(a,[F(int(i==j)) for i in range(3)],T) for j in range(3)]
    M=[[columns[j][i] for j in range(3)]+[F(int(i==0))] for i in range(3)]
    for k in range(3):
        i=next(i for i in range(k,3) if M[i][k]);M[i],M[k]=M[k],M[i]
        c=M[k][k];M[k]=[x/c for x in M[k]]
        for i in range(3):
            if i!=k:
                c=M[i][k];M[i]=[x-c*y for x,y in zip(M[i],M[k])]
    return [M[i][-1] for i in range(3)]


def fingerprints(polys,cubic,roster):
    sigs=[0]*len(polys);col=0;used=[]
    for p in roster:
        if any(c.denominator%p==0 for f in polys+[cubic] for c in f):continue
        coeff=[r.mod(c,p) for c in cubic]
        roots=[x for x in range(p) if sum(c*pow(x,i,p) for i,c in enumerate(coeff))%p==0]
        if len(roots)!=3:continue
        vals=[[r.mod(evaluate(f,F(x)),p) for x in roots] for f in polys]
        if any(not x for row in vals for x in row):continue
        used.append(p)
        for j,row in enumerate(vals):
            for k,x in enumerate(row):sigs[j]|=int(pow(x,(p-1)//2,p)==p-1)<<(col+k)
        col+=3
    return sigs,used


def coordinates(sigs,target):
    answers=[]
    for mask in range(1<<len(sigs)):
        v=0
        for i,s in enumerate(sigs):
            if mask>>i&1:v^=s
        if v==target:answers.append([(mask>>i)&1 for i in range(len(sigs))])
    assert len(answers)==1
    return answers[0]


def verify():
    label=r.read(LABEL);inp=r.read(LABEL_INPUT);loc=r.read(LOCAL);primes=r.read(PRIMES)
    for data in (label,loc):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    done=set()
    def prime(p):
        if p in done:return
        cert=primes['proofs'][str(p)];assert cert['prime']==p
        if p!=2:
            factors=cert['predecessor_factorization'];a=cert['witness']
            assert len({q for q,e in factors})==len(factors) and all(e>0 for q,e in factors)
            assert prod(q**e for q,e in factors)==p-1
            for q,e in factors:prime(q)
            assert pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q,e in factors)
        done.add(p)
    for row in loc['local_evaluations']:prime(row['prime'])
    assert primes['local_result_sha256']==r.digest(LOCAL.read_bytes())
    cubic=list(map(F,label['cubic_coefficients']));basis=[list(map(F,f)) for f in label['Selmer_basis_coefficients']]
    beta=list(map(F,label['carrier_beta_coefficients']))
    assert cubic==gp_polynomial(inp['internal_cubic_GP'])
    assert basis==[gp_polynomial(t) for t in inp['basis_GP']]
    norms=[resultant(cubic,f) for f in basis+[beta]]
    assert all(n and sqrtq(n) is not None for n in norms)
    model=label['internal_model'];assert cubic==[F(model[4]),F(model[3]),F(model[1]),F(1)]
    old=r.read(r.OUT/'rank_jump_carrier_sha_class_v1.json')['mapping'];lam=F(label['raw_to_internal_short_scaling'])
    short,_=r.short(model,[])
    assert F(old['Jacobian_model'][3])==lam**4*F(short[3]) and F(old['Jacobian_model'][4])==lam**6*F(short[4])
    b0,b1=map(F,old['beta']);assert beta==[b0+b1*lam*lam*F(model[1])/3,b1*lam*lam]
    sigs,used=fingerprints(basis+[beta],cubic,label['proof_primes'])
    assert used==label['proof_primes'] and all(p in r.primes(503) for p in used)
    assert sigs[:5]==label['Selmer_fingerprints'] and sigs[5]==label['carrier_fingerprint']
    assert r.rank(sigs[:5])==r.rank(sigs)==5
    coeff=coordinates(sigs[:5],sigs[5]);assert coeff==label['carrier_coordinates']==[0,0,1,0,0]
    C=label['Cassels_matrix'];raw=[[int(x.strip()) for x in row.split(',')] for row in inp['Cassels_matrix_GP'].strip('[]').split(';')]
    assert C==raw and len(C)==5 and all(len(row)==5 for row in C)
    assert all(C[i][j]==C[j][i] for i in range(5) for j in range(5)) and all(C[i][i]==0 for i in range(5))
    assert r.rank([r.pack(row) for row in C])==2
    pairing=[sum(row[j]*coeff[j] for j in range(5))%2 for row in C]
    assert pairing==label['carrier_pairing_row']==[1,1,0,0,0]
    # Identify the actual quartics used in the local formula in the same Selmer basis.
    qs=[list(map(F,q)) for q in loc['quartics']];I,J=map(F,loc['invariants_I_J'])
    zclasses=[];scale2=sqrtq(-27*I/F(short[3]));assert scale2 is not None and sqrtq(scale2) is not None
    assert -27*J==scale2**3*F(short[4])
    for q in qs:
        e,d,c,b,a=q
        assert I==12*a*e-3*b*d+c*c and J==72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
        zclasses.append([b*b-8*a*c/3-4*a*scale2*F(model[1])/27,-4*a*scale2/9])
    ss,rr=fingerprints(basis+zclasses,cubic,r.primes(503));assert r.rank(ss[:5])==5
    qc=[coordinates(ss[:5],s) for s in ss[5:]]
    assert qc==[[0,0,1,0,0],[1,0,0,0,0],[1,0,1,0,0]]
    # Fisher's cubic square and auxiliary-quadratic calculation, without Sage.
    T=[J,-3*I,F(0),F(1)];m=list(map(F,loc['cubic_square_root']))
    zs=[[q[3]**2-8*q[4]*q[2]/3,4*q[4]/3,F(0)] for q in qs]
    assert productK(m,m,T)==productK(productK(zs[0],zs[1],T),zs[2],T)
    multiplier=productK(m,inverseK(zs[0],T),T)
    e,d,c,b,a=qs[0]
    H=[[8*c*c-12*(b*d+8*a*e)+8*I,8*c,F(-8)],
       [24*(b*c-6*a*d),24*b,F(0)],
       [12*(3*b*b-8*a*c),48*a,F(0)]]
    gamma_raw=[productK(multiplier,h,T)[2] for h in H]
    gamma=list(map(F,loc['auxiliary_quadratic']));ratio=F(loc['auxiliary_quadratic_scalar'])
    assert gamma_raw==[ratio*x for x in gamma] and ratio
    aa=F(loc['a']);assert aa==qs[1][4] and aa>0
    support=[row['prime'] for row in loc['local_evaluations']]
    required=abs(int(aa*(4*I**3-J**2)));assert required
    for p in support:
        while required%p==0:required//=p
    assert required==1 and {2,3,5,7}<=set(support)
    product=1;negative=[]
    for row in loc['local_evaluations']:
        p=row['prime'];x=F(row['local_coordinate']);gv=evaluate(gamma,x);qv=evaluate(qs[0],x)
        assert str(gv)==row['gamma_value'] and str(qv)==row['quartic_value'] and gv and qv
        v,u=valuation_unit(qv,p)
        assert v%2==0 and (r.mod(u,8)==1 if p==2 else pow(r.mod(u,p),(p-1)//2,p)==1)
        symbol=hilbert(aa,gv,p);assert symbol==row['Hilbert_symbol'];product*=symbol
        if symbol==-1:negative.append(p)
    assert product==-1 and negative==[2] and loc['real_Hilbert_symbol']==1
    assert label['carrier_global_solubility']=='NO' and label['carrier_nonzero_Sha_class']
    return {'schema':'rank-jump.labelled-carrier-ct-verification.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (LABEL,LABEL_INPUT,LOCAL,PRIMES,Path(__file__),HERE/'retrospective.py',HERE/'cover_experiment.py',HERE/'carrier_ct_local_witness.py')},
            'Selmer_norm_square_roots':list(map(str,(sqrtq(n) for n in norms))),
            'carrier_coordinates':coeff,'nonzero_pairing_row':pairing,'local_formula_quartic_coordinates':qc,
            'Hilbert_product':product,'negative_symbol_primes':negative,'primality_certificates_checked':len(done),
            'carrier_global_solubility':'NO',
            'dependency_boundary':'The full Selmer basis and completeness of Fisher place support inherit the pinned PARI computation. Class coordinates, model transport, norms, cubic formula, local evaluations, Hilbert product and primality are replayed arithmetically without that runtime.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['primes','build','check']);a=p.parse_args()
    if a.mode=='primes':capture_primes()
    else:
        result=verify()
        if a.mode=='build':r.write_new(OUTPUT,result)
        else:assert r.read(OUTPUT)==result
        print('PASS native carrier Sha obstruction; Hilbert product',result['Hilbert_product'])
