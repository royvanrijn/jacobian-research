#!/usr/bin/env sage-python
"""Two finite-field point counts bound Picard rank of an existing K3 parent."""
import argparse,sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
import numpy as np
from sage.all import QQ,ZZ,AA,GF,PolynomialRing,EllipticCurve,prime_range,cyclotomic_polynomial
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
surface=SourceFileLoader('mestre_surface_count_source',str(CAS/'audit_mestre_parent_surfaces.sage')).load_module()
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/mestre-picard-count-v1'

def legendre_table(p):
    chars=np.full(p,-1,dtype=np.int64);chars[0]=0
    for a in range(1,p):chars[a*a%p]=1
    return chars

def count(A,B,p,extension,progress):
    aa,bb,delta,_=surface.good_model(A,B,p);chars=legendre_table(p)
    d=next(n for n in range(2,p) if chars[n]==-1)
    if extension==1:
        F=GF(p);q=p;x0=np.arange(p,dtype=np.int64);x1=None;z=None
    else:
        R=PolynomialRing(GF(p),'z');F=GF(p*p,'z',modulus=R.gen()**2-d);z=F.gen();q=p*p
        index=np.arange(q,dtype=np.int64);x0=index%p;x1=index//p
        cube0=(x0*x0*x0+3*d*x0*x1*x1)%p;cube1=(3*x0*x0*x1+d*x1*x1*x1)%p
    R=PolynomialRing(F,'T');a=R(aa);b=R(bb);disc=-16*(4*a**3+27*b**2)
    def pair_eval(coeffs,t0,t1):
        r0=r1=0
        for c in reversed(coeffs):r0,r1=(r0*t0+d*r1*t1+int(c))%p,(r0*t1+r1*t0)%p
        return r0,r1
    counts=[];repairs=[];smooth=0;count_total=0
    for index in range(q+1):
        if index==q:
            t=None;a0,a1=int(aa[8]),0;b0,b1=int(bb[12]),0;av,bv=F(a0),F(b0)
        else:
            t0=index%p;t1=0 if extension==1 else index//p
            t=F(t0) if extension==1 else F(t0)+F(t1)*z
            a0,a1=pair_eval(aa.list(),t0,t1);b0,b1=pair_eval(bb.list(),t0,t1)
            av,bv=a(t),b(t)
            if extension==1:
                if av!=F(a0) or bv!=F(b0) or a1 or b1:raise ArithmeticError('prime-field coefficient evaluation differs')
            elif av!=F(a0)+F(a1)*z or bv!=F(b0)+F(b1)*z:raise ArithmeticError('quadratic-field coefficient evaluation differs')
        if extension==1:
            values=(x0*x0*x0+a0*x0+b0)%p;naive=q+1+int(chars[values].sum())
        else:
            f0=(cube0+a0*x0+d*a1*x1+b0)%p;f1=(cube1+a0*x1+a1*x0+b1)%p
            norm=(f0*f0-d*f1*f1)%p;naive=q+1+int(chars[norm].sum())
        if 4*av**3+27*bv**2:
            independent=int(EllipticCurve(F,[av,bv]).cardinality(algorithm='pari'))
            if independent!=naive:raise ArithmeticError('independent PARI cardinality differs from norm-character count')
            smooth+=1;correction=0
        elif index==q:correction=3*q
        else:
            multiplicity=0;temp=disc;T=R.gen()
            while temp(t)==0:temp=temp//(T-t);multiplicity+=1
            if multiplicity not in (1,2):raise ArithmeticError('semistable reduction type changed')
            correction=(multiplicity-1)*q
        if correction:repairs.append({'base_index':index,'correction':correction})
        counts.append(naive);count_total+=naive+correction
        if (index+1)%256==0:progress(index+1,q+1)
    if len(repairs)!=3 or sorted(r['correction'] for r in repairs)!=[q,q,3*q]:raise ArithmeticError('two I2 and one split I4 correction required')
    return {'prime':p,'extension_degree':extension,'field_order':q,'quadratic_nonsquare':None if extension==1 else d,
      'base_order':'index=a+p*b represents a+b*z; infinity has index q',
      'weierstrass_fibre_counts':counts,'resolution_corrections':repairs,'surface_point_count':count_total,
      'smooth_fibres_independently_counted':smooth,'all_smooth_counts_verified_by':'PARI via Sage elliptic cardinality'}

def picard_bound(p,n1,n2):
    R=PolynomialRing(QQ,'z');z=R.gen();a=ZZ(n1-1-p*p-18*p);s2=ZZ(n2-1-p**4-18*p*p);b=QQ(a*a-s2)/2
    if b.denominator()!=1:raise ArithmeticError('integral residual Frobenius coefficient required')
    candidates=[]
    for sign in (1,-1):
        if sign==-1 and b:continue
        f=z**4-a*z**3+b*z*z-sign*p*p*a*z+sign*p**4
        if sign==1:
            g=z*z-QQ(a)/p*z+QQ(b)/(p*p)-2;roots=g.roots(AA)
            if sum(m for _,m in roots)!=2 or any(abs(x)>2 for x,_ in roots):continue
        elif abs(a)>2*p:continue
        normalized=R(f(p*z)/p**4);remainder=normalized;cycles=[]
        for n in (1,2,3,4,5,6,8,10,12):
            phi=R(cyclotomic_polynomial(n));multiplicity=0
            while remainder and remainder%phi==0:remainder=remainder//phi;multiplicity+=1
            if multiplicity:cycles.append({'order':n,'multiplicity':multiplicity,'degree':int(phi.degree())})
        cycle_degree=sum(r['multiplicity']*r['degree'] for r in cycles)
        candidates.append({'orthogonal_determinant_sign':sign,'residual_characteristic_polynomial_coefficients':list(map(str,f.list())),
          'cyclotomic_factors_after_dividing_eigenvalues_by_p':cycles,'geometric_picard_upper_bound':18+cycle_degree})
    if not candidates:raise ArithmeticError('no weight-two reciprocal residual polynomial matches counts')
    upper=max(r['geometric_picard_upper_bound'] for r in candidates)
    return {'known_rational_divisor_rank':18,'residual_trace_1':int(a),'residual_trace_2':int(s2),'residual_second_coefficient':str(b),
      'frobenius_candidates':candidates,'geometric_picard_upper_bound':upper,'geometric_picard_rank':18 if upper==18 else None,
      'generic_K3_MW_rank_over_Qbar_T':11 if upper==18 else None,
      'scope':'Smooth proper good reduction injects the characteristic-zero NS group. The18 rational classes contribute18 copies of eigenvalue p. Orthogonal reciprocity leaves a degree4 factor determined by two traces unless its middle coefficient is zero; all remaining Weil-compatible determinant signs are retained. Only p times roots of unity can represent further geometric divisor classes. This gives an upper bound without assuming the Tate conjecture.'}

def worker(index):
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen sources or proofs changed')
    row=p['rows'][index];out=D/row['id']/'counts.json'
    if out.exists():raise FileExistsError('preserve finite-field count')
    parent=next(r for r in cert.read(ART/'mestre_parent_portfolio_intake_v1.json')['rows'] if r['outer_u']==row['outer_u'])
    A,B=surface.model(parent);prime=row['prime']
    if surface.good_model(A,B,prime) is None:raise ArithmeticError('common smooth-reduction gate failed')
    data={'status':'RUNNING','protocol_hash':digest(p),'outer_u':row['outer_u'],'prime':prime,'counts':[]};checkpoint(out,data)
    for extension in (1,2):
        c=count(A,B,prime,extension,lambda i,n:print(row['id'],extension,i,'of',n,flush=True));data['counts'].append(c);checkpoint(out,data)
    data.update(status='PASS',picard=picard_bound(prime,*[r['surface_point_count'] for r in data['counts']]));checkpoint(out,data)
    print('PASS',row['id'],'PICARD UPPER',data['picard']['geometric_picard_upper_bound'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);a=p.parse_args();worker(a.index)
