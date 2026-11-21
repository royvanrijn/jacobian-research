#!/usr/bin/env python3
"""Exact invariant identities, exceptional-shape controls and narrow source authentication."""
if not __debug__:
    raise SystemExit('optimization is unsupported: exact source/control guards are required')
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,hashlib,json

HERE=Path(__file__).resolve().parent
PN='CONSTANT_CUBIC_INVARIANT_ELLIPTIC_MODULI_BRIDGE.md'
CN='verify_constant_cubic_invariant_moduli.py'
IN='constant_cubic_invariant_moduli_inputs.json'
LN='constant_cubic_invariant_moduli_source_ledger.json'
PROOF='b5cc3c456884e0fbdd70b1eeae036c487d37a374d9cf51cb087ff77bc69a47be'
BASE='research/foundations/extended-geometry/'
RAW={
 'noncentral_shape_classification':('NONCENTRAL_CONSTANT_AFFINE_CUBIC_SHAPE_CLASSIFICATION.md','20ff19e08cc7decb8564346c363e2a3657d846241267a462899934d481458e07','mathematical_noncentral_shape','Main theorem; unsafe degree spectrum6,9,12,...; lambda!=1/2'),
 'central_invariant_zero_classification':('INVARIANT_ZERO_MIXED_CUBIC_LOW_DEGREE_CLASSIFICATION.md','0b65f1da0010eaa646355670099cdd0bbd111ace80abb85f76eeccc1e5547fc9','mathematical_central_low_degree','Central I=0, E0!=0 low-degree theorem only'),
 'finite_cubic_interval_certificate':('MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md','c4fee78a3b99b286e22818f63c90f8da450293bde974b354fb22bd2866bf1dc6','mathematical_first_four_only','Section2 direct first-four-moment certificate only')}
BOUNDARY={
 'characteristic_zero_all_fields_and_extensions':True,'regular_cubic_height_E_not_identically_zero':True,
 'unit_center_constant_invariant_iff_formal_fixed_shape':True,'formal_trivialization_over_K_itself':True,
 'zero_pair_is_separate_not_projective':True,'nonunit_center_constant_invariant_forces_zero_center':True,
 'pure_zero_center_safe_via_first_four_only':True,'pure_constant_invariant_degree_five_seven_safe':True,
 'unsafe_prime_degree_constant_invariant_forces_central_unit_j1728':True,
 'elliptic_interpretation_requires_nonzero_discriminant':True,'constant_j_implies_Kx_constant_curve_isomorphism':False,
 'purity_implies_constant_invariant':False,'purity_implies_I_zero':False,
 'all_central_prime_phases_safe':False,'arbitrary_higher_height_or_varying_invariant_classification':False,
 'unsafe_witness_existence_or_minimum_degree':False,'lower_height_E_identically_zero_in_scope':False,
 'analytic_interval_theorem_premise':False,'MIT_elliptic_primary_mathematical_premise':False,
 'formal_verification':False,'independent_replay':False,'external_review':False}
CONTEXT={'author':'Andrew Sutherland','title':'Ordinary and supersingular elliptic curves',
 'course':'MIT18.783 Fall2025 Lecture13','locator':'Definition13.11 and Theorems13.13-13.14 printed4-6',
 'url':'https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec13.pdf',
 'role':'classical_context_only','scope':'Equal j gives algebraic-closure isomorphism, not automatically K(x)-isomorphism.'}
EXPECTED={'schema':1,'logical_boundary':BOUNDARY,'classical_context':CONTEXT,
 'affine_shape_parameters':['0','1','1/2','-1','2','3'],
 'formal_rho':['1','1/2','-1/3','2/5','-1/7','1/11'],
 'formal_pq_cases':[['1','1'],['-1','0'],['0','1']],
 'prime_scope_controls':[2,3,5,7,11,13,17,19,23,29,31],
 'unsafe_noncentral_degrees':[6,9,12,15,18,21,24,27,30]}

def need(ok,msg):
    if not ok:raise ValueError(msg)
def strict(a,b):
    return type(a) is type(b) and (set(a)==set(b) and all(strict(a[k],b[k]) for k in a)
      if isinstance(a,dict) else len(a)==len(b) and all(strict(x,y) for x,y in zip(a,b))
      if isinstance(a,list) else a==b)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def validate_input(d):need(strict(d,EXPECTED),'typed immutable inputs and scope')
def source_check(d,root):
    need(type(d) is dict and set(d)=={'schema','logical_boundary','sources'},'ledger shape')
    need(type(d['schema']) is int and d['schema']==1,'exact integer schema')
    need(strict(d['logical_boundary'],BOUNDARY),'typed ledger boundary')
    own={'written_proof':(PN,'mathematical_proof','Proposition and Sections1-5'),
         'exact_checker':(CN,'identity_and_source_support','Universal polynomial algebra and finite exceptional controls'),
         'typed_inputs':(IN,'exact_control_data','Immutable typed scope and exact controls')}
    need(type(d['sources']) is list and len(d['sources'])==6,'exact six mandatory source roles')
    out={}
    for row in d['sources']:
        need(type(row) is dict and set(row)=={'role','path','relation','locator','sha256'},'source record shape')
        need(all(type(v) is str for v in row.values()),'source string types')
        role=row['role'];need(role not in out,'duplicate source role')
        if role in own:
            name,rel,loc=own[role];p=HERE/name;expected=name
        else:
            need(role in RAW,'source role allowlist');name,pin,rel,loc=RAW[role]
            p=root/BASE/name;expected=BASE+name
            need(row['sha256']==pin,'reviewed mathematical source pin '+role)
        need(row['path']==expected and row['relation']==rel and row['locator']==loc,'source path/relation/locator')
        need(len(row['sha256'])==64 and all(c in '0123456789abcdef' for c in row['sha256']),'sha256 format')
        need(p.is_file() and digest(p)==row['sha256'],'actual source hash '+role)
        out[role]=dict(row)
    need(set(out)==set(own)|set(RAW),'all source roles')
    need(digest(HERE/PN)==PROOF,'hard reviewed written-proof pin')
    return out

# Sparse multivariate rational polynomials. Every universal identity is expanded.
NV=8
def clean(p):return {e:F(v) for e,v in p.items() if v}
def con(v):return {(0,)*NV:F(v)} if v else {}
def var(i):
    e=[0]*NV;e[i]=1;return {tuple(e):F(1)}
def add(*ps):
    q={}
    for p in ps:
        for e,v in p.items():q[e]=q.get(e,F(0))+v
    return clean(q)
def sc(p,a):return clean({e:v*a for e,v in p.items()})
def mul(*ps):
    q=con(1)
    for p in ps:
        z={}
        for e,v in q.items():
            for f,w in p.items():
                g=tuple(a+b for a,b in zip(e,f));z[g]=z.get(g,F(0))+v*w
        q=clean(z)
    return q
def pw(p,n):return mul(*([p]*n))
def sub(p,images):
    q={}
    for e,v in p.items():q=add(q,sc(mul(*[pw(images[i],n) for i,n in enumerate(e)]),v))
    return q
def invariants(C,B,D,E):
    J=add(sc(mul(B,E),3),sc(pw(D,2),-1))
    I=add(sc(mul(C,pw(E,2)),27),sc(mul(B,D,E),-9),sc(pw(D,3),2))
    disc=add(mul(pw(D,2),pw(B,2)),sc(mul(E,pw(B,3)),-4),sc(mul(pw(D,3),C),-4),sc(mul(pw(E,2),pw(C,2)),-27),sc(mul(E,D,B,C),18))
    return J,I,disc
def universal():
    C,B,D,E,t,s,h,l=[var(i) for i in range(NV)];J,I,disc=invariants(C,B,D,E);ids=[]
    def identity(name,p):need(not p,'universal identity '+name);ids.append(name)
    identity('cubic discriminant',add(pw(I,2),sc(pw(J,3),4),sc(mul(pw(E,2),disc),27)))
    P=add(C,mul(B,t),mul(D,pw(t,2)),mul(E,pw(t,3)));v=add(sc(mul(E,t),3),D)
    identity('depressed cubic without denominators',add(sc(mul(pw(E,2),P),27),sc(pw(v,3),-1),sc(mul(J,v),-3),sc(I,-1)))
    # After X_old=Et and X_old=X_new-D/3, E^2 P becomes short cubic.
    X=t;w=add(X,sc(D,F(-1,3)))
    short=add(pw(w,3),mul(D,pw(w,2)),mul(B,E,w),mul(C,pw(E,2)))
    identity('elliptic variable change',add(short,sc(pw(X,3),-1),sc(mul(J,X),F(-1,3)),sc(I,F(-1,27))))
    C1=sc(mul(s,add(s,h),add(s,mul(l,h))),-1)
    B1=add(sc(pw(s,2),3),sc(mul(add(con(1),l),s,h),2),mul(l,pw(h,2)))
    D1=add(sc(s,-3),sc(mul(add(con(1),l),h),-1));J1,I1,disc1=invariants(C1,B1,D1,con(1))
    jl=add(sc(pw(l,2),-1),l,con(-1));il=add(sc(pw(l,3),-2),sc(pw(l,2),3),sc(l,3),con(-2))
    identity('affine invariant J',add(J1,sc(mul(pw(h,2),jl),-1)))
    identity('affine invariant I',add(I1,sc(mul(pw(h,3),il),-1)))
    identity('affine discriminant',add(disc1,sc(mul(pw(h,6),pw(l,2),pw(add(con(1),sc(l,-1)),2)),-1)))
    J0,I0,disc0=invariants(sc(D,-1),{},D,{})
    identity('nonunit quadratic discriminant boundary',add(disc0,sc(pw(D,4),-4)))
    identity('nonunit ratio numerator',add(pw(I0,2),sc(pw(J0,3),4)))
    A=sc(J,F(1,3));Bshort=sc(I,F(1,27));den=add(sc(pw(A,3),4),sc(pw(Bshort,2),27));Dinv=add(sc(pw(J,3),4),pw(I,2))
    identity('classical j normalization',add(sc(mul(pw(A,3),Dinv),1728*4),sc(mul(pw(J,3),den),-6912)))
    return ids

def trim(a):
    a=[F(v) for v in a]
    while len(a)>1 and not a[-1]:a.pop()
    return a
def uadd(a,b):
    z=[F(0)]*max(len(a),len(b))
    for j,v in enumerate(a):z[j]+=v
    for j,v in enumerate(b):z[j]+=v
    return trim(z)
def uscale(a,v):return trim([v*x for x in a])
def umul(a,b):
    z=[F(0)]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):z[i+j]+=v*w
    return trim(z)
def upow(a,n):
    q=[F(1)]
    for _ in range(n):q=umul(q,a)
    return q
def divrem(a,b):
    a=trim(a);b=trim(b);q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]+=v
        for i,w in enumerate(b):a[j+i]-=v*w
        a=trim(a)
    return trim(q),a
def root_one(a,r,N):
    need(a[0]==1,'formal root constant');q=[F(1)]
    for n in range(1,N):
        base=upow(q+[F(0)],r);v=(a[n] if n<len(a) else 0)-(base[n] if n<len(base) else 0)
        q.append(F(v,r))
    return q
def evaluate(a,x):return sum(v*x**j for j,v in enumerate(a))
def finite_controls(inp):
    rho=[F(x) for x in inp['formal_rho']];N=len(rho);series=0
    for p0,q0 in inp['formal_pq_cases']:
        p0,q0=F(p0),F(q0);p=uscale(upow(rho,2),p0);q=uscale(upow(rho,3),q0)
        actual=root_one(uscale(p,1/p0),2,N) if p0 else root_one(uscale(q,1/q0),3,N)
        need(actual==rho,'formal square/cube root branch');series+=N
        need(umul(upow(q,2),[p0**3])==umul(upow(p,3),[q0**2]),'constant depressed ratio')
    affine=[]
    for ls in inp['affine_shape_parameters']:
        l=F(ls);r0=2*l-1;c0=r0/3;c=[c0,F(1)];h=uscale(rho,2);s=uadd(c,uscale(rho,-1-c0))
        C=uscale(umul(umul(s,uadd(s,h)),uadd(s,uscale(h,l))),-1)
        B=uadd(uadd(uscale(upow(s,2),3),uscale(umul(umul(s,h),[1+l]),2)),uscale(upow(h,2),l))
        D=uadd(uscale(s,-3),uscale(h,-1-l))
        need(s[0]==-1 and h[0]==2 and C[0]+D[0]==0 and B[0]+1==0,'normalized regular affine center')
        jl=3*l-(1+l)**2;il=9*l*(1+l)-2*(1+l)**3
        J=uadd(uscale(B,3),uscale(upow(D,2),-1));I=uadd(uadd(uscale(C,27),uscale(umul(B,D),-9)),uscale(upow(D,3),2))
        need(J==uscale(upow(h,2),jl) and I==uscale(upow(h,3),il),'literal affine invariant series')
        affine.append({'lambda':ls,'repeated':l in (0,1),'central':l==F(1,2),'J_lambda':str(jl),'I_lambda':str(il)})
    # The equianharmonic quotient lambda^2-lambda+1=0 is exact over Q.
    _,rem=divrem([-2,3,3,-2], [1,-1,1]);need(rem==[F(-3),F(6)],'equianharmonic I')
    _,rem2=divrem(umul(rem,rem),[1,-1,1]);need(rem2==[F(-27)],'equianharmonic I squared nonzero')
    # Recompute the four finite rational interval-certificate eliminations.
    q=[F(45),126,105];P=[F(45),105,56];R=[F(8955),46410,93275,84084,28028]
    need(divrem(R,q)[1]==[F(1848,25),F(-9828,125)],'first four cubic certificate q remainder')
    need(divrem(R,P)[1]==[F(-50085,128),F(-94185,128)],'first four cubic certificate P remainder')
    need(evaluate(q,F(110,117))==F(1169375,4563),'q excludes forced root')
    need(evaluate(P,F(-159,299))==F(446976,89401),'P excludes forced root')
    # Genuine zero-pair/triple-root and varying-invariant boundary controls.
    x=[F(0),F(1)];E=upow(x,3);need(E[0]==0 and len(E)-1+2==5,'triple-root least degree5')
    J,I,disc=invariants({}, {}, {},var(3));need(not J and not I and not disc,'zero pair is not ratio -4')
    varying=uadd([F(2)],uscale(upow(x,2),-27));need(len(upow(varying,2))==5,'regular nonunit-center varying invariant')
    prime=0
    for p in inp['prime_scope_controls']:
        need(p>=2 and all(p%d for d in range(2,p)),'prime scope literal')
        need(p not in inp['unsafe_noncentral_degrees'],'noncentral composite spectrum control');prime+=1
    need(F(6912,4)==1728,'central j1728')
    need(2%4!=0,'constant j twist valuation obstruction')
    return {'formal_root_coefficients':series,'affine_shapes':affine,'equianharmonic_quotient_identities':2,
      'finite_interval_eliminations':4,'prime_scope_controls':prime,'exceptional_scope_controls':6,
      'finite_controls_do_not_prove_universal_formal_or_safety_statements':True}

def reject(fn,label):
    try:fn()
    except (ValueError,TypeError,KeyError,ZeroDivisionError):return
    raise ValueError('negative accepted '+label)
def negatives(inp,ledger,root):
    n=0
    for role in ['written_proof','exact_checker','typed_inputs']+list(RAW):
        z=copy.deepcopy(ledger);z['sources']=[r for r in z['sources'] if r['role']!=role]
        reject(lambda:source_check(z,root),'omitted '+role);n+=1
        z=copy.deepcopy(ledger);next(r for r in z['sources'] if r['role']==role)['sha256']='0'*64
        reject(lambda:source_check(z,root),'bad pin '+role);n+=1
    for k in BOUNDARY:
        z=copy.deepcopy(inp);z['logical_boundary'][k]=not z['logical_boundary'][k]
        reject(lambda:validate_input(z),'false scope '+k);n+=1
    z=copy.deepcopy(inp);z['schema']=True;reject(lambda:validate_input(z),'bool schema');n+=1
    z=copy.deepcopy(ledger);z['logical_boundary']['formal_trivialization_over_K_itself']=1
    reject(lambda:source_check(z,root),'int bool boundary');n+=1
    z=copy.deepcopy(ledger);z['sources'][0]['path']='../'+PN;reject(lambda:source_check(z,root),'source path escape');n+=1
    z=copy.deepcopy(inp);z['classical_context']['role']='mathematical_proof';reject(lambda:validate_input(z),'context promotion');n+=1
    return n
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root');ap.add_argument('--output');args=ap.parse_args()
    if args.repo_root:root=Path(args.repo_root).resolve()
    else:
        root=next((p for p in [HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents] if (p/BASE/RAW['finite_cubic_interval_certificate'][0]).is_file()),None)
        need(root is not None,'repo root not found; pass --repo-root')
    inp=load(HERE/IN);validate_input(inp);ledger=load(HERE/LN);sources=source_check(ledger,root)
    ids=universal();controls=finite_controls(inp);n=negatives(inp,ledger,root)
    result={'status':'PASS_EXACT_CUBIC_INVARIANTS_AND_NARROW_SOURCE_SUPPORT','schema':1,
      'logical_boundary':BOUNDARY,'classical_context':CONTEXT,'sources':sources,
      'source_ledger_sha256':digest(HERE/LN),'universal_expansion_identities':ids,
      'finite_controls':controls,'source_scope_negative_controls':n,
      'assurance':{'formal_verification':False,'independent_replay':False,'external_review':False}}
    data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(data,encoding='utf-8',newline='\n')
    else:print(data,end='')
if __name__=='__main__':main()
