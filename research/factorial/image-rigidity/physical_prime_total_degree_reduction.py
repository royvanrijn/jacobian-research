#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,math,platform
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
COUNTS={}
def count(key,n=1):COUNTS[key]=COUNTS.get(key,0)+n
def vf(n,p):
    ans=0
    while n:n//=p;ans+=n
    return ans
def vp(n,p):
    if not n:return None
    ans=0
    while n%p==0:n//=p;ans+=1
    return ans
def raw(n,p,k):
    return math.comb(n,k)*math.factorial(n+(p-1)*k-1)//math.factorial(n-1)
def hcoeff(n,p,mod=None):
    out=[p];u=1
    for k in range(1,n+1):
        for j in range(p*(k-1)+1,p*k+1):
            if j%p:u*=j
        if mod:u%=mod
        N=n+(p-1)*k-1
        A=math.comb(N,p*k) if N>=p*k else 0
        B=math.comb(N,p*k-1)
        z=u*(p*A+B)
        out.append(z if mod is None else z%mod)
    return out
def shift(c,a,mod):
    return [sum(c[k]*math.comb(k,j)*pow(a,k-j,mod)
                for k in range(j,len(c)))%mod
            for j in range(len(c))]
def kpoly(n,p,Y):
    if not n:return s.Integer(p)
    return p+sum((p*(math.comb(n+(p-1)*k-1,p*k)
                      if n+(p-1)*k-1>=p*k else 0)
                 +math.comb(n+(p-1)*k-1,p*k-1))*Y**k
                for k in range(1,n+1))
Y=s.Symbol('Y')
for p in (2,3,5,7,11):
    for n in range(1,15):
        h=hcoeff(n,p)
        assert h[0]==p and h[-1]%p
        count('integral_constant_and_leading_controls')
        for k,z in enumerate(h):
            assert z*p**k==p*raw(n,p,k)
            count('literal_original_scaled_coefficient_identities')
    polys=[]
    for n in range(1,26):
        h=hcoeff(n,p,p)
        f=s.Poly.from_list(h[::-1],Y,modulus=p)
        target=-(s.Poly(Y,Y,modulus=p)*
                 (polys[n-2] if n>1 else s.Poly(0,Y,modulus=p)))
        if n>p:target+=polys[n-p-1]
        if n==1:target-=s.Poly(Y,Y,modulus=p)
        assert f==target
        assert all(k%p==n%p for k,c in enumerate(h) if c%p)
        polys.append(f)
        count('finite_field_generating_function_coefficient_controls')

indices=[]
for p,maxr in ((2,5),(3,4),(5,3),(7,3),(11,2)):
    for r in range(maxr):
        for b in range(1,p):
            n=b*p**r
            if n==1:continue
            for k in range(1,n):
                v=vf(n,p)-vf(k,p)-vf(n-k,p)+vf(n+(p-1)*k-1,p)-vf(n-1,p)
                assert v>=k
                count('Eisenstein_interior_valuation_controls')
            top=vf(p*n-1,p)-vf(n-1,p)
            assert top==n-1
            indices.append([p,b,r,n])
            count('Eisenstein_index_instances')
p=7;n=343
for k in range(1,n):
    v=vf(n,p)-vf(k,p)-vf(n-k,p)+vf(n+6*k-1,p)-vf(n-1,p)
    assert v>=k
    count('Eisenstein_interior_valuation_controls')
assert vf(7*n-1,7)-vf(n-1,7)==n-1
indices.append([7,1,3,343]);count('Eisenstein_index_instances')

for p in (2,3,5,7):
    K=[kpoly(n,p,Y) for n in range(13)]
    D=[s.Integer(1),-p-Y]+[s.Integer((-1)**j*math.comb(p,j))
                               for j in range(2,p+1)]
    for n in range(13):
        lhs=s.expand(sum(D[j]*K[n-j] for j in range(min(n,p)+1)))
        rhs=(p-n)*D[n] if n<=p else 0
        assert s.expand(lhs-rhs)==0
        count('exact_power_sum_rational_generating_function_coefficients')
for p in (3,5,7,11):
    u=1;w=math.factorial(p-1)
    for k in range(1,51):
        for j in range(p*(k-1)+1,p*k+1):
            if j%p:u=u*j%(p*p)
        assert u==pow(w,k,p*p)
        count('mod_prime_squared_trace_unit_controls')

c15=hcoeff(15,7,343)
b15=shift(c15,3,343)
expected15=[49,210,182,14,70,280,273,153]
def check_displayed_15(data):
    assert data[:8]==expected15
check_displayed_15(b15)
assert [vp(z,7) for z in b15[:8]]==[2,1,1,1,1,1,1,0]
f15=s.Poly.from_list([z%7 for z in c15][::-1],Y,modulus=7)
assert f15==s.Poly(6*Y*(Y-3)**7*(Y-5)**7,Y,modulus=7)
count('exact_n15_unit_cluster_controls',3)
bad=b15[:];bad[0]=(bad[0]+1)%343
try:check_displayed_15(bad)
except AssertionError:count('corrupted_display_rejected')
else:raise AssertionError('corruption was accepted')
c145=hcoeff(145,7,2401);b145=shift(c145,4,2401)
assert b145[:2]==[1470,1176] and [vp(z,7) for z in b145[:2]]==[2,2]
assert [z%7 for z in b145[:15]]==[0]*14+[5]
count('exact_n145_transversality_obstruction_controls',2)

out={
    'passed':True,
    'scope':'Uniform written prime-total-degree scaling/Eisenstein/radius reductions; bounded controls are not an all-index simplicity proof.',
    'ordinary_gap_six_all_indices_claimed':False,
    'all_integral_Gamma_shapes_claimed':False,
    'Eisenstein_theorem':'physical gamma=n, gap=p-1, n=b*p^r, 1<=b<p',
    'critical_radius_theorem':'any hypothetical multiple root u=epsilon/t has v_p(u)=1',
    'controls':COUNTS,
    'Eisenstein_index_instances':indices,
    'n15_shift_first8_mod343':b15[:8],
    'n145_shift_first17_mod2401':b145[:17],
    'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in (HERE/'PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex',Path(__file__).resolve())},
    'python':platform.python_version(),'sympy':s.__version__,
    'command':'python research/factorial/image-rigidity/physical_prime_total_degree_reduction.py',
    'independent_replay':False,'formal_verification':False,'external_review':False,
}
(HERE/'physical_prime_total_degree_reduction_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
