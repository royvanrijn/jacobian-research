#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,platform
HERE=Path(__file__).resolve().parent
PROOF=HERE/"GAP_SIX_BOUNDARY_COUPLED_FROBENIUS.tex"
IMPORTS={'ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex': '9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07', 'PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex': '42c5fc615400a49dd3ae4f19a9a832381607e69ce4f6ff6d255b57ca4cfd156f', 'gap_three_hurwitz_structure.tex': '22c690a45976c9b72b210ff76bbcd3a7084ddba0aae4af0c7f40ce2b109f705c'}
for name,digest in IMPORTS.items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,("Changed exact proof input",name)
COUNTS={}
def count(k,n=1):COUNTS[k]=COUNTS.get(k,0)+n
def addscale(a,b,c):return [x+c*y for x,y in zip(a,b)]
def conv(a,b,N):
    ans=[Q(0)]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]):
                if y:ans[i+j]+=x*y
    return ans
def action(f,c,F):
    N=len(f)-1;ans=[(c-k)*x for k,x in enumerate(f)]
    for d,a in F.items():
        for k in range(N+1-d):ans[k+d]+=a*f[k]
    return ans
def affine(f,c,F,m,b):
    return addscale([b*x for x in f],action(f,c,F),m)
def power(f,c,F,m,b,r):
    for _ in range(r):f=affine(f,c,F,m,b)
    return f
def product(f,c,F,m,constants):
    for b in constants:f=affine(f,c,F,m,b)
    return f
def inverse_AH(N):
    h=[Q(1)]
    for n in range(1,N+1):
        z=1;tot=Q(0)
        while z<=n:tot+=h[n-z];z*=7
        h.append(-tot/n)
    return h
def vp(a):
    if not a:return None
    v=0
    while a.numerator%7==0:a/=7;v+=1
    while a.denominator%7==0:a*=7;v-=1
    return v
L=8;N=7*L+6;E=inverse_AH(N)
assert all(vp(x) is None or vp(x)>=0 for x in E)
count('integral_inverse_Artin_Hasse_coefficients',N+1)
cases=[(1,0,False),(2,0,False),(6,0,False),(7,0,False),
       (8,0,False),(15,0,False),(49,0,False),
       (8,-6,True),(15,1,True)]
bad_control=False
for n,eta,deform in cases:
    c=n+eta+5;a1=Q(7 if deform else 0);a2=Q(-14 if deform else 0)
    F={1:Q(1),7:a1,49:a2}
    b=a1-1
    h=[Q(0)]*(N+1);h[0]=1
    if deform:h[1]=7;h[2]=14
    eh=conv(E,h,N)
    for r in range(5):
        f=power(h,c,F,1,0,r)
        f=product(f,c,F,1,[6*j for j in range(n)])
        gauged=conv(E,f,N)
        for j in range(7):
            e=(c-j)%7;nj=sum(ell%7==e for ell in range(n))
            assert nj==n//7+int(e<n%7)
            cj=Q(c-j+6*e,7)
            Fs={1:Q(1),7:(a2-1)*7**6/b**7}
            hj=[eh[7*k+j]*(7/b)**k for k in range(L+1)]
            assert vp(hj[0])==0
            assert all(not x or vp(x)>=1 for x in hj[1:])
            constants=[6*(ell-e) for ell in range(n) if ell%7!=e]
            hh=product(hj,cj,Fs,7,constants)
            assert vp(hh[0])==0
            assert all(not x or vp(x)>=1 for x in hh[1:])
            right=power(hh,cj,Fs,7,-6*e,r)
            right=product(right,cj,Fs,1,[6*k for k in range(nj)])
            right=[7**nj*x for x in right]
            left=[gauged[7*k+j]*(7/b)**k for k in range(L+1)]
            assert left==right,(n,eta,deform,r,j)
            count('exact_residue_section_coefficient_identities',L+1)
            count('actual_descended_unit_seed_controls')
            if not bad_control and r==1 and e:
                wrong=power(hh,cj,Fs,7,6*e,r)
                wrong=product(wrong,cj,Fs,1,[6*k for k in range(nj)])
                wrong=[7**nj*x for x in wrong]
                assert wrong!=left
                bad_control=True;count('corrupted_affine_shift_rejected')
assert bad_control
physical={}
for eta in range(-6,2):
 for r in range(7):
  match=[]
  for j in range(7):
   e=(r+eta+5-j)%7;d=(r+eta+5-j)//7
   difference=e+d-5-int(e<r)
   assert -6<=difference<=1
   if not eta and not difference:match.append(j)
   count('complete_offset_residue_window_states')
  if not eta:physical[r]=match
assert physical=={0:[0,6],1:[1],2:[2],3:[3],4:[4],5:[5],6:[]}
count('complete_physical_section_table_rows',7)
# u=7v,w^6=7v^7, b=-1: two NONZERO sections cancel.
for v in (1,2,3,8):
 assert -7*Q(v)+Q(1,v**6)*(7*v**7)==0
 assert -7*v and Q(1,v**6)
 count('non_splitting_orbit_controls')

import math
def binomial_fraction(x,n):
    a=Q(1)
    for k in range(n):a*=x-k
    return a/math.factorial(n)
def digitsum(j):
    s=0
    while j:s+=j%7;j//=7
    return s
for n in range(1,11):
    for exponent in range(5):
        D=n+exponent;c=n+5
        f=[Q(0)]*(D+1);f[0]=1
        f=power(f,c,{1:Q(1)},1,0,exponent)
        f=product(f,c,{1:Q(1)},1,[6*k for k in range(n)])
        values=[binomial_fraction(Q(c-ell,6)+n-1,n)*(c-ell)**exponent
                for ell in range(D+1)]
        for j in range(D+1):
            dj=sum((-1)**ell*math.comb(j,ell)*values[ell] for ell in range(j+1))
            assert dj==Q(math.factorial(j)*f[j],6**n*math.factorial(n))
            assert not dj or vp(dj)>=0
            a=j%6;ell=(j-a)//6
            weight=Q(7**ell*math.factorial(a),math.factorial(j))
            assert vp(weight)==(digitsum(j)-a)//6>=0
            count('exact_finite_difference_lattice_coefficient_controls')
        assert (7*n-1)*pow(6,-1,7)%7==1
        count('physical_binomial_first_digit_controls')
n=3;D=7;f=[Q(0)]*(D+1);unit=f[:];unit[0]=1
for a,z in enumerate([360,-342,119,-18,1]):
    column=power(unit[:],8,{1:Q(1)},1,0,a)
    column=product(column,8,{1:Q(1)},1,[0,6,12])
    f=addscale(f,column,z)
for a in range(6):
    for ell in range((D-a)//6+1):
        z=Q(math.factorial(a)*7**ell*f[a+6*ell],6**3*math.factorial(3))
        assert not z or vp(z)>=1
        count('primitive_first_lattice_zero_controls')


# The first-layer zero is a complete polynomial-vector relation.
# Mutating one constant coefficient must not be accepted as that relation.
mutated=f[:]
column=product(unit[:],8,{1:Q(1)},1,[0,6,12])
mutated=addscale(mutated,column,1)
assert any(z and vp(z)==0 for z in
    [Q(math.factorial(a)*7**ell*mutated[a+6*ell],6**3*math.factorial(3))
     for a in range(6) for ell in range((D-a)//6+1)])
count('corrupted_first_layer_relation_rejected')

# Generic Laurent divided differences in the free two-variable Laurent ring.
def ladd(p,q,scale=1):
    ans=p.copy()
    for mon,c in q.items():
        ans[mon]=ans.get(mon,Q(0))+scale*c
        if not ans[mon]:del ans[mon]
    return ans
def lmul(p,q):
    ans={}
    for (a,b),c in p.items():
        for (e,f),d in q.items():
            mon=(a+e,b+f);ans[mon]=ans.get(mon,Q(0))+c*d
    return {m:c for m,c in ans.items() if c}
xy={(1,0):Q(1),(0,1):Q(-1)}
for exponent in range(-7,8):
    if exponent>0:
        divided={(exponent-1-j,j):Q(1) for j in range(exponent)}
    elif exponent<0:
        m=-exponent
        divided={(-1-j,j-m):Q(-1) for j in range(m)}
    else:divided={}
    target=ladd({(exponent,0):Q(1)},{(0,exponent):Q(1)},-1)
    assert lmul(xy,divided)==target
    count('formal_Laurent_divided_difference_identities')
def eval_laurent(coeffs,x):
    return sum(a*x**j for j,a in coeffs.items())
affine={0:Q(3),1:Q(2),2:Q(7),-1:Q(14),-4:Q(-21)}
for xx,yy in [(1,2),(1,8),(2,9),(1,-6),(8,15)]:
    xx,yy=Q(xx),Q(yy)
    assert vp(xx)==vp(yy)==0
    divided=(eval_laurent(affine,xx)-eval_laurent(affine,yy))/(xx-yy)
    assert vp(divided)==0
    count('affine_Laurent_unit_disk_controls')
assert (Q(1)**2-Q(-1)**2)/(Q(1)-Q(-1))==0
count('nonaffine_false_injectivity_rejected')
# A complete symbolic integral column lift with determinant exactly one.
vv={(1,0):Q(1)}
entry01={(1,0):Q(7)}
entry10=vv
entry11={(0,0):Q(1),(2,0):Q(7)}
det=ladd(entry11,lmul(entry01,entry10),-1)
assert det=={(0,0):Q(1)}
count('pointwise_unit_column_operation_identity')
# Actual characteristic-zero moment bridge coefficients.
def rising(a,k):
    ans=1
    for j in range(k):ans*=a+j
    return ans
for n in range(2,9):
    for j in range(n+1):
        full=math.comb(n,j)*rising(n,6*j)
        lower=math.comb(n-1,j)*rising(n,6*j) if j<n else 0
        upper=math.comb(n,j)*rising(n+1,6*j)
        assert Q(n-j,n)*full==lower
        assert Q(n+6*j,n)*full==upper
        count('literal_actual_moment_bridge_coefficients',2)

out={
 'passed':True,
 'scope':{
  'exact_coupled_stage_all_indices':True,
  'closed_eight_offset_window':True,
  'integral_finite_difference_lattice':True,
  'primitive_first_layer_obstruction':True,
  'conditional_affine_saturation_criterion':True,
  'uniform_saturation_certificate':False,
  'all_index_gap_six_SFC':False,
  'independent_section_vanishing':False,
  'all_binomial_SFC':False,
  'unrestricted_FC_or_SFC':False,
 },
 'controls':COUNTS,'total_controls':sum(COUNTS.values()),
 'physical_section_table':physical,
 'source_sha256':{
   PROOF.name:hashlib.sha256(PROOF.read_bytes()).hexdigest(),
   Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   **IMPORTS,
 },
 'exact_imports':{
  'Artin_Hasse_closed_class':'ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex lem:artin-hasse-closure',
  'all_conjugate_boundary_radius':'PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex thm:radius',
  'critical_rank_identity':'gap_three_hurwitz_structure.tex lem:critical-ideal; rational polynomial transfer',
 },
 'python':platform.python_version(),
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'novelty_claim':False,'priority_claim':False,
 'command':'.venv/bin/python research/factorial/image-rigidity/gap_six_boundary_coupled_frobenius.py',
}
(HERE/'gap_six_boundary_coupled_frobenius_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
