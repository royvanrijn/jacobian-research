#!/usr/bin/env python3
"""Exact original-contraction support for the prime balanced nonpower theorem."""
if not __debug__:
    raise SystemExit('optimization is unsupported: exact source/control guards are required')
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
from math import factorial, prod
import argparse, copy, hashlib, json

HERE=Path(__file__).resolve().parent
PN='PRIME_BALANCED_IMAGE_FAILURES_ARE_ABSOLUTELY_NONPOWER.md'
CN='verify_prime_balanced_image_nonpower.py'
IN='prime_balanced_image_nonpower_inputs.json'
LN='prime_balanced_image_nonpower_source_ledger.json'
PROOF='46669b6dc70e08b197607a5f10ca70f08016e35ca8e7f3540501346344d15642'
BOUNDARY={
 'all_characteristic_zero_extensions':True,'original_full_polynomial_markers':True,
 'fixed_r_multiples_and_eventual_multiples':True,'nilpotence_safety_both_directions':True,
 'prime_counterexample_absolute_scalar_nonpower':True,'radial_primitivity_required':False,
 'purity_alone_nonpower':False,'composite_automatic_nonpower':False,
 'all_nonzero_input_weights_annihilated':False,'arbitrary_shifted_progression':False,
 'positive_characteristic':False,'unconditional_witness_existence':False,
 'squarefree_or_irreducible_conclusion':False,'two_pair_five_seven_resolution':False,
 'external_mathematical_premise':False,'external_review':False,
 'formal_verification':False,'independent_replay':False}
CONTEXT={'title':'A combinatorial approach to matrix algebra','author':'Doron Zeilberger',
 'year':1985,'journal':'Discrete Mathematics','volume':56,'pages':'61-72',
 'url':'https://sites.math.rutgers.edu/~zeilberg/mamarimY/DM85.pdf',
 'sha256':'6378bd7078732326b6b228d2d27e490f3cde33a278bbec40ba6d3a79640b9e5f',
 'locator':'Section2 printed63 equation2','role':'classical_context_only'}
MATRICES=[[[2]],[[1,2],[3,4]],[[0,1],[0,0]],[[1,1],[-1,-1]],
 [[2,1,-1],[3,0,2],[0,1,1]],[[0,1,2],[0,0,3],[0,0,0]]]
EXPECTED={'schema':1,'logical_boundary':BOUNDARY,'classical_context':CONTEXT,
 'matrix_controls':MATRICES,'moment_orders':list(range(8)),
 'filter_orders':list(range(1,8)),'filter_degrees':list(range(22)),
 'sharp_marker_degrees':list(range(6)),'prime_scope_controls':[2,3,5,7],
 'composite_seed_orders':[1,2,3]}

def need(ok,msg):
    if not ok: raise ValueError(msg)
def strict(a,b):
    return type(a) is type(b) and (set(a)==set(b) and all(strict(a[k],b[k]) for k in a)
      if isinstance(a,dict) else len(a)==len(b) and all(strict(x,y) for x,y in zip(a,b))
      if isinstance(a,list) else a==b)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding='utf-8'))
def validate_input(d): need(strict(d,EXPECTED),'typed immutable inputs/scope mismatch')
def source_check(ledger):
    need(type(ledger) is dict and set(ledger)=={'schema','logical_boundary','sources'},'ledger shape')
    need(type(ledger['schema']) is int and ledger['schema']==1,'ledger schema')
    need(strict(ledger['logical_boundary'],BOUNDARY),'typed ledger boundary')
    roles={'written_proof':(PN,'mathematical_proof'),'exact_checker':(CN,'identity_and_source_support'),
           'typed_inputs':(IN,'exact_control_data')}
    rows=ledger['sources']; need(type(rows) is list and len(rows)==3,'exact three source records')
    found={}
    for row in rows:
        need(type(row) is dict and set(row)=={'role','path','relation','sha256'},'source record shape')
        need(all(type(v) is str for v in row.values()),'source record strings')
        role=row['role']; need(role in roles and role not in found,'source role')
        name,relation=roles[role]
        need(row['path']==name and row['relation']==relation,'source path/relation')
        need(len(row['sha256'])==64 and all(c in '0123456789abcdef' for c in row['sha256']),'sha type')
        need((HERE/name).is_file() and digest(HERE/name)==row['sha256'],'source hash '+role)
        found[role]=dict(row)
    need(set(found)==set(roles),'all source roles')
    need(digest(HERE/PN)==PROOF,'reviewed proof pin')
    return found

# Multivariate sparse Fraction polynomials, with exponent tuples.
def clean(p): return {e:F(c) for e,c in p.items() if c}
def add(a,b):
    o=dict(a)
    for e,c in b.items(): o[e]=o.get(e,F(0))+c
    return clean(o)
def scale(p,c): return clean({e:v*c for e,v in p.items()})
def mul(a,b):
    o={}
    for e,c in a.items():
        for f,d in b.items():
            v=tuple(x+y for x,y in zip(e,f)); o[v]=o.get(v,F(0))+c*d
    return clean(o)
def power(p,k,nv):
    q={(0,)*nv:F(1)}
    for _ in range(k): q=mul(q,p)
    return q
def variable(i,nv):
    e=[0]*nv;e[i]=1;return {tuple(e):F(1)}
def contraction(p,n):
    o={}
    for e,c in p.items():
        a,b=e[:n],e[n:]
        if any(x>y for x,y in zip(a,b)): continue
        f=tuple(y-x for x,y in zip(a,b))
        v=c*prod(factorial(y)//factorial(y-x) for x,y in zip(a,b))
        o[f]=o.get(f,F(0))+v
    return clean(o)
def bilinear(B):
    n=len(B);o={}
    for i in range(n):
        for j in range(n):
            if B[i][j]:
                e=[0]*(2*n);e[i]=1;e[n+j]=1;o[tuple(e)]=F(B[i][j])
    return o
def upmul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def determinant_I_minus_t(B):
    n=len(B);out=[F(0)]*(n+1)
    for pp in permutations(range(n)):
        sign=(-1)**sum(pp[i]>pp[j] for i in range(n) for j in range(i+1,n))
        term=[F(sign)]
        for i,j in enumerate(pp):term=upmul(term,[F(i==j),-F(B[i][j])])
        for i,x in enumerate(term):out[i]+=x
    need(out[0]==1,'determinant constant')
    return out
def inverse_series(D,N):
    c=[F(1)]
    for k in range(1,N+1):c.append(-sum(D[j]*c[k-j] for j in range(1,min(k,len(D)-1)+1)))
    return c
def mm(A,B):
    return [[sum(F(A[i][k])*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def inv(A):
    n=len(A);a=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        u=a[j][j];a[j]=[x/u for x in a[j]]
        for i in range(n):
            if i!=j:
                u=a[i][j];a[i]=[x-u*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]
def substitute(p,images,nv):
    o={}
    for e,c in p.items():
        term={(0,)*nv:c}
        for k,pow_ in enumerate(e):term=mul(term,power(images[k],pow_,nv))
        o=add(o,term)
    return o
def linear_image(row,offset,nv):
    return clean({next(iter(variable(offset+j,nv))):F(v) for j,v in enumerate(row) if v})
def transform(p,S):
    n=len(S);I=inv(S)
    images=[linear_image([I[j][i] for j in range(n)],0,2*n) for i in range(n)]
    images += [linear_image(row,n,2*n) for row in S]
    return substitute(p,images,2*n)
def output_transform(p,S):return substitute(p,[linear_image(row,0,len(S)) for row in S],len(S))
def encoded(p):return [[list(e),str(c)] for e,c in sorted(p.items())]
def p_hash(p):return hashlib.sha256(json.dumps(encoded(p),separators=(',',':')).encode()).hexdigest()

# Exact univariate cyclotomic reduction; no floating roots of unity.
def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a
def divrem(a,b):
    a=trim(a);b=trim(b);q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        k=len(a)-len(b);c=a[-1]/b[-1];q[k]+=c
        for j,v in enumerate(b):a[k+j]-=c*v
        a=trim(a)
    return trim(q),a
def cyclotomic(r):
    p=[F(-1)]+[F(0)]*(r-1)+[F(1)]
    for d in range(1,r):
        if r%d==0:p,rem=divrem(p,cyclotomic(d));need(not any(rem),'cyclotomic division')
    return p
def expect_failure(fn,label):
    try:fn()
    except (ValueError,TypeError,KeyError):return
    raise ValueError('negative did not reject '+label)

def run_math(inp):
    rows=[];count=0
    for B in inp['matrix_controls']:
        n=len(B);H=bilinear(B);D=determinant_I_minus_t(B);c=inverse_series(D,max(inp['moment_orders']))
        for k in inp['moment_orders']:
            E=contraction(power(H,k,2*n),n);want={(0,)*n:F(factorial(k))*c[k]} if c[k] else {}
            need(E==want,'literal determinant coefficient');count+=1
        rows.append({'matrix':B,'determinant':[str(v) for v in D],'scalar_coefficients':[str(v) for v in c]})
    filters=0
    for r in inp['filter_orders']:
        phi=cyclotomic(r)
        for k in inp['filter_degrees']:
            a=[F(0)]*(max(1,(r-1)*k+1))
            for j in range(r):a[j*k]+=F(1,r)
            _,rem=divrem(a,phi);need(trim(rem)==[F(int(k%r==0))],'exact root-unity filter');filters+=1
    H=bilinear([[0,1],[0,0]]);sharp=0
    for M in inp['sharp_marker_degrees']:
        Q=power(variable(2,4),M,4)
        for k in range(M+3):
            actual=contraction(mul(Q,power(H,k,4)),2)
            want={(M-k,k):F(factorial(M)//factorial(M-k))} if k<=M else {}
            need(actual==want,'sharp original marker boundary');sharp+=1
    S=[[1,0],[1,1]];GL=0
    markers=[power(variable(2,4),3,4),mul(variable(1,4),power(variable(3,4),2,4)),
             mul(variable(0,4),power(variable(2,4),2,4)),mul(variable(1,4),mul(variable(2,4),variable(3,4)))]
    for Q in markers:
        for k in range(5):
            P=mul(Q,power(H,k,4));need(contraction(transform(P,S),2)==output_transform(contraction(P,2),S),'full GL polynomial equivariance');GL+=1
    Q=mul(variable(1,4),power(variable(3,4),2,4))
    need(contraction(transform(Q,S),2)!=contraction(Q,2),'full E is not scalar invariant')
    # BOTH grading extremes are computed in the full original ring.
    Gs=[add(variable(0,4),mul(variable(0,4),variable(2,4))),
        add(mul(variable(0,4),variable(2,4)),mul(power(variable(0,4),2,4),variable(2,4)))]
    grading=0
    for G in Gs:
        for r in range(2,6):
            P=power(G,r,4)
            for sl in [slice(0,2),slice(2,4)]:
                ds=[sum(e[sl]) for e in G];es=[sum(e[sl]) for e in P]
                need(min(es)==r*min(ds) and max(es)==r*max(ds),'full grading extrema');grading+=1
    # Scope controls: pure prime powers may be safe; composite unsafe powers exist.
    pure_controls=0
    for p in inp['prime_scope_controls']:
        for m in range(1,4):need(not contraction(power(H,p*m,4),2),'pure safe prime power');pure_controls+=1
    tau,w,v,t,z,y=[variable(i,6) for i in range(6)]
    seed=mul(tau,mul(add(t,scale(y,-1)),add(mul(w,z),mul(v,t))))
    mark=mul(tau,y);composite=[]
    for m in inp['composite_seed_orders']:
        pure=power(seed,2*m,6);need(not contraction(pure,3),'composite square pure')
        marked=contraction(mul(mark,pure),3);want={(0,0,0):F(-factorial(2*m+1)*factorial(2*m))}
        need(marked==want,'composite unsafe proper square');composite.append({'m':m,'marked':str(next(iter(marked.values())))})
    # Actual mathematical countermodels, not sampled universal inferences.
    need(contraction(mul(variable(2,4),H),2)=={(0,1):F(1)},'positive input weight countermodel')
    D=determinant_I_minus_t([[1,0],[0,-1]]);need(D==[F(1),F(0),F(-1)],'shifted progression nonnilpotent')
    cs=inverse_series(D,9);need(all(cs[k]==0 for k in [1,3,5,7,9]) and cs[2]==1,'shifted progression zero samples')
    need(all(factorial(2*m)%2==0 for m in range(1,5)),'characteristic-two sampled countermodel')
    return {'determinant_literal_controls':count,'root_unity_fraction_controls':filters,
      'sharp_full_marker_controls':sharp,'full_GL_polynomial_controls':GL,
      'full_two_grading_controls':grading,'pure_safe_prime_power_controls':pure_controls,
      'composite_proper_square_controls':composite,'actual_math_scope_countermodels':6,
      'determinant_rows':rows,'actual_composite_seed_sha256':p_hash(seed)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root');ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    inp=load(HERE/IN);validate_input(inp);sources=source_check(load(HERE/LN));result=run_math(inp)
    negative=0
    for key,value in BOUNDARY.items():
        bad=copy.deepcopy(inp);bad['logical_boundary'][key]=not value
        expect_failure(lambda bad=bad:validate_input(bad),'false scope '+key);negative+=1
    bad=copy.deepcopy(inp);bad['schema']=True;expect_failure(lambda:validate_input(bad),'typed schema');negative+=1
    bad=copy.deepcopy(inp);bad['logical_boundary']['all_characteristic_zero_extensions']=1
    expect_failure(lambda:validate_input(bad),'boolean versus integer');negative+=1
    bad=copy.deepcopy(inp);bad['matrix_controls']=[];expect_failure(lambda:validate_input(bad),'empty literal controls');negative+=1
    report={'schema':1,'status':'PASS_PRIME_BALANCED_NONPOWER_AND_LINEAR_SAMPLE_NILPOTENCE_SUPPORT',
      'logical_boundary':BOUNDARY,'source_records':sources,'support':result,
      'typed_source_scope_internal_negatives':negative,
      'finite_controls_prove_universal_theorem':False,'classical_context_is_mathematical_premise':False,
      'external_review':False,'formal_verification':False,'independent_replay':False}
    encoded_report=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(encoded_report,encoding='utf-8')
    else:print(encoded_report,end='')
if __name__=='__main__':main()
