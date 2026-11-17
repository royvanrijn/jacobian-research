#!/usr/bin/env python3
"""Complete finite polynomial/matrix controls for the fixed three-torsion field."""
from hashlib import sha256
from itertools import permutations,product
from pathlib import Path
import json
import platform

HERE=Path(__file__).resolve().parent
A=-10758382;B=6600759744

def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return trim(r)
def rem(a,b,p):
    a=trim([x%p for x in a]);b=trim([x%p for x in b]);ib=pow(b[-1],-1,p)
    while a!=[0] and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]*ib%p
        for j,x in enumerate(b):a[j+k]=(a[j+k]-c*x)%p
        a=trim(a)
    return a
def pgcd(a,b,p):
    while trim(b)!=[0]:a,b=b,rem(a,b,p)
    return trim([x*pow(a[-1],-1,p)%p for x in a])
def mm(a,b):
    return ((a[0]*b[0]+a[1]*b[2])%3,(a[0]*b[1]+a[1]*b[3])%3,
            (a[2]*b[0]+a[3]*b[2])%3,(a[2]*b[1]+a[3]*b[3])%3)
def det(a):return (a[0]*a[3]-a[1]*a[2])%3
def pc(a,b):return tuple(a[b[i]] for i in range(4))
def pinv(a):return tuple(a.index(i) for i in range(4))
def closure(gs):
    seen={(0,1,2,3)};todo=list(seen)
    while todo:
        x=todo.pop()
        for g in gs:
            y=pc(x,g)
            if y not in seen:seen.add(y);todo.append(y)
    return seen
def sign(p):return 1 if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0 else 2
def line(v):
    if v[0]%3:return (1,v[1]*pow(v[0],-1,3)%3)
    return (0,1)

def replay():
    checks=[]
    def check(n,v):
        if not v:raise AssertionError(n)
        checks.append(n)
    f=[B,A,0,1];df=[A,0,3]
    psi=[-A*A,12*B,6*A,0,3]
    duplication=add(mul(df,df),scale([0]+f,-12))
    check('unrestricted duplication polynomial identity',duplication==scale(psi,-1))
    g=[-27*A*A,108*B,18*A,0,1]
    check('literal monic integral transform',g==[-3125055147963948,712882052352,-193650876,0,1])
    check('scaled variable coefficient identity',all(g[i]*3**i==27*psi[i] for i in range(5)))
    g11=[10,8,10,0,1]
    check('exact reduction modulo eleven',[x%11 for x in g]==g11)
    divisors=[]
    for d in [1,2]:
        for coeff in product(range(11),repeat=d):
            if rem(g11,list(coeff)+[1],11)==[0]:divisors.append(coeff)
    check('all 132 monic linear and quadratic divisors excluded',not divisors)
    derivative=[i*g[i] for i in range(1,len(g))]
    check('squarefree reduction eleven',pgcd(g,derivative,11)==[1])
    lin=[-6,1];cub=[5,-6,6,1]
    check('complete factor product thirteen',[(x-y)%13 for x,y in zip(mul(lin,cub),g)]==[0]*5)
    roots=[x for x in range(13) if sum(c*pow(x,i,13) for i,c in enumerate(cub))%13==0]
    check('all thirteen cubic roots excluded',not roots)
    check('squarefree reduction thirteen',pgcd(g,derivative,13)==[1])
    S4=list(permutations(range(4)));odd=[p for p in S4 if sign(p)==2]
    check('every odd permutation normally generates S4',all(len(closure([pc(pc(q,p),pinv(q)) for q in S4]))==24 for p in odd))
    mats=[m for m in product(range(3),repeat=4) if det(m)]
    check('complete GL2F3 matrix enumeration',len(mats)==48)
    lines=[(1,0),(0,1),(1,1),(1,2)]
    actions={m:tuple(lines.index(line(((m[0]*v[0]+m[1]*v[1])%3,(m[2]*v[0]+m[3]*v[1])%3))) for v in lines) for m in mats}
    check('all twenty four projective permutations occur',set(actions.values())==set(S4))
    I=(1,0,0,1);minusI=(2,0,0,2)
    check('projective action kernel is exactly plus minus identity',{m for m in mats if actions[m]==(0,1,2,3)}=={I,minusI})
    check('determinant equals projective sign for every matrix',all(det(m)==sign(actions[m]) for m in mats))
    check('only determinant one nonidentity involution is minus identity',{m for m in mats if m!=I and mm(m,m)==I and det(m)==1}=={minusI})
    detminus=[m for m in mats if det(m)==2]
    sigma=[m for m in detminus if det(((1-m[0])%3,-m[1]%3,-m[2]%3,(1-m[3])%3))!=0]
    fixed=[m for m in detminus if m not in sigma]
    check('determinant minus one matrix count',len(detminus)==24)
    check('no fixed nonzero vector matrix count',len(sigma)==12)
    check('fixed vector complement count',len(fixed)==12)
    diagonal=(1,0,0,2)
    central=[m for m in mats if mm(m,diagonal)==mm(diagonal,m)]
    check('reflection centralizer order',len(central)==4)
    check('reflection complement consists exactly of involutions',all(mm(m,m)==I for m in fixed))
    check('exact joint-density matrix factor',len(sigma)*4==len(mats))
    check('curve good at three',(-16*(4*A**3+27*B**2))%3!=0)
    return dict(passed=True,exact_controls=len(checks),checks=checks,
                finite_polynomial_primes=[11,13],monic_linear_quadratic_divisors_checked=132,
                GL2_matrices=48,projective_permutations=24,sigma3_matrices=12,
                python=platform.python_version(),arithmetic='standard-library exact integers and finite fields',
                finite_controls_certify_imported_theorems=False,
                formal_verification=False,external_review=False,independent_replay=False)

if __name__=='__main__':
    report=replay()
    report['input_sha256']={name:sha256((HERE/name).read_bytes()).hexdigest() for name in
        ['ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md',Path(__file__).name]}
    (HERE/'elliptic_exception_three_torsion_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
