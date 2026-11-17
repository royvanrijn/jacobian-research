#!/usr/bin/env python3
"""Exact induced-action and cyclotomic controls for written general proofs."""
from pathlib import Path
from itertools import permutations,product
from math import lcm
import hashlib,importlib.util,json,platform,time
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('retained_permutation_helpers',HERE.parent/'universal-keller-homogeneous-receivers/verify.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
mul,inv,generated,primitive=old.mul,old.inv,old.generated,old.primitive
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def order(g):
    e=tuple(range(len(g)));h=e
    for i in range(1,1+len(g)*len(g)):
        h=mul(g,h)
        if h==e:return i
    raise AssertionError('unexpected group order')
def induced(G,h,modulus):
    n=len(h);identity=tuple(range(n));p=order(h);powers=[identity]
    for j in range(1,p):powers.append(mul(h,powers[-1]))
    H=set(powers);cosets=[];reps=[];left=set(G)
    while left:
        a=min(left);C=frozenset(mul(a,b) for b in H);cosets.append(C);reps.append(a);left-=C
    where={g:i for i,C in enumerate(cosets) for g in C};data={}
    for g in G:
        pi=[];t=[0]*len(reps)
        for j,a in enumerate(reps):
            v=mul(g,a);i=where[v];pi.append(i);local=mul(inv(reps[i]),v);t[i]=powers.index(local)*(modulus//p)
        data[g]=(tuple(pi),tuple(t))
    return data
def join_actions(parts,G):
    out={}
    for g in G:
        pi=[];t=[];offset=0
        for data in parts:
            perm,shift=data[g];pi.extend(offset+j for j in perm);t.extend(shift);offset+=len(perm)
        out[g]=(tuple(pi),tuple(t))
    return out
def affine_mul(g,h,m):
    a,t=g;b,u=h;ai=inv(a)
    return mul(a,b),tuple((t[j]+u[ai[j]])%m for j in range(len(t)))
def has_fixed_point(action,m):
    pi,t=action;left=set(range(len(pi)))
    while left:
        j=min(left);cycle=[]
        while j not in cycle:cycle.append(j);j=pi[j]
        left-=set(cycle)
        if sum(t[j] for j in cycle)%m:return False
    return True

def main():
    started=time.monotonic();base_controls=[]
    examples=[('S3',frozenset(permutations(range(3))),[(1,0,2),(1,2,0)],5),('S4',frozenset(permutations(range(4))),[(1,0,2,3),(1,0,3,2),(1,2,0,3)],32)]
    A5=generated([(1,2,0,3,4),(1,2,3,4,0)],5)
    examples.append(('A5',A5,[(1,0,3,2,4),(1,2,0,3,4),(1,2,3,4,0)],62))
    for name,G,hs,dimension in examples:
        modulus=lcm(*(order(h) for h in hs));parts=[induced(G,h,modulus) for h in hs];actions=join_actions(parts,G)
        ck(name+' complete induced action multiplication',all(actions[mul(g,h)]==affine_mul(actions[g],actions[h],modulus) for g in G for h in G))
        ck(name+' claimed split-base dimension',len(next(iter(actions.values()))[0])==dimension)
        identity=tuple(range(len(next(iter(G)))))
        ck(name+' whole torus free on every nonidentity element',all(not has_fixed_point(actions[g],modulus) for g in G if g!=identity))
        base_controls.append({'group':name,'order':len(G),'split_dimension':dimension,'scalar_modulus':modulus,'nonidentity_fixed_loci_empty':len(G)-1})
    # Exact seventh-cyclotomic self-isogeny, not a numerical norm estimate.
    z=s.symbols('z');Phi=s.cyclotomic_poly(7,z);a=z+z*z+z**4;abar=z**3+z**5+z**6
    def multiplication(poly):
        return s.Matrix(6,6,lambda i,j:s.Poly(s.rem(poly*z**j,Phi,z),z).coeff_monomial(z**i))
    A=multiplication(a);B=multiplication(abar);R=multiplication(z)
    ck('Gaussian period trace identity',s.rem(a+abar+1,Phi,z)==0)
    ck('Gaussian period exact norm identity',A*B==2*s.eye(6) and A+B==-s.eye(6))
    ck('degree8 determinant',A.det()==8)
    smith=[abs(int(v)) for v in smith_normal_form(A,domain=s.ZZ).diagonal()]
    ck('degree8 elementary kernel exact Smith factors',smith==[1,1,1,2,2,2])
    ck('isogeny commutes with faithful cyclic holonomy',A*R==R*A and R**7==s.eye(6) and all(R**i!=s.eye(6) for i in range(1,7)))
    q=s.gcd(s.Poly(Phi,z,modulus=2),s.Poly(a,z,modulus=2)).as_expr()
    ck('exact residue quotient cubic',q==z**3+z+1 and s.Poly(q,z,modulus=2).is_irreducible)
    ck('cyclic holonomy preserves global top logarithmic form',R.det()==1)
    # AGL1(F8) full actual permutation monodromy.
    def ff_mul(a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&8:a^=11
        return out
    T=tuple(v^1 for v in range(8));rot=tuple(ff_mul(2,v) for v in range(8));G=generated([T,rot],8)
    expected={tuple(ff_mul(a,v)^b for v in range(8)) for a in range(1,8) for b in range(8)}
    ck('degree8 monodromy equals AGL1(F8)',G==expected and len(G)==56)
    ck('degree8 full monodromy primitive',primitive(G,8))
    ck('Singer stabilizer acts transitively off origin',len({g[1] for g in G if g[0]==0})==7)
    # Real volume Jacobian on the torus and canonical compensating line.
    xx=s.symbols('x:6',nonzero=True);w=s.symbols('w');F=s.Matrix([s.prod(xx[j]**int(A[i,j]) for j in range(6)) for i in range(6)])
    jac=s.factor(F.jacobian(xx).det()*s.prod(xx)/s.prod(F))
    ck('degree8 torus logarithmic Jacobian exactly8',jac==8)
    ck('canonical-line compensation volume exactly1',jac*s.diff(w/8,w)==1)
    # Signed-permutation irreducibility mechanism: projections isolate every coordinate.
    for e in [1,2,4]:
        I=s.eye(e)
        for j in range(e):
            flip=s.eye(e);flip[j,j]=-1;project=(I-flip)*s.Rational(1,2)
            ck(f'signed rank{e} coordinate{j} projection over Z[1/2]',project==s.diag(*[int(i==j) for i in range(e)]))
    # General Sylow norm obstruction at a noncyclic Sylow, using A5 degree6.
    H=generated([(1,2,3,4,0),(0,4,3,2,1)],5);P=generated([(1,0,3,2,4),(2,3,0,1,4)],5)
    cosets=[];left=set(A5)
    while left:
        g=min(left);C=frozenset(mul(g,h) for h in H);cosets.append(C);left-=C
    acts=[tuple(next(i for i,C in enumerate(cosets) if mul(g,min(D)) in C) for D in cosets) for g in P]
    left=set(range(6));sizes=[]
    while left:
        j=min(left);orb={g[j] for g in acts};sizes.append(len(orb));left-=orb
    ck('general norm obstruction uses whole Sylow orbit degrees',sorted(sizes)==[2,2,2])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'free_base_controls':base_controls,'cyclotomic_degree8':{'matrix':[[int(v) for v in row] for row in A.tolist()],'determinant':8,'smith_factors':smith,'residue_polynomial':'z^3+z+1 over F2','actual_monodromy_order':56,'primitive':True},'general_sylow_norm_control':{'group':'A5 on D10 cosets','sylow':'V4','orbit_degrees':sizes},'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact induced finite-group actions and whole-torus freeness in three controls, cyclotomic identities/Smith quotient, full degree8 monodromy and volume, plus finite module/Sylow controls. No formalization of general Hilbert90/cohomology, ideal-class/stable-isomorphism arguments, canonical descent or the imported versality theorems. No independent whole-proof replay, external review or publication-priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
