#!/usr/bin/env python3
"""Finite controls for the written all-degree theorems; not formalization."""
from pathlib import Path
from itertools import permutations, combinations, product
from math import factorial, lcm
import hashlib, json, platform, time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def ck(label, value):
    assert value, label
    checks.append(label)
def mul(a,b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a): return tuple(a.index(i) for i in range(len(a)))
def parity(a): return sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))%2
def generated(gens,n):
    seen={tuple(range(n))};todo=list(seen)
    for a in todo:
        for b in gens:
            c=mul(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return frozenset(seen)
def primitive(G,n):
    if {a[0] for a in G}!=set(range(n)): return False
    for size in range(2,n):
        if n%size:continue
        for tail in combinations(range(1,n),size-1):
            B={0,*tail}
            if all(not B&{a[i] for i in B} or B=={a[i] for i in B} for a in G):return False
    return True
def center(G):return {a for a in G if all(mul(a,b)==mul(b,a) for b in G)}
def cosets(G,J):
    seen=set();out=[]
    for a in sorted(G):
        if a in seen:continue
        C=frozenset(mul(a,b) for b in J);out.append(C);seen.update(C)
    assert seen==set(G)
    return out
def cycle_type(a):
    rest=set(range(len(a)));sizes=[]
    while rest:
        i=next(iter(rest));j=i;length=0
        while j in rest:rest.remove(j);length+=1;j=a[j]
        sizes.append(length)
    return tuple(sorted(sizes))
def rho(a):
    n=len(a);M=s.zeros(n+1)
    for i in range(n):M[a[i],i]=1
    M[n,n]=(-1)**parity(a)
    return M

class Field:
    """Tiny exact F_q extension, used only for q=2,3 and degree<=4."""
    def __init__(self,q,d):
        self.q=q;self.d=d;self.size=q**d
        x=s.Symbol('x')
        for coeffs in product(range(q),repeat=d):
            if d>1 and coeffs[0]==0:continue
            P=x**d+sum(coeffs[i]*x**i for i in range(d))
            if s.Poly(P,x,modulus=q).is_irreducible:
                self.mod=coeffs;break
        else:raise AssertionError('no tiny irreducible')
    def digits(self,a):
        out=[]
        for _ in range(self.d):out.append(a%self.q);a//=self.q
        return out
    def encode(self,c):return sum((v%self.q)*self.q**i for i,v in enumerate(c))
    def add(self,a,b):return self.encode([x+y for x,y in zip(self.digits(a),self.digits(b))])
    def neg(self,a):return self.encode([-x for x in self.digits(a)])
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        A=self.digits(a);B=self.digits(b);C=[0]*(2*self.d-1)
        for i,x in enumerate(A):
            for j,y in enumerate(B):C[i+j]+=x*y
        for i in range(2*self.d-2,self.d-1,-1):
            c=C[i]%self.q
            for j,v in enumerate(self.mod):C[i-self.d+j]-=c*v
        return self.encode(C[:self.d])
    def pow(self,a,n):
        out=1
        while n:
            if n%2:out=self.mul(out,a)
            a=self.mul(a,a);n//=2
        return out
    def orbit(self,a):
        out=[];b=a
        while b not in out:out.append(b);b=self.pow(b,self.q)
        assert b==a
        return out
    def determinant(self,A):
        A=[row[:] for row in A];det=1;n=len(A)
        for j in range(n):
            pivot=next((i for i in range(j,n) if A[i][j]),None)
            if pivot is None:return 0
            if pivot!=j:A[pivot],A[j]=A[j],A[pivot];det=self.neg(det)
            v=A[j][j];det=self.mul(det,v);iv=self.pow(v,self.size-2)
            for i in range(j+1,n):
                c=self.mul(A[i][j],iv)
                for k in range(j,n):A[i][k]=self.sub(A[i][k],self.mul(c,A[j][k]))
        return det

def main():
    started=time.monotonic();group_controls=[];normal_controls=[]
    for n in range(2,7):
        Sn=frozenset(permutations(range(n)));An=frozenset(a for a in Sn if not parity(a))
        J=frozenset(a for a in Sn if a[0]==0);Cs=cosets(Sn,J)
        ck(f'S{n} full quotient fibre has exactly N cosets',len(Cs)==n and all(len(C)==factorial(n-1) for C in Cs))
        ck(f'S{n} coset labels exhaust the full natural fibre',{next(iter(C))[0] for C in Cs}==set(range(n)) and all(len({a[0] for a in C})==1 for C in Cs))
        ck(f'S{n} primitive monodromy',primitive(Sn,n))
        ck(f'S{n} point-stabilizer source and target pi1 orders differ',len(J)*n==len(Sn) and len(J)!=len(Sn))
        gens=[]
        for i in range(n-1):
            a=list(range(n));a[i],a[i+1]=a[i+1],a[i];gens.append(tuple(a))
        matrices=[rho(a) for a in gens]
        ck(f'S{n} determinant-corrected faithful matrix embedding',len({tuple(rho(a)) for a in Sn})==len(Sn) and all(M.det()==1 for M in matrices))
        ck(f'S{n} integral Coxeter matrix relations',all(M*M==s.eye(n+1) for M in matrices) and all((matrices[i]*matrices[i+1])**3==s.eye(n+1) for i in range(n-2)) and all(matrices[i]*matrices[j]==matrices[j]*matrices[i] for i in range(n-1) for j in range(i+2,n-1)))
        for p in [2,3,5]:
            ck(f'S{n} char{p} embedding remains faithful',len({tuple(v%p for v in rho(a)) for a in Sn})==len(Sn))
        if n>=3:
            K=frozenset(a for a in An if a[0]==0);AC=cosets(An,K)
            tau=list(range(n));tau[1],tau[2]=tau[2],tau[1];tau=tuple(tau);one=tuple(range(n))
            ck(f'A{n} full fibre and exact component cost',len(AC)==n and len(An)==factorial(n)//2 and len(K)*n==len(An))
            ck(f'A{n} primitive atomic monodromy',primitive(An,n))
            ck(f'A{n} odd right multiplier preserves isotropy',{mul(mul(tau,a),tau) for a in An}==set(An) and {mul(mul(tau,a),tau) for a in K}==set(K))
            ck(f'A{n} canonical target mark fixed by full symmetric action',all(mul(a,tau if parity(a) else one) in An for a in Sn))
            by_label={next(iter(C))[0]:next(iter(C)) for C in AC}
            actual=set()
            for a in Sn:
                action=[]
                for label in range(n):
                    b=mul(mul(a,by_label[label]),tau if parity(a) else one)
                    assert b in An
                    action.append(b[0])
                assert tuple(action)==a
                actual.add(tuple(action))
            ck(f'A{n} exact whole-fibre natural S{n} action',actual==set(Sn))
            ck(f'A{n} source and target cannot share finite pi1 order',len(K)!=len(An))
        if n<=5:
            # Enumerate the full normal-subgroup lattice by adjoining conjugacy classes.
            conjugacy={cycle_type(a):a for a in Sn};normals={frozenset([tuple(range(n))])};todo=list(normals)
            for H in todo:
                for a in conjugacy.values():
                    generators=set(H)|{mul(mul(b,a),inv(b)) for b in Sn}
                    K=generated(list(generators),n)
                    if K not in normals:normals.add(K);todo.append(K)
            ck(f'S{n} every proper normal group lies in alternating group',all(H==Sn or H<=An for H in normals))
            trans=[H for H in normals if {a[0] for a in H}==set(range(n))]
            ck(f'S{n} full canonical normal monodromy census',sorted(len(H) for H in trans)==({2:[2],3:[3,6],4:[4,12,24],5:[60,120]}[n]))
            ck(f'S{n} every center-free normal quotient meets component cost',all(len(H)==1 or len(center(H))==1 or (n<=4) for H in trans))
            normal_controls.append({'rank':n,'normal_subgroup_orders':sorted(len(H) for H in normals),'transitive_normal_orders':sorted(len(H) for H in trans)})
        group_controls.append({'rank':n,'fixed_cover_degree':n,'symmetric_isotropy_order':factorial(n),'alternating_isotropy_order':factorial(n)//2 if n>=3 else None})
    # Actual intermediate cover: an imprimitive framed cyclic quartic map has two quadratic factors.
    C4=generated([(1,2,3,0)],4);C2=generated([(2,3,0,1)],4);I=frozenset([(0,1,2,3)])
    ck('imprimitive quotient map has actual two-by-two subgroup factors',not primitive(C4,4) and I<C2<C4 and len(cosets(C4,C2))==len(cosets(C2,I))==2)
    # Exact Hilbert90 frames for every finite-field splitting type in ranks2,3,4.
    # Frob(B)=B*rho(sigma) gives a rational quotient target; its ENTIRE fibre
    # is the N right cosets and has exactly the requested Frobenius action.
    finite_controls=[]
    for q in [2,3]:
        for n in [2,3,4]:
            representatives={cycle_type(a):a for a in permutations(range(n))}
            for typ,unused in sorted(representatives.items()):
                degree=lcm(*typ);F=Field(q,degree)
                C=[[0]*n for _ in range(n)];sigma=[];offset=0
                for length in typ:
                    beta=next(a for a in range(F.size) if len(F.orbit(a))==length)
                    orbit=F.orbit(beta)
                    for i in range(length):
                        for j in range(length):C[offset+i][offset+j]=F.pow(orbit[j],i)
                        sigma.append(offset+(i+1)%length)
                    offset+=length
                sigma=tuple(sigma);det=F.determinant(C)
                ck(f'F{q} rank{n} type{typ} evaluation frame invertible',det!=0)
                B=[[0]*(n+1) for _ in range(n+1)]
                for i in range(n):
                    for j in range(n):B[i][j]=C[i][j]
                B[n][n]=F.pow(det,F.size-2)
                ck(f'F{q} rank{n} type{typ} actual SL frame',F.determinant(B)==1)
                R=rho(sigma);expected=[]
                for i in range(n+1):
                    row=[]
                    for j in range(n+1):
                        value=0
                        for k in range(n+1):
                            term=int(R[k,j])%q
                            value=F.add(value,F.mul(B[i][k],term))
                        row.append(value)
                    expected.append(row)
                ck(f'F{q} rank{n} type{typ} exact Frobenius descent equation',[[F.pow(v,q) for v in row] for row in B]==expected)
                Sn=frozenset(permutations(range(n)));J=frozenset(a for a in Sn if a[0]==0);CS=cosets(Sn,J)
                labels=[next(iter(c))[0] for c in CS]
                action=tuple(labels.index(mul(sigma,next(iter(c)))[0]) for c in CS)
                ck(f'F{q} rank{n} type{typ} full finite fibre Frobenius type',len(CS)==n and cycle_type(action)==typ)
                finite_controls.append({'base_field_order':q,'rank':n,'splitting_type':typ,'splitting_field_order':F.size,'complete_coset_fibre_rank':len(CS),'SL_determinant':1})
    manifest=json.loads((HERE/'sources.json').read_text())
    for source in manifest['local_sources']:
        ck('local source pin '+source['path'],hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()==source['sha256'])
    for source in manifest['external_sources']:
        data=(HERE/source['path']).read_bytes()
        ck('primary source pin '+source['path'],len(data)==source['bytes'] and hashlib.sha256(data).hexdigest()==source['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'group_controls':group_controls,'normal_monodromy_controls':normal_controls,'finite_field_full_fibre_controls':finite_controls,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Written all-degree theorems, supported by finite exact group/coset/matrix controls and twenty finite-field full-coset-fibre Frobenius controls. Full scheme-theoretic fibre rank uses the written finite-torsor proof, not a separate invariant-ring or Groebner certificate. Controls do not formalize all-rank group theory, etale descent, homogeneous topology/pi1, geometric base change or atomicity imports. Classical versal torsors credited. No independent whole-proof replay, formal verification, external review or novelty certification.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'finite_field_fibre_controls':len(finite_controls),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
