#!/usr/bin/env python3
"""Small exact controls for the written structural proofs; not a whole-proof replay."""
from pathlib import Path
from itertools import product, permutations
import hashlib, importlib.util, json, platform, time
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('retained_permutation_helpers',HERE.parent/'universal-keller-homogeneous-receivers/verify.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
mul,inv,generated,cosets=old.mul,old.inv,old.generated,old.cosets
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def derived(G,n):
    return generated([mul(mul(a,b),mul(inv(a),inv(b))) for a in G for b in G],n)
def transitive(G,n):return len({a[0] for a in G})==n

def main():
    started=time.monotonic()
    R=s.diag(1,1,s.Matrix([[-1,-1],[1,0]]))
    S=s.diag(-1,1,s.Matrix([[0,1],[1,0]]))
    ck('quartic linear holonomy S3 relations',R**3==s.eye(4) and S**2==s.eye(4) and S*R*S==R**2)
    ck('quartic faithful linear holonomy',len({tuple(R**i*S**j) for i in range(3) for j in range(2)})==6)
    variables=s.symbols('m:16');M=s.Matrix(4,4,variables)
    solution=s.linsolve(list(M*R-R*M)+list(M*S-S*M),variables)
    expected=(variables[0],0,0,0,0,variables[5],0,0,0,0,variables[15],0,0,0,0,variables[15])
    ck('exact rational centralizer is sign plus trivial plus scalar standard',solution==s.FiniteSet(expected))
    # Affine centralizer over sixth roots of unity, where theta=2 and -1=3.
    rr=(2,0,0,0);ss=(0,3,0,0)
    def mv(A,t):return tuple(sum(int(A[i,j])*t[j] for j in range(4))%6 for i in range(4))
    actual=set()
    for signs in product([1,-1],repeat=3):
        E=s.diag(signs[0],signs[1],signs[2],signs[2])
        Er,Es=mv(E,rr),mv(E,ss)
        for t in product(range(6),repeat=4):
            Rt,St=mv(R,t),mv(S,t)
            if all((t[i]+Er[i]-rr[i]-Rt[i])%6==0 and (t[i]+Es[i]-ss[i]-St[i])%6==0 for i in range(4)):
                actual.add((signs,t))
    predicted={((1,d,e),(x,z,a,a)) for d,e,x,z,a in product([1,-1],[1,-1],[0,3],range(6),[0,2,4])}
    ck('complete finite affine centralizer agrees with formula',actual==predicted)
    ck('finite affine centralizer order',len(actual)==144)
    def amul(g,h):
        a,t=g;b,u=h;eps=(a[0],a[1],a[2],a[2])
        return (tuple(x*y for x,y in zip(a,b)),tuple((t[i]+eps[i]*u[i])%6 for i in range(4)))
    def ainv(g):
        a,t=g;eps=(a[0],a[1],a[2],a[2]);return a,tuple(-eps[i]*t[i]%6 for i in range(4))
    comm={amul(amul(a,b),amul(ainv(a),ainv(b))) for a in actual for b in actual}
    ck('finite centralizer commutator is abelian order9',len(comm)==9 and all(amul(a,b)==amul(b,a) for a in comm for b in comm))
    # The natural symmetric second derived subgroup is transitive.
    derived_controls=[]
    for n in [4,5]:
        G=frozenset(permutations(range(n)));D=derived(G,n);DD=derived(D,n)
        ck(f'S{n} second derived transitive',transitive(DD,n))
        ck(f'S{n} second derived nontrivial',len(DD)==(4 if n==4 else 60))
        derived_controls.append({'rank':n,'derived_order':len(D),'second_derived_order':len(DD),'second_derived_orbit':len({a[0] for a in DD})})
    # Whole Sylow fixed loci genuinely differ from single-element fixed loci.
    A5=generated([(1,2,0,3,4),(1,2,3,4,0)],5)
    H=generated([(1,2,3,4,0),(0,4,3,2,1)],5)
    P=generated([(1,0,3,2,4),(2,3,0,1,4)],5)
    C=cosets(A5,H)
    action=[tuple(next(j for j,c in enumerate(C) if mul(a,next(iter(b))) in c) for b in C) for a in P]
    ck('A5 degree6 coset and Sylow orders',(len(A5),len(H),len(P),len(C))==(60,10,4,6))
    ck('A5 coset Sylow has no common fixed label',not any(all(a[i]==i for a in action) for i in range(6)))
    individual=sorted(sum(a[i]==i for i in range(6)) for a in action)
    ck('A5 coset each involution fixes two labels',individual==[2,2,2,6])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'finite_centralizer':{'translation_modulus':6,'order':len(actual),'commutator_order':len(comm),'second_commutator_order':1},'symmetric_derived_controls':derived_controls,'sylow_coset_control':{'group':'A5','subgroup_order':10,'rank':6,'sylow_order':4,'individual_fixed_counts':individual,'common_fixed_labels':0},'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact matrix, finite affine centralizer and permutation controls for the written all-degree proofs. No formalization of topology, Smith theory, normalization, lifting, cylinder rigidity or JC reductions; no independent whole-proof replay, external review or publication-priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'finite_centralizer_order':len(actual),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
