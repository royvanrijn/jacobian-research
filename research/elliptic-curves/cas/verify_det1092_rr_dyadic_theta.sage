#!/usr/bin/env sage-python
"""Independent Newton-sum resolvents and exhaustive local root-tree replay."""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,binomial
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_theta_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def power_sums(g,limit):
    n=g.degree();a=[g[n-i] for i in range(n+1)];p=[QQ(n)]
    for k in range(1,limit+1):
        value=-sum(a[j]*p[k-j] for j in range(1,min(k,n+1)))
        if k<=n:value-=k*a[k]
        p.append(value)
    return p
def from_sums(p,R):
    n=len(p)-1;a=[QQ(1)]
    for k in range(1,n+1):a.append(-sum(a[j]*p[k-j] for j in range(k))/k)
    return R(list(reversed(a)))
def resolvents(g):
    R=g.parent();p=power_sums(g,20)
    pair=[sum(binomial(k,j)*p[j]*p[k-j] for j in range(k+1))/2-2**(k-1)*p[k]
          for k in range(16)]
    triple=[]
    for k in range(21):
        total=sum(binomial(k,a)*binomial(k-a,b)*p[a]*p[b]*p[k-a-b]
                  for a in range(k+1) for b in range(k-a+1))
        diagonal=sum(binomial(k,a)*2**a*p[a]*p[k-a] for a in range(k+1))
        triple.append((total-3*diagonal+2*3**k*p[k])/6)
    assert pair[0]==15 and triple[0]==20
    theta=[sum(binomial(2*k,j)*2**j*(-p[1])**(2*k-j)*triple[j] for j in range(2*k+1))/2
           for k in range(11)]
    assert theta[0]==10
    return from_sums(pair,R),from_sums(theta,R)
def check_tree(F,tree):
    R=F.parent();x=R.gen();records={(ZZ(row['a']),row['depth']):row for row in tree['nodes']}
    assert len(records)==len(tree['nodes'])<=4096
    expected={(ZZ(0),0)};visited=set();roots=[]
    # Rebuild each node directly from the original polynomial and its ball,
    # rather than propagating the constructor's successive substitutions.
    for a,depth in sorted(records,key=lambda v:(v[1],v[0])):
        assert (a,depth) in expected and depth<=64 and 0<=a<2**depth
        visited.add((a,depth));row=records[a,depth]
        P=R(F(a+2**depth*x));content=min(c.valuation(2) for c in P if c);P/=2**content
        assert hashlib.sha256(str(P).encode()).hexdigest()==row['normalized_polynomial_sha256']
        assert all(c.denominator()==1 for c in P)
        assert [c['residue'] for c in row['children']]==[0,1]
        for child in row['children']:
            r=child['residue'];status=child['status']
            # Constant/linear coefficients of P(r+x) independently certify
            # the modulo-two root and simple Hensel condition.
            coefficients=R(P(r+x)).list();value=coefficients[0];derivative=coefficients[1]
            if value%2:assert status=='NO_ROOT_MOD_TWO'
            elif derivative%2:
                assert status=='UNIQUE_HENSEL_ROOT';roots.append((a+2**depth*r,depth+1))
            else:
                assert status=='SUBDIVIDE';expected.add((a+2**depth*r,depth+1))
    assert expected==visited
    assert sorted(roots)==sorted((ZZ(r['a']),r['depth']) for r in tree['root_balls'])
    assert len(roots)==tree['root_count'] and len(records)==tree['node_count']
    assert max(d for a,d in records)==tree['maximum_depth']
    # Leaf balls are pairwise disjoint, so the count has no multiplicity.
    for j,(a,d) in enumerate(roots):
        for b,e in roots[:j]:assert (a-b)%2**min(d,e)!=0
    return len(roots)
def sturm_real_count(f):
    seq=[f,f.derivative()]
    while seq[-1].degree()>0:seq.append(-(seq[-2]%seq[-1]))
    def variations(side):
        signs=[a.leading_coefficient().sign()*(side**a.degree()) for a in seq]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    return variations(-1)-variations(1)
def verify():
    proto=read(OUT/'protocol.json')
    for p,h in proto['inputs'].items():assert sha(ROOT/p)==h
    paths=[OUT/'protocol.json'];results=[];R=PolynomialRing(QQ,'x');x=R.gen()
    for i in range(10):
        path=OUT/('case-%02d.json'%i);d=read(path);paths.append(path)
        for p,h in d['inputs'].items():assert sha(ROOT/p)==h
        source=read(ROOT/d['source']);q=R(source['q']).monic();f=R(d['integral_monic_model'])
        s=d['T_equals_two_power_times_x'];assert f==q(2**s*x)/2**(6*s)
        assert all(c.valuation(2)>=0 for c in f if c) and f.is_monic()
        assert any(c.valuation(2)<6-j for j,c in enumerate(f.list()[:-1]) if c)
        g=R(d['proxy_polynomial']);N=d['coefficient_precision'];vd=d['discriminant_v2']
        assert 64<=N<=512 and N>2*vd and f.discriminant().valuation(2)==vd
        assert g.is_monic() and all(c.denominator()==1 for c in g)
        assert all((a-b).valuation(2)>=N for a,b in zip(f,g) if a!=b)
        assert g.discriminant().valuation(2)==vd
        pairs,theta=resolvents(g)
        assert pairs==R(d['pair_sum_resolvent']) and theta==R(d['theta_squared_difference_resolvent'])
        assert pairs.gcd(pairs.derivative())==theta.gcd(theta.derivative())==1
        counts={}
        for name,F in [('branch',g),('pair',pairs),('theta',theta)]:
            assert F.is_monic() and all(c.denominator()==1 for c in F)
            tp=OUT/('case-%02d-%s-tree.json'%(i,name));paths.append(tp);tree=read(tp)
            counts[name]=check_tree(F,tree)
            assert {k:tree[k] for k in ['root_count','node_count','maximum_depth']}==d['tree_summaries'][name]
        t=d['local_rational_2_torsion_dimension'];assert 2**t==1+counts['pair']
        divisible=bool(counts['branch'] or counts['theta']);kernel=int(not divisible)
        assert d['local_fake_kernel_dimension']==kernel
        assert d['inherited_canonical_minus_two_basepoint_is_locally_divisible_by_two']==divisible
        assert d['local_true_Kummer_dimension']==2+t and d['local_fake_Kummer_dimension']==2+t-kernel
        real_roots=sturm_real_count(q)
        pre=ART/'det1092_rr_dyadic_preflight_v1'/('case-%02d.json'%i)
        assert read(pre)['real_root_count_by_sturm']==real_roots;paths.append(pre)
        result={'case_index':i,'status':'PASS_INDEPENDENT_DYADIC_THETA_AND_TORSION',
                'root_counts':counts,'local_true_Kummer_dimension':2+t,
                'local_fake_Kummer_dimension':2+t-kernel,'D0_locally_divisible_by_two':divisible,
                'real_branch_points':real_roots,'real_Kummer_dimension':max(0,real_roots//2-1),
                'source':d['source'],'checker_sha256':sha(Path(__file__))}
        retain(OUT/('case-%02d-replay.json'%i),result);results.append(result)
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_TEN_DYADIC_THETA_AND_TORSION_CERTIFICATES',
            'cases':results,'proxy_equivalence_proof':'For each integral root alpha of f, v(fprime(alpha))<=D=v(disc(f)). Coefficient closeness N>2D gives a unique root of g within valuation>N-D>D of alpha by strong Hensel. Distinct f roots cannot have this distance, giving a Galois-equivariant bijection of the six branch sets.',
            'scope':'Exact Q2 torsion and local Kummer dimensions, and localization of the inherited rational D0 class. No full local image generators, global Selmer group or other Sha class is computed. At R only the dimension is certified here.',
            'limits':{**proto['limits'],'independent_panel_wall_seconds':25},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();retain(OUT/'replay.json',result)
    print(result['status'],flush=True)
