#!/usr/bin/env python3
"""Stdlib: reconstruct complete rational identities and literal integer Bezout."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib

def add(a,b,scale=1):
    c=dict(a)
    for m,v in b.items():c[m]=c.get(m,0)+scale*v
    return {m:v for m,v in c.items()if v}
def mul(a,b):
    c={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(x+y for x,y in zip(m,n));c[k]=c.get(k,0)+v*w
    return {m:v for m,v in c.items()if v}
def power(a,n):
    out={(0,)*len(next(iter(a))):1}
    for _ in range(n):out=mul(out,a)
    return out
def derivative(a):
    out={}
    for m,v in a.items():
        if m[0]:out[(m[0]-1,)+m[1:]]=m[0]*v
    return out
def flip(a):return {m:(-1)**m[0]*v for m,v in a.items()}
def product(items,dim):
    out={(0,)*dim:1}
    for a in items:out=mul(out,a)
    return out
def term(dim,powers,coefficient=1):
    m=[0]*dim
    for i,v in powers.items():m[i]=v
    return {tuple(m):coefficient}
def decode(terms):return {tuple(m):v for m,v in terms if v}

def reconstruct(row,variables):
    local=row['local'];pairs=local['paired_multiplicities'];ks=local['pair_pole_orders']
    dim=len(variables)+1;one=term(dim,{})
    params=[variables.index('a')+1]
    if len(pairs)==2:params.append(variables.index('b')+1)
    Ls=[];plus_factors=[];powers=[]
    for (m,n),k,j in zip(pairs,ks,params):
        L=add(one,term(dim,{0:2,j:2}),-1);Ls.append(L)
        exponent=max(m,n,k+1 if k else 0);powers.append(exponent)
        minus=add(one,term(dim,{0:1,j:1}),-1)
        plus=add(one,term(dim,{0:1,j:1}))
        plus_factors.extend([power(minus,exponent-m),power(plus,exponent-n)])
    N=add(term(dim,{0:1},2),term(dim,{0:3},2),-1)
    if 'q'in variables:
        N=mul(N,add(one,term(dim,{0:2,variables.index('q')+1:1})))
    first=mul(derivative(N),product([power(L,e-k)for L,e,k in zip(Ls,powers,ks)],dim))
    for i,k in enumerate(ks):
        if not k:continue
        other=product([power(L,e-t-int(j==i))for j,(L,e,t)in enumerate(zip(Ls,powers,ks))],dim)
        first=add(first,mul(mul(N,derivative(Ls[i])),other),-k)
    shifted=mul(power(add(one,term(dim,{0:1})),5),product(plus_factors,dim))
    numerator=add(first,add(shifted,flip(shifted)),-1)
    by_w={}
    for m,v in numerator.items():
        degree=m[0];tail=m[1:]
        by_w.setdefault(degree,{})[tail]=v
    assert all(k%2==0 for k in by_w)and 0 not in by_w
    return [by_w[k]for k in sorted(by_w,reverse=True)]

def verify(entry,row):
    inputs=[decode(x)for x in entry['inputs']]
    actual=reconstruct(row,entry['variables'])
    assert len(actual)==len(inputs)
    ratios=[]
    for coeff,p in zip(actual,inputs):
        assert set(coeff)==set(p)
        first=next(iter(p));ratio=Q(coeff[first],p[first])
        assert ratio and all(coeff[m]==ratio*v for m,v in p.items())
        ratios.append(str(ratio))
    total={}
    for multiplier,p in zip(entry['multipliers'],inputs):total=add(total,mul(decode(multiplier),p))
    assert entry['constant']!=0 and total=={(0,)*len(entry['variables']):entry['constant']}
    return {'index':entry['index'],'variables':entry['variables'],
            'coefficient_scales':ratios,'literal_constant':entry['constant'],
            'multiplier_terms':sum(len(p)for p in entry['multipliers'])}


def univariate(values):return {(i,):v for i,v in enumerate(values)if v}
def coeff_w(poly,k):return {m[1:]:v for m,v in poly.items()if m[0]==k}
def scale(poly,v):return {m:v*w for m,w in poly.items()if v*w}
def reflection_check(case):
    # Ordered variables w,r; normalize s=1 by its nonzero scalar.
    one=term(2,{});w=term(2,{0:1});r=term(2,{1:1})
    r2=power(r,2);w2=power(w,2)
    wr=add(w2,r2,-1);ws=add(w2,one,-1)
    k,l=case['pole_orders'];A=case['A'];B=case['B']
    D=add(add(mul(wr,ws),mul(w2,ws),-2*k),mul(w2,wr),-2*l)
    P=mul(power(add(w,r,-1),A),power(add(w,one,-1),B))
    even={m:v for m,v in P.items()if m[0]%2==0}
    constant=coeff_w(even,0)
    E=add(mul({(0,)+m:v for m,v in constant.items()},D),mul(r2,even),-1)
    assert coeff_w(E,0)=={}
    q=univariate(case['quadratic']);high=univariate(case['high'])
    if A==1:
        assert coeff_w(E,4)==mul({(1,):1},high)
        assert coeff_w(E,2)==mul({(1,):-1},q)
    else:
        assert coeff_w(E,4)==mul({(2,):2},high)
        assert coeff_w(E,2)==mul({(2,):-2},q)
    total=add(mul(univariate(case['q_multiplier']),q),
              mul(univariate(case['high_multiplier']),high))
    assert total=={(0,):case['constant']}and case['constant']
    return {'pole_orders':[k,l],'literal_constant':case['constant']}

def integer_branches():
    pole=[(g,c,r)for g in range(1,17)for c in range(2,20)for r in range(1,20)
          if g*(3*c+2*r-3)==16]
    assert pole==[(1,3,5),(1,5,2),(2,3,1)]
    moving=[(g,c,3*c-4)for g in range(1,17)for c in range(2,20)if g*(3*c-5)==16]
    assert moving==[(1,7,17),(4,3,5),(16,2,2)]
    fixed=[(g,1+16//g)for g in range(1,17)if 16%g==0]
    assert fixed==[(1,17),(2,9),(4,5),(8,3),(16,2)]
    edge=[]
    for g in range(1,19):
        if 18%g:continue
        c=1+18//g
        for e in range(0,c-1):
            if(c-e)%2:continue
            for H in range(4,c+1,2):
                if 5*e<3*c+6*H:continue
                s=(3*c-e-2)//2
                for n in range(s+1,2*(c-H)+1):
                    if n%2 or 6*c<5*n+3*H:continue
                    edge.append((g,c,e,H,s,n))
    assert edge==[(1,19,17,4,19,20)]
    return {'D_pole_cases':pole,'D_moving_zero_cases':moving,
            'D_fixed_zero_cases':fixed,'C_reflection_edge':edge}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    parser.add_argument('--source-root',type=Path)
    args=parser.parse_args();here=Path(__file__).resolve().parent
    proof=here/'ODD_CUBIC_DEGREE_NINETEEN_SAFETY.md'
    data_path=here/'odd_cubic_degree_nineteen_safety_certificates.json'
    ledger_path=here/'odd_cubic_degree_nineteen_safety_source_ledger.json'
    data=json.loads(data_path.read_text());ledger=json.loads(ledger_path.read_text())
    canonical='ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md'
    source_root=args.source_root or(here if here.name=='extended-geometry'else Path.cwd()/'research/foundations/extended-geometry')
    pins=[]
    for resource in ledger['retained_resources']:
        p=(here if resource.get('location')=='package'else source_root)/resource['retained_path']
        assert sha(p)==resource['sha256'];pins.append({'path':resource['retained_path'],'sha256':sha(p)})
    moving=[verify(entry,entry['normal_form'])for entry in data['moving_cases']]
    assert len(moving)==10
    reflection=[reflection_check(entry)for entry in data['reflection_cases']]
    assert len(reflection)==3
    mutant=dict(data['moving_cases'][0]);mutant['constant']+=1
    try:verify(mutant,mutant['normal_form'])
    except AssertionError:pass
    else:raise AssertionError('Corrupt literal constant accepted')
    out={'passed':True,'scope':'Exact global integer reductions, complete product-rule/reflection numerators, and complete moving-case rational-purity numerators with106 literal integer Bezout multiplier terms. Universal exclusion follows from the separately written exhaustive field/normal-form proof, not a moment prefix.',
         'integer_branches':integer_branches(),'moving_cases':moving,
         'reflection_cases':reflection,'negative_controls':{'corrupt_constant_rejected':True},
         'pins':pins+[{'path':Path(__file__).name,'sha256':sha(Path(__file__))},
                      {'path':ledger_path.name,'sha256':sha(ledger_path)}]}
    out.update({'schema_version':1,'status':'passed','independent_replay':False,'external_review':False,'formal_verification':False})
    output=args.output or here/'odd_cubic_degree_nineteen_safety_verification.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'moving_full_identities':len(moving),
        'literal_moving_multiplier_terms':sum(r['multiplier_terms']for r in moving),
        'reflection_full_identities':len(reflection),'output':str(output)}))
if __name__=='__main__':main()
