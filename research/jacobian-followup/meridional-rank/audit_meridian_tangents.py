#!/usr/bin/env python3
"""Exact first-order KCH deformations at the retained H3 points.

The frozen campaign uses seed propagation. A separate verifier evaluates
the oriented relator words directly with dual numbers and checks minors.
"""
import hashlib
from itertools import permutations
import json
from pathlib import Path
import time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from parabolic_lifts import oriented_diagram
from quadratic_integer_lifts import ZERO, ONE, IDENTITY, add, neg, mul, matmul, h3_model

HERE=Path(__file__).resolve().parent


def ma(*xs):
    return tuple(tuple_sum(a) for a in zip(*xs))


def tuple_sum(xs):
    out=ZERO
    for x in xs:
        out=add(out,x)
    return out


def mn(a):
    return tuple(neg(x) for x in a)


def product(*xs):
    out=IDENTITY
    for x in xs:
        out=matmul(out,x)
    return out


def conjugation_derivative(o,a,b,do,da,positive):
    if positive:
        return ma(product(do,a,o),product(o,da,o),mn(product(b,do,o)))
    return ma(mn(product(o,do,b)),product(o,da,o),product(o,a,do))


def base_matrices(arc_images, elements, model):
    images=[tuple(elements[i]) for i in arc_images]
    def pm(a,b):
        return tuple(a[b[i]] for i in range(5))
    identity=tuple(range(5))
    def pi(a):
        return tuple(a.index(i) for i in range(5))
    deck=[c for c in permutations(range(5))
          if all(pm(pm(c,a),pi(c))==pi(a) for a in images)]
    assert len(deck)==1
    c=deck[0]
    assert pm(c,c)==identity and sum(c[i]!=i for i in range(5))==4
    keys=[(pm(a,c),1) for a in images]
    assert all(k in model for k in keys)
    return c,[model[k] for k in keys]


def seed_derivatives(matrices,crossings,seeds,trace,variable):
    zero=(ZERO,)*9
    derivatives=[None]*len(matrices)
    derivatives[seeds[0]]=ma(IDENTITY,mn(matrices[seeds[0]])) if variable==0 else zero
    for j,seed in enumerate(seeds[1:]):
        derivatives[seed]=tuple(ONE if variable==1+9*j+k else ZERO for k in range(9))
    for o,a,b in trace:
        assert derivatives[o] is not None and derivatives[a] is not None and derivatives[b] is None
        matches=[row for row in crossings if row[0]==o and set(row[1:3])=={a,b}]
        assert matches and len(set(row[3]^(a!=row[1]) for row in matches))==1
        sign=matches[0][3]^(a!=matches[0][1])
        derivatives[b]=conjugation_derivative(matrices[o],matrices[a],matrices[b],
                                              derivatives[o],derivatives[a],not sign)
    assert all(d is not None for d in derivatives)
    return derivatives


def seed_matrix(matrices,crossings,seeds,trace):
    columns=[]
    for variable in range(1+9*(len(seeds)-1)):
        ds=seed_derivatives(matrices,crossings,seeds,trace,variable)
        col=[]
        for o,a,b,sign in crossings:
            assert product(matrices[o],matrices[a],matrices[o])==matrices[b]
            db=conjugation_derivative(matrices[o],matrices[a],matrices[b],ds[o],ds[a],not sign)
            col.extend(ma(ds[b],mn(db)))
        for seed in seeds[1:]:
            plus=ma(IDENTITY,matrices[seed])
            col.extend(product(plus,ds[seed],plus))
        columns.append(col)
    return [list(r) for r in zip(*columns)]


def modular_rank(rows,prime,tau):
    a=[[(x+y*tau)%prime for x,y in row] for row in rows]
    n=len(a[0]);k=0
    for j in range(n):
        i=next((i for i in range(k,len(a)) if a[i][j]),None)
        if i is None:
            continue
        a[k],a[i]=a[i],a[k]
        c=pow(a[k][j],-1,prime)
        a[k]=[x*c%prime for x in a[k]]
        for i in range(k+1,len(a)):
            c=a[i][j]
            if c:
                a[i]=[(x-c*y)%prime for x,y in zip(a[i],a[k])]
        k+=1
    return k


def exact_certificate(rows):
    tau=(1+sp.sqrt(5))/2
    field=sp.QQ.algebraic_field(tau)
    converted=[[field.from_sympy(x+y*tau) for x,y in row] for row in rows]
    dm=DomainMatrix(converted,(len(rows),len(rows[0])),field)
    ns=dm.nullspace().to_list()
    def encode(x):
        coeff=x.to_list()
        if not coeff:
            a,b=field.dom.zero,field.dom.zero
        elif len(coeff)==1:
            a,b=coeff[0],field.dom.zero
        else:
            assert len(coeff)==2
            b,a=coeff
        return [[int(a.numerator),int(a.denominator)],[int(b.numerator),int(b.denominator)]]
    basis=[[encode(x) for x in row] for row in ns]
    controls=[{'prime':p,'tau':t,'rank':modular_rank(rows,p,t)} for p in [11,19,29,31]
              for t in range(p) if (t*t-t-1)%p==0]
    rank=len(rows[0])-len(basis)
    assert max(c['rank'] for c in controls)==rank
    return {'rows':len(rows),'columns':len(rows[0]),'exact_rank':rank,'nullity':len(basis),
            'basis':basis,'modular_rank_certificates':controls,
            'eigenvalue_derivative_available':any(row[0]!=[[0,1],[0,1]] for row in basis)}


def dual_word(matrices,derivatives,word):
    base,derivative=IDENTITY,(ZERO,)*9
    for letter in word:
        i=abs(letter)-1
        m,d=matrices[i],derivatives[i]
        if letter<0:
            # All base matrices in this helper's controls need not be involutions.
            m=inverse_matrix(m)
            d=mn(product(m,d,m))
        derivative=ma(product(derivative,m),product(base,d))
        base=product(base,m)
    return base,derivative


def inverse_matrix(a):
    # Cofactors over Z[tau], for controls with determinant +/-1.
    co=[]
    for i in range(3):
        for j in range(3):
            rr=[r for r in range(3) if r!=i]; cc=[c for c in range(3) if c!=j]
            d=add(mul(a[3*rr[0]+cc[0]],a[3*rr[1]+cc[1]]),neg(mul(a[3*rr[0]+cc[1]],a[3*rr[1]+cc[0]])))
            co.append(d if (i+j)%2==0 else neg(d))
    det=tuple_sum(mul(a[j],co[j]) for j in range(3))
    assert det in [ONE,neg(ONE)]
    return tuple(mul(det,co[3*j+i]) for i in range(3) for j in range(3))


def controls():
    # The exact T35 family from the preceding turn, specialized at s=-1.
    f=tuple((v,0) for v in [0,0,1,1,0,0,0,1,0])
    cc=[(0,0),(0,0),(-1,0),(1,0),(0,0),(-1,1),(0,0),(1,0),(-1,1)]
    dd=mn(cc)
    x=mn(f); y=mn(product(dd,dd,dd))
    m,_=dual_word([x,y],[(ZERO,)*9]*2,[-1,2,2])
    assert product(m,m)==IDENTITY
    columns=[]
    for variable in range(19):
        ds=[tuple(ONE if variable==1+9*j+k else ZERO for k in range(9)) for j in range(2)]
        c,dc=dual_word([x,y],ds,[1,1,1,-2,-2,-2,-2,-2])
        assert c==IDENTITY
        _,dm=dual_word([x,y],ds,[-1,2,2])
        target=ma(IDENTITY,mn(m)) if variable==0 else (ZERO,)*9
        columns.append(list(dc)+list(ma(dm,mn(target))))
    matrix=[list(r) for r in zip(*columns)]
    free=exact_certificate(matrix)
    assert free['eigenvalue_derivative_available']
    # Adding m^2=I linearizes, after normalization, to s=0 in characteristic0.
    extra=[[ONE]+[ZERO]*18]
    fixed=exact_certificate(matrix+extra)
    assert not fixed['eigenvalue_derivative_available']
    return {'torus_matrices':[x,y],'meridian':m,'varying_T35':free,'involutive_T35':fixed}


def main():
    plan=json.loads((HERE/'meridian_tangent_plan.json').read_text())
    assert plan['frozen_before_execution']
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases'][:16]
    maps=json.loads((HERE/'joint_a5_results.json').read_text())['results'][:16]
    catalogue=json.loads((HERE/'joint_a5_catalogue.json').read_text())
    model,_,_=h3_model()
    paths=['audit_meridian_tangents.py','meridian_tangent_plan.json','joint_a5_panel.json',
           'joint_a5_results.json','joint_a5_catalogue.json','quadratic_integer_lifts.py','parabolic_lifts.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    output=HERE/'meridian_tangent_results.json'
    report=json.loads(output.read_text()) if output.exists() else {'status':'RUNNING','results':[], 'input_sha256':hashes}
    assert report['input_sha256']==hashes
    started=time.monotonic()
    if 'controls' not in report:
        report['controls']=controls()
        output.write_text(json.dumps(report,indent=2)+'\n')
    done={(r['name'],r['orbit_index']) for r in report['results']}
    for case,source in zip(panel,maps):
        assert case['name']==source['name']
        crossings=oriented_diagram(case['pd'])
        for mapping in source['epimorphism_representatives']:
            key=source['name'],mapping['orbit_index']
            if key in done:
                continue
            assert time.monotonic()-started<plan['limits']['seconds_total']
            stamp=time.monotonic()
            deck,matrices=base_matrices(mapping['arc_images'],catalogue['elements'],model)
            seeds,trace=source['seeds'],source['seed_trace']
            assert len(seeds)==4
            matrix=seed_matrix(matrices,crossings,seeds,trace)
            cert=exact_certificate(matrix)
            assert cert['rows']<=144 and cert['columns']==28
            row={'name':key[0],'orbit_index':key[1],'source_group':mapping['source_group'],
                 'deck':deck,'matrices':matrices,'seeds':seeds,'trace':trace,
                 'matrix':matrix,'certificate':cert,'seconds':time.monotonic()-stamp}
            report['results'].append(row)
            output.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'name':key[0],'orbit_index':key[1],**{k:cert[k] for k in ['exact_rank','nullity','eigenvalue_derivative_available']},'seconds':row['seconds']}),flush=True)
            assert row['seconds']<plan['limits']['seconds_per_case']
    assert len(report['results'])==16
    report['status']='COMPLETE'
    output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
