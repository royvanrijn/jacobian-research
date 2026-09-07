#!/usr/bin/env python3
"""Interpolation verification of norm divisors, their node factor and genera."""
import argparse
from itertools import combinations
from math import gcd,lcm
from pathlib import Path
import json
import zipfile
import retrospective as r
import package_additive_branch_evidence as archive
import fresh_governing_panel as panel

OUTPUT=r.OUT/'rank_jump_additive_branch_geometry_verification_v1.json'


def evaluate(cs,t):
    result=0
    for c in reversed(cs):result=result*t+c
    return result


def determinant_norm(A,B,a,b,c):
    # Expand the six permutations in the cubic multiplication determinant.
    return c*(c-A*a)*(c-A*a)+(-B*a)*(-A*b-B*a)*a+(-B*b)*b*b-(-B*b)*(c-A*a)*a-(-B*a)*b*(c-A*a)-c*(-A*b-B*a)*b


def rational_sections(source,R):
    F=R.fraction_field();out=[]
    for section in source['sections']:
        coordinates=[]
        for key in ('x','y'):
            value=section[key]
            if isinstance(value,list):coordinates.append(F(R(value)))
            else:coordinates.append(F(R(value['numerator_coefficients_low_to_high']))/R(value['denominator_coefficients_low_to_high']))
        out.append(tuple(coordinates))
    return out


def norm_check(source,row,R):
    from sage.all import QQ,lcm as polynomial_lcm,gcd as polynomial_gcd
    A=R(source['A']);B=R(source['B']);sections=rational_sections(source,R)
    assert all(y*y==x**3+A*x+B for x,y in sections)
    assert all(c.denominator()==1 for c in A.list()+B.list())
    aa,bb=[list(map(int,p.list())) for p in (A,B)];count=0
    triples=list(combinations(range(len(sections)),3));assert len(triples)==len(row['norms'])
    for triple,saved in zip(triples,row['norms']):
        assert list(triple)==saved['generic_indices']
        x,y,z=[sections[i][0] for i in triple]
        raw=[-3+0*x,2*(x+y+z),(x+y-z)**2-4*x*y]
        common=polynomial_lcm([q.denominator() for q in raw]);coeff=[R(q*common) for q in raw]
        content=polynomial_gcd(coeff);coeff=[q//content for q in coeff]
        den=lcm(*(int(v.denominator()) for q in coeff for v in q.list()))
        ac,bc,cc=[list(map(int,(q*den).list())) for q in coeff]
        a,b,c=coeff
        # Every term of the symbolic determinant has degree at most this bound.
        terms=[(c,c-A*a,c-A*a),(-B*a,-A*b-B*a,a),(-B*b,b,b),
               (-B*b,c-A*a,a),(-B*a,b,c-A*a),(c,-A*b-B*a,b)]
        bound=max(sum(int(q.degree()) for q in term) for term in terms)
        N=R(saved['norm_ascending']);assert N.degree()==saved['degree']<=bound
        nc=list(map(int,N.list()));assert gcd(*nc)==1 and nc[-1]>0
        anchor=None
        for t in range(bound+1):
            value=determinant_norm(evaluate(aa,t),evaluate(bb,t),evaluate(ac,t),evaluate(bc,t),evaluate(cc,t))
            target=evaluate(nc,t)
            if anchor is None and target:anchor=(value,target);assert value
            if anchor is not None:assert value*anchor[1]==target*anchor[0]
            else:assert value==0
            count+=1
        assert anchor is not None
    return sections,count


def modular_checks(row,polys,D,degree_key,record_key,R):
    from sage.all import GF,PolynomialRing
    records=row[record_key];default=row['default_pair_witness_prime'];exceptions={tuple(x['pair']):x['prime'] for x in row['pair_witness_exceptions']}
    assert len(exceptions)==len(row['pair_witness_exceptions']) and not row['unresolved_pairs']
    primes={default}|set(exceptions.values())|{x['squarefree_witness_prime'] for x in records}|{x['discriminant_coprime_witness_prime'] for x in records}
    assert None not in primes
    reductions={p:[PolynomialRing(GF(p),'t')(q) for q in polys] for p in primes}
    for i,saved in enumerate(records):
        p=saved['squarefree_witness_prime'];q=reductions[p][i]
        assert q.degree()==saved[degree_key] and q.gcd(q.derivative())==1
        p=saved['discriminant_coprime_witness_prime'];q=reductions[p][i];dp=q.parent()(D)
        assert q.degree()==saved[degree_key] and dp.degree()==D.degree() and q.gcd(dp)==1
    count=0
    for i in range(len(polys)):
        for j in range(i):
            p=exceptions.get((j,i),default);a,b=reductions[p][j],reductions[p][i]
            assert a.degree()==records[j][degree_key] and b.degree()==records[i][degree_key] and a.gcd(b)==1
            count+=1
    assert count==row['pair_count']
    degrees=sorted(int(p.degree()) for p in polys);assert all(d%2==0 for d in degrees)
    genera=[{'classes':k,'genus_range':[1+2**k*(sum(ds)-4)//4 for ds in (degrees[:k],degrees[-k:])]} for k in range(1,5)]
    return count,genera


def verify():
    from sage.all import QQ,PolynomialRing
    manifest=r.read(archive.MANIFEST);assert r.digest(archive.BUNDLE.read_bytes())==manifest['bundle_sha256']
    assert manifest['packager_sha256']==r.digest(Path(archive.__file__).read_bytes())
    with zipfile.ZipFile(archive.BUNDLE) as z:raw={x['path']:z.read(x['member']) for x in manifest['members']}
    for entry in manifest['members']:assert len(raw[entry['path']])==entry['bytes'] and r.digest(raw[entry['path']])==entry['sha256']
    objects=[json.loads(raw[str(p.relative_to(r.ROOT))]) for p in archive.PATHS]
    for obj in objects+[r.read(archive.first.INPUT)]:
        for name,sha in obj['bindings'].items():assert r.digest(raw[name] if name in raw else (r.ROOT/name).read_bytes())==sha,name
    initial,complete,node=objects;assert node['status']=='PASS'
    R=PolynomialRing(QQ,'t');sources={x['family']:x for x in r.read(archive.first.INPUT)['families']}
    rows=[];evaluations=0;pairs=0;all_sections={}
    for row in complete['rows']:
        family=row['family'];source=sources[family];sections,count=norm_check(source,row,R);all_sections[family]=sections;evaluations+=count
        A=R(source['A']);B=R(source['B']);D=-4*A**3-27*B**2
        if family=='a1-fibration-01':
            tau=QQ(node['node_parameter']);e=QQ(node['node_x']);assert tau==-2
            assert A(tau)==-3*e*e and B(tau)==2*e**3 and e!=0
            assert D(tau)==D.derivative()(tau)==0 and D.derivative(2)(tau)!=0
            through=[i for i,(x,y) in enumerate(sections) if x(tau)==e and y(tau)==0]
            assert through==node['sections_through_node'] and len(through)==10
            wanted=set(combinations(through,3));polys=[]
            for old,new in zip(row['norms'],node['rows']):
                assert old['generic_indices']==new['generic_indices'];N=R(old['norm_ascending']);v=4 if tuple(old['generic_indices']) in wanted else 0
                assert new['node_multiplicity']==v
                assert all(N.derivative(j)(tau)==0 for j in range(v)) and N.derivative(v)(tau)!=0
                q=R(new['branch_ascending']);assert N==QQ(new['norm_decomposition_scalar'])*(R.gen()-tau)**v*q
                assert q.degree()==new['branch_degree'];polys.append(q)
            pc,genera=modular_checks(node,polys,D,'branch_degree','rows',R)
            assert genera==node['simultaneous_norm_cover_genus_ranges_first_four']
        else:
            assert row['status']=='PASS';polys=[R(x['norm_ascending']) for x in row['norms']]
            pc,genera=modular_checks(row,polys,D,'degree','norms',R)
            if family=='published-r17':assert [x['genus_range'][0] for x in genera]==row['summary']['simultaneous_norm_cover_genera_first_four']
            else:assert genera==row['summary']['simultaneous_norm_cover_genus_ranges_first_four']
        pairs+=pc;rows.append({'family':family,'norms':len(polys),'interpolation_evaluations':count,
            'coprime_branch_pairs':pc,'branch_degree_range':[int(min(q.degree() for q in polys)),int(max(q.degree() for q in polys))],
            'simultaneous_norm_cover_genus_ranges_first_four':genera})
        print(family,'PASS',len(polys),pc,flush=True)
    # Post-computation transport check only: these old parameter values never
    # enter generation of norms or selection of branch certificates.
    bridge=[]
    for label in r.read(panel.MANIFEST)['rows']:
        family=label['family']
        if label['parameter'] is None:continue
        t0=QQ(label['parameter']);source=sources[family];A=R(source['A'])(t0);B=R(source['B'])(t0)
        f,pts,_=panel.model_data(label['token']);u2=f[0]*A/(f[1]*B)
        assert u2>0 and u2.is_square() and f[1]==u2*u2*A and f[0]==u2**3*B
        xs=[u2*x(t0) for x,y in all_sections[family]];target=[x for x,y in pts]
        assert len(set(xs))==len(xs) and set(xs)==set(target)
        bridge.append({'token':label['token'],'family':family,'parameter':str(t0),'x_scale':str(u2),
                       'atlas_to_frozen_generic_indices':[target.index(x) for x in xs]})
    assert len(bridge)==13 and pairs==1079960
    return {'schema':'rank-jump.additive-branch-geometry-verification.v1','status':'PASS','rows':rows,
        'norm_polynomials':sum(x['norms'] for x in rows),'exact_interpolation_evaluations':evaluations,
        'pairwise_branch_coprimality_certificates':pairs,'node_triples_verified':120,'node_multiplicity':4,
        'frozen_specialization_bridge':bridge,
        'method':'Reconstruct primitive quadratics from rational-function sections; prove each norm identity by degree-bounded exact integer interpolation; replay every modular branch gcd; verify the node multiplicity by derivatives; check all thirteen fresh model and generic-x transports after arithmetic.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),archive.BUNDLE,archive.MANIFEST,archive.first.INPUT,panel.INPUT,panel.MANIFEST,Path(panel.__file__),Path(archive.__file__),Path(r.__file__))}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=verify()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS exact branch geometry and frozen generic transports')
