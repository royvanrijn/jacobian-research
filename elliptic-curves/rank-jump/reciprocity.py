#!/usr/bin/env python3
"""Explain correlations of the fixed-field local cuts by rational reciprocity."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc

PROTOCOL=Path(__file__).resolve().parent/'RECIPROCITY_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_reciprocity_v1.json'


def valuation_unit(value,p):
    value=r.F(value)
    if not value:raise ValueError('zero has no Hilbert unit')
    a,b=value.numerator,value.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v,r.F(a,b)


def hilbert(a,b,p):
    a,b=r.F(a),r.F(b)
    if not a or not b:raise ValueError('zero Hilbert argument')
    if p=='infinity':return -1 if a<0 and b<0 else 1
    v,x=valuation_unit(a,p);w,y=valuation_unit(b,p)
    if p==2:
        x=r.mod(x,8);y=r.mod(y,8)
        bit=((x-1)*(y-1)//4+v*(y*y-1)//8+w*(x*x-1)//8)%2
    else:
        bit=(v*w*((p-1)//2)+w*int(pow(r.mod(x,p),(p-1)//2,p)==p-1)+v*int(pow(r.mod(y,p),(p-1)//2,p)==p-1))%2
    return -1 if bit else 1


def build(check=False):
    inp=r.read(lc.INPUT);prior=r.read(lc.OUTPUT)
    model=inp['anchor']['short_model_ainvariants'];pts=inp['anchor']['known_points_on_short_model'];A,B=map(r.F,model[3:])
    assert A.denominator==B.denominator==1
    rows=[]
    for row in inp['rows']:
        u=r.F(row['parameter_u'])
        if not u:continue
        s=1/u;f=s**3+A*s+B
        finite=[x['prime'] for x in row['finite_local_conditions']]
        # Prove the support contains 2 and all primes dividing f(s) or denominator(s).
        for n in (abs(f.numerator),f.denominator,s.denominator):
            for p in finite:
                while n%p==0:n//=p
            assert n==1
        assert 2 in finite
        vectors=[]
        for p in finite+['infinity']:
            v=r.pack([int(hilbert(f,r.F(P[0])-s,p)==-1) for P in pts])
            vectors.append({'place':p,'vector':v,'new':p in row['newly_bad_primes_relative_to_anchor']})
        allsum=0;newsum=0;boundary=0
        for v in vectors:
            allsum^=v['vector']
            if v['new']:newsum^=v['vector']
            else:boundary^=v['vector']
        assert allsum==0 and newsum==boundary
        local=next(x for x in prior['local_rows'] if x['u']==int(u))
        matches=[]
        for c in local['root_character_checks']:
            if c['eligible']:
                h=next(x['vector'] for x in vectors if x['place']==c['prime'])
                assert h==c['predicted_constraint']
                matches.append(c['prime'])
        newvectors=[v['vector'] for v in vectors if v['new']]
        inherited=[v for p in local['local_places'] if not p['new'] for v in p['constraint_basis']]
        inherited_basis=r.basis(inherited)
        assert r.reduce(boundary,inherited_basis)==0
        relative=[r.reduce(v,inherited_basis) for v in newvectors]
        relations=lc.orthogonal([r.pack([(v>>j)&1 for v in relative]) for j in range(20)],len(relative))
        assert relations==[(1<<len(relative))-1]
        predicted_codimension=r.rank(inherited)+len(relative)-1+local['real_constraint_rank']
        assert predicted_codimension==20-local['dimension']
        rows.append({'u':int(u),'evaluation_s':str(s),'f_at_s':str(f),'hilbert_vectors':vectors,
          'new_vector_sum':newsum,'boundary_vector_sum':boundary,'product_relation_verified':True,
          'root_character_identifications_verified_at':matches,'new_prime_count':len(newvectors),'new_vector_rank':r.rank(newvectors),
          'new_relation_dimension':len(newvectors)-r.rank(newvectors),
          'pure_new_product_relation':boundary==0,
          'nonzero_boundary_places':[v['place'] for v in vectors if not v['new'] and v['vector']],
          'new_prime_order':[v['place'] for v in vectors if v['new']],
          'relation_basis_mod_inherited_conditions':relations,
          'local_codimension_from_old_plus_new_minus_one_plus_real':predicted_codimension,
          'additional_new_relations_beyond_product':len(newvectors)-r.rank(newvectors)-int(boundary==0 and bool(newvectors))})
    out={'schema':'rank-jump.reciprocity.v1','input_sha256':r.digest(lc.INPUT.read_bytes()),'collision_report_sha256':r.digest(lc.OUTPUT.read_bytes()),
      'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'script_sha256':r.digest(Path(__file__).read_bytes()),'rows':rows,
      'omitted_place_argument':'For odd p outside support, s,A,B,f(s) are integral and f(s) is a unit. If x(P) is nonintegral then its valuation is even, so the Hilbert symbol of the unit f(s) with x(P)-s is 1. If x(P)-s is a unit both arguments are units. Otherwise f(s) is congruent to y(P)^2 modulo p and is a nonzero residue square. Thus the symbol is again 1.',
      'scope':'Rational Hilbert reciprocity among local incidence bits only. No Cassels-Tate form, soluble point, Selmer complement or rank bound is inferred.'}
    if check:
        if r.read(OUTPUT)!=out:raise ValueError('reciprocity replay mismatch')
        print('PASS reciprocity replay')
    else:r.write_new(OUTPUT,out)
    for x in rows:print('u',x['u'],'new',x['new_prime_count'],'rank',x['new_vector_rank'],'product pure',x['pure_new_product_relation'],'boundary',x['nonzero_boundary_places'],'extra relations',x['additional_new_relations_beyond_product'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('build','check'));a=p.parse_args();build(a.mode=='check')
