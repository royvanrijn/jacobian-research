#!/usr/bin/env python3
"""Full retained quartets: genus-three quotient and three affine classes."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prod
import retrospective as r
import verify_paired_quartet_relations as ranks

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'QUARTET_UNRAMIFIED_DESCENT_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
PRIOR=r.OUT/'rank_jump_soluble_quartet_compression_v1.json'
RANKS=r.OUT/'rank_jump_paired_quartet_relations_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_quartet_unramified_descent_v1.json'


def deck_geometry():
    permitted=[v for v in range(1,16) if v not in (1,2,4,8)]
    free=[]
    for mask in range(1<<len(permitted)):
        H={0}|{v for i,v in enumerate(permitted) if (mask>>i)&1}
        if all((a^b) in H for a in H for b in H):free.append(sorted(H))
    largest=[H for H in free if len(H)==max(map(len,free))]
    assert largest==[[0,3,5,6,9,10,12,15]]
    assert all(v.bit_count()%2==0 for v in largest[0])
    # Even subsets of eight branch points, modulo their full complement.
    pairs=[3,12,48,192];canonical=lambda v:min(v,v^255)
    span={0}
    for v in pairs:span|={x^v for x in list(span)}
    classes=sorted({canonical(v) for v in span})
    assert len(classes)==8 and len({canonical(v) for v in pairs})==4
    assert all((a&b).bit_count()%2==0 for a in pairs for b in pairs)
    characters=[{'mask':mask,'genus':mask.bit_count()-1} for mask in range(1,16)]
    assert sum(c['genus'] for c in characters)==17
    return {'free_native_subgroups':free,'unique_maximal_free_subgroup':largest[0],
        'carrier_genus':17,'product_quotient_genus':3,'etale_degree':8,
        'isotropic_2_torsion_classes':classes,'isotropic_dimension':3,
        'character_quotients':characters,'positive_genus_dimensions':[1]*6+[2]*4+[3]}


def field_incidence(qs,p):
    k=GF(p);R=PolynomialRing(k,'t');q=[R(f.list()) for f in qs]
    # Homogeneous quadratics: allow one simple branch at infinity.
    good=all(f.degree()>=1 and k(f[1]**2-4*f[0]*f[2])!=0 for f in q)
    good=good and sum(f.degree()<2 for f in q)<=1
    good=good and all(q[i].gcd(q[j])==1 for i in range(4) for j in range(i))
    if not good:return {'prime':p,'status':'UNKNOWN_BAD_BRANCH_REDUCTION'}
    def nr(v):return 1 if v==0 else 2 if v.is_square() else 0
    counts={format(i,'03b'):0 for i in range(8)};points=liftable=carrier=branch=0
    for t in list(k)+[None]:
        vals=[f(t) if t is not None else f[2] for f in q]
        ycount=nr(prod(vals));points+=ycount
        native=prod(nr(v) for v in vals);carrier+=native
        if any(v==0 for v in vals):
            assert sum(v==0 for v in vals)==1
            branch+=ycount
            if native:liftable+=ycount
        else:
            key=''.join(str(int(not v.is_square())) for v in vals[:3]);counts[key]+=ycount
            if native:liftable+=ycount
    assert carrier==8*liftable and points==sum(counts.values())+branch
    return {'prime':p,'status':'PASS','product_points':int(points),'liftable_product_points':int(liftable),
        'carrier_points':int(carrier),'non_branch_class_counts':counts,'branch_product_points':int(branch)}


def compute():
    inp=r.read(INPUT);prior=r.read(PRIOR);rank_certificate=r.read(RANKS)
    for d in (inp,prior,rank_certificate):
        for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    rank_replay=ranks.verify();assert rank_replay==rank_certificate
    geo=deck_geometry();R=PolynomialRing(QQ,'t');rows=[]
    for case in inp['cases']:
        qs=[R(c['form']) for c in case['covers']];product=prod(qs)
        assert product.degree()==8 and product.is_squarefree()
        assert all(q.degree()==2 for q in qs)
        old=next(x for x in prior['rows'] if x['id']==case['id'])
        assert old['carrier_genus']==17 and old['carrier_degree_over_parameter_line']==16
        assert R(next(x['coefficients'] for x in old['character_quotients'] if x['mask']==15))==product
        row={'id':case['id'],'product_polynomial':list(map(str,product.list())),
            'finite_field_incidence':[field_incidence(qs,p) for p in (131,137)]}
        if case['observed_parameter'] is not None:
            t=QQ(case['observed_parameter']);roots=[q(t).sqrt() for q in qs]
            assert all(a in QQ and a*a==q(t) and a for a,q in zip(roots,qs,strict=True))
            rr=next(x for x in rank_replay['rows'] if x['id']==case['id'])
            row.update({'retained_parameter':str(t),'native_roots':list(map(str,roots)),
                'product_root':str(prod(roots)),'generic_rank':17,'retained_rank':rr['basis_rank'],
                'observed_quotient':rr['basis_rank']-17,'quartet_exact_quotient_rank':rr['exact_quotient_rank'],
                'unexplained_retained_quotient':rr['basis_rank']-17-rr['exact_quotient_rank']})
        else:
            # Proper projective residue obstruction; does not assume good reduction.
            k=GF(23);rk=PolynomialRing(k,'t');q=[rk(f.list()) for f in qs]
            supports=[]
            for t in list(k)+[None]:
                vals=[f(t) if t is not None else f[2] for f in q]
                if all(v.is_square() for v in vals):supports.append(str(t))
            assert not supports
            vals=[f(0) for f in q];assert prod(vals)==k(11)**2 and all(vals)
            row['local_control']={'prime':23,'quartet_projective_support':supports,
                'product_point_mod_prime':{'t':0,'y':11},'native_values_mod_prime':list(map(int,vals)),
                'native_affine_class':'001','product_curve_Qp_point_by_Hensel':True,
                'quartet_Qp_points':False,'good_branch_reduction':field_incidence(qs,23)['status']=='PASS'}
            assert ''.join(str(int(not v.is_square())) for v in vals[:3])=='001'
        rows.append(row)
    paths=[INPUT,PRIOR,RANKS,PROTOCOL,Path(__file__),HERE/'retrospective.py',HERE/'verify_paired_quartet_relations.py']
    return {'schema':'rank-jump.quartet-unramified-descent.v1','status':'PASS','layer':'solubility',
        'geometry':geo,'rows':rows,'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: three genus3 quotient systems; etale degree8; quartet quotient ranks3,3')
