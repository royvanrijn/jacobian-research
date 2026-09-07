#!/usr/bin/env python3
"""Verify norm-projection candidates and factor-free isolated-norm exclusions."""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import gcd,isqrt,lcm,prod
from pathlib import Path
import sys
import retrospective as r
import fresh_norm_projection as source
from verify_unpointed_governing_norm import Algebra

PROTOCOL=Path(__file__).with_name('FRESH_NORM_ISOLATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_norm_projection_verification_v1.json'


def isolation(row, survivor, discriminant):
    m=row['generic_dimension'];selected=[j for j in range(12) if survivor['coefficient_mask']>>(m+j)&1]
    assert selected
    bad=abs(discriminant)*2;contents=[]
    for j in selected:
        coeff=list(map(Q,row['candidates'][j]['alpha_ascending']));den=lcm(*(c.denominator for c in coeff))
        ints=[int(c*den) for c in coeff];content=gcd(*ints);assert content
        bad*=den*content;contents.append({'candidate_index':j,'denominator':den,'integer_content':content})
    for i in selected:
        norm=abs(int(row['candidates'][i]['norm_alpha']));avoid=bad*prod(abs(int(row['candidates'][j]['norm_alpha'])) for j in selected if j!=i)
        remaining=norm;steps=[]
        while (g:=gcd(remaining,avoid))>1:remaining//=g;steps.append(str(g))
        assert gcd(remaining,avoid)==1 and prod(map(int,steps))*remaining==norm
        root=isqrt(remaining)
        if root*root!=remaining:
            assert root*root<remaining<(root+1)**2
            return {'status':'PROVED_RAMIFIED_OUTSIDE_S','coefficient_mask':survivor['coefficient_mask'],
                'selected_candidate_indices':selected,'isolated_candidate_index':i,'coefficient_denominators_and_contents':contents,
                'norm_alpha_absolute':str(norm),'excluded_support_integer':str(avoid),'gcd_removal_steps':steps,
                'isolated_coprime_remainder':str(remaining),'floor_square_root':str(root),
                'witness':'Some odd prime dividing the isolated remainder to odd order is a good unramified cubic prime. Only this alpha has nonunit norm there; its polynomial has degree<3 and is nonzero modulo that prime, so one cubic component is a unit. Norm projection gives odd valuation on that component. Generic point classes cannot cancel it.'}
    return {'status':'UNKNOWN','coefficient_mask':survivor['coefficient_mask'],'reason':'No isolated nonsquare norm remainder'}


def compute():
    from sage.all import QQ,ZZ,AA,pari,GF,matrix
    sys.path.insert(0,str(source.base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    pari.allocatemem(64000000,r.read(source.PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    data=r.read(source.OUTPUT)
    for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    rows=[];norm_checks=0;local_checks=0;character_checks=0
    for row in data['rows']:
        assert row['status']=='PASS';token=row['token'];old,f,pts,primes,nf=source.setup(token)
        m=len(pts);k=old['local']['strict_generic_dimension'];S=old['local']['S_finite'];K=Algebra(list(map(str,f.list())))
        enc=lambda x:[str(pari.lift(x).polcoef(i)) for i in range(3)]
        th=pari.Mod('z',pari(f));basis=[K.elt(enc(pari.Mod(b,pari(f)))) for b in nf.nf_get_zk()]
        assert [list(map(str,b)) for b in basis]==row['maximal_order_basis_ascending']
        coords=[v for v in product((-1,0,1),repeat=3) if any(v) and next(x for x in v if x)>0 and v!=(1,0,0)]
        betas=[K.elt([str(x),'-1']) for x,y in pts]
        for v,c in zip(coords,row['candidates']):
            assert list(v)==c['basis_coordinates']
            alpha=tuple(sum(v[j]*basis[j][i] for j in range(3)) for i in range(3))
            assert alpha==K.elt(c['alpha_ascending']);N=K.norm(alpha)
            assert N.denominator==1 and str(N)==c['norm_alpha'] and N!=0
            beta=tuple(N*x for x in alpha);assert beta==K.elt(c['beta_ascending']) and K.norm(beta)==N**4
            betas.append(beta);norm_checks+=1
        assert len(betas)==m+12
        elements=[pari(f.parent()([QQ(str(c)) for c in b]))(th) for b in betas]
        constraints=[]
        for loc in row['local']:
            if loc['place']=='infinity':
                roots=f.roots(AA,multiplicities=False)
                sigs=[[int(f.parent()([QQ(str(c)) for c in b])(a)<0) for a in roots] for b in betas]
            else:
                chars=LocalSquareclasses(nf,loc['place']);sigs=[list(chars.signature(b)) for b in elements];local_checks+=len(betas)
            assert sigs==loc['signatures'];constraints.extend(list(zip(*sigs)))
        assert matrix(GF(2),constraints).right_kernel().dimension()==row['local_only_kernel_dimension']
        for probe in row['outside_S_probes']:
            p=probe['prime'];assert p not in S and p<=r.read(source.PROTOCOL)['limits']['probe_prime_bound']
            P=pari.idealprimedec(nf,p)[probe['prime_index']];vals=[int(pari.idealval(nf,b,P)) for b in elements]
            assert vals==probe['valuations'] and all(v%2==0 for v in vals[:m]);constraints.append([v%2 for v in vals])
        kernel=matrix(GF(2),constraints).right_kernel();generic_kernel=matrix(GF(2),[c[:m] for c in constraints]).right_kernel()
        assert generic_kernel.dimension()==k and kernel.dimension()==row['final_kernel_dimension']
        masks=row['final_kernel_masks'];assert len(masks)==kernel.dimension()
        for mask in masks:assert matrix(GF(2),[[mask>>i&1 for i in range(m+12)]]).row(0) in kernel
        assert r.rank(masks)==len(masks)
        survivors=row['survivors'];assert len(survivors)==kernel.dimension()-k==row['candidate_coefficient_quotient_dimension']
        assert r.rank([s['coefficient_mask']>>m for s in survivors])==len(survivors)
        classes=betas[:m]
        for s in survivors:
            mask=s['coefficient_mask'];assert matrix(GF(2),[[mask>>i&1 for i in range(m+12)]]).row(0) in kernel
            beta=K.elt([1])
            for i,b in enumerate(betas):
                if mask>>i&1:beta=K.mul(beta,b)
            assert beta==K.elt(s['beta_ascending']);classes.append(beta)
        signatures=[0]*len(classes)
        for j,(p,x) in enumerate(row['character_places']):
            assert p!=2 and int(f.discriminant())%p and f(x)%p==0
            for i,b in enumerate(classes):
                assert all(c.denominator%p for c in b)
                v=sum(c.numerator*pow(c.denominator,-1,p)*pow(x,n,p) for n,c in enumerate(b))%p;assert v
                signatures[i]|=int(pow(v,(p-1)//2,p)==p-1)<<j;character_checks+=1
        assert signatures==row['character_signatures'] and r.rank(signatures[:m])==m
        increment=r.rank(signatures)-m;assert increment==row['certified_global_squareclass_increment_beyond_G']
        assert len(survivors)<=1 # basis-wise ramification alone would not suffice above dimension one
        exclusions=[isolation(row,s,int(f.discriminant())) for s in survivors]
        assert all(x['status']=='PROVED_RAMIFIED_OUTSIDE_S' for x in exclusions)
        if survivors:assert increment==1 and not survivors[0]['is_global_square']
        rows.append({'token':token,'generic_dimension':m,'generic_strict_dimension':k,'norm_projection_candidates':12,
            'coefficient_quotient_dimension_after_finite_gates':len(survivors),'global_squareclass_increment':increment,
            'outside_S_exclusions':exclusions,'additional_strict_dimension_of_generated_span':0,
            'total_strict_dimension_of_generated_span':k,'status':'PASS'})
        print(token,'PASS','candidate quotient',len(survivors),'additional strict image 0',flush=True)
    files=(Path(__file__),PROTOCOL,source.OUTPUT,Path(source.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return {'schema':'rank-jump.fresh-norm-projection-verification.v1','status':'PASS',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},'rows':rows,
        'norm_projection_identities':norm_checks,'replayed_finite_local_signatures':local_checks,
        'independent_finite_field_character_evaluations':character_checks,
        'boundary':'Zero additional strict image for this fixed twelve-class dictionary on all eleven fibres. Two independent global classes fail outside-S ramification by factor-free certificates. No full Selmer dimension or rank bound.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['norm_projection_identities'],'norm identities',flush=True)
