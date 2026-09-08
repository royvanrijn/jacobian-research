#!/usr/bin/env sage-python
"""Independent finite-field-norm and branch-subset replay of local images.

No constructor imports, point searches, number fields, or class groups.
The inherited-divisor theorem is a hash-bound dependency, not re-proved here.
"""
import hashlib,json,signal
from collections import Counter
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,prod
from sage.version import version as sage_version

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_good_local_images_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def bindings(d):
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
def binary_rank(rows):
    pivots={}
    for row in rows:
        n=sum(int(b)<<i for i,b in enumerate(row))
        while n:
            j=n.bit_length()-1
            if j in pivots:n^=pivots[j]
            else:pivots[j]=n;break
    return len(pivots)
def subset_gate(degrees):
    permutation=[];offset=0
    for d in degrees:
        permutation.extend([offset+(j+1)%d for j in range(d)]);offset+=d
    assert offset==6
    def action(n):return sum(((n>>i)&1)<<permutation[i] for i in range(6))
    def canon(n):return min(n,n^63)
    torsion={canon(n) for n in range(64) if n.bit_count()%2==0}
    fixed=sorted(n for n in torsion if canon(action(n))==n)
    assert len(torsion)==16
    dim=len(fixed).bit_length()-1;assert 2**dim==len(fixed)
    assert dim==len(degrees)-1-int(any(d%2 for d in degrees))
    # Six odd theta characteristics are branch points; ten even ones are
    # unordered 3+3 partitions. A fixed theta lifts at good odd reduction.
    theta_odd=[i for i in range(6) if permutation[i]==i]
    triples={canon(n) for n in range(64) if n.bit_count()==3}
    theta_even=sorted(n for n in triples if canon(action(n))==n)
    assert len(triples)==10 and (theta_odd or theta_even)
    return {'degrees':degrees,'permutation':permutation,'fixed_torsion_masks':fixed,
            'dimension':dim,'fixed_odd_theta_indices':theta_odd,
            'fixed_even_theta_masks':theta_even,'local_fake_kernel_zero':True}
def verify():
    protocol=read(OUT/'protocol.json');manifest=read(OUT/'manifest.json')
    bindings(protocol);bindings(manifest)
    generic_path=ART/'det1092_rr_full_inherited_jacobian_v1.json'
    generic=read(generic_path)['generic_r_functions']
    old_protocol=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
    assert protocol['primes']==old_protocol['primes'] and len(protocol['primes'])==64
    dependency_paths=[ART/'det1092_rr_full_inherited_jacobian_replay_v1.json',
                      ART/'det1092_rr_generic_point_controls_v2/panel-replay.json']
    for path in dependency_paths:bindings(read(path))
    assert read(dependency_paths[0])['inherited_NS_restriction_rank']==17
    assert read(dependency_paths[1])['status']=='PASS_NINE_GENERIC_ELLIPTIC_POINTS_GIVE_NON_GENERIC_RATIONAL_JACOBIAN_CLASSES'
    R=PolynomialRing(QQ,'T');totals=Counter();results=[];patterns={};newly_usable=0
    all_paths=[OUT/'protocol.json',OUT/'manifest.json',generic_path]+dependency_paths
    for i in range(10):
        path=OUT/('case-%02d.json'%i);d=read(path);all_paths.append(path)
        src_path=ROOT/d['source_path'];source=read(src_path)
        assert protocol['inputs'][d['source_path']]==sha(src_path)
        old_trials=read(generic_path)['trials'] if i==9 else source['trials']
        old_by_p={row['p']:row for row in old_trials}
        curve=source['curve'] if i==9 else source
        q=R(curve['q']);scale=QQ(curve['scale']);u=QQ(curve['u'])
        x0=QQ(source['inherited_x_values'][source['base_pair_index']] if i==9 else source['base_x'])
        assert d['q']==list(map(str,q.list())) and QQ(d['scale'])==scale
        assert QQ(d['base_x'])==x0 and QQ(d['u'])==u
        assert (scale*q(x0)).is_square() and scale*q(x0)!=0
        gs=[]
        for row,saved in zip(generic,d['generic_divisors']):
            raw=R(row['numerator'])-u*R(row['denominator']);g=R(saved)
            assert g.is_monic() and g.degree()==row['degree']
            assert raw==raw.leading_coefficient()*g
            gs.append(g)
        assert len(gs)==len(d['generic_divisors'])==17
        if i<9:
            dependency=ART/'det1092_rr_generic_point_controls_v2'/('case-%02d-replay.json'%i)
            bindings(read(dependency));all_paths.append(dependency)
            assert read(dependency)['status']=='PASS_GENERIC_ELLIPTIC_POINT_GIVES_NON_GENERIC_RATIONAL_JACOBIAN_CLASS'
        counts=Counter();dimensions=Counter();new_case=0
        assert [row['p'] for row in d['trials']]==protocol['primes']
        for trial in d['trials']:
            p=trial['p'];assert ZZ(p).is_prime() and p>2
            F=GF(p);S=PolynomialRing(F,'X');qm=QQ(trial['q_multiplier'])
            # Check integrality and primitiveness directly, not a replay of
            # the constructor's minimum-valuation choice.
            assert qm>0 and qm.numerator().prime_to_m_part(p)==1 and qm.denominator().prime_to_m_part(p)==1
            normalized=q*qm
            assert all(a.valuation(p)>=0 for a in normalized if a)
            assert any(a.valuation(p)==0 for a in normalized if a)
            f=S(normalized.list());twist_v=(scale/qm).valuation(p)
            assert trial['adjusted_scale_valuation']==twist_v
            counts[trial['status']]+=1
            if f.degree()!=6 or f.gcd(f.derivative())!=1:
                assert trial['status']=='UNRESOLVED_NOT_SQUAREFREE_DEGREE_SIX_REDUCTION';continue
            if twist_v%2:
                assert trial['status']=='UNRESOLVED_ODD_SCALAR_TWIST';continue
            factors=[S(row) for row in trial['factors']]
            assert prod(factors)==f.monic()
            assert all(h.is_monic() and h.is_irreducible() for h in factors)
            degrees=[int(h.degree()) for h in factors];assert trial['degrees']==degrees
            gate=subset_gate(degrees);patterns[tuple(degrees)]=gate
            dim=gate['dimension'];assert trial['local_Kummer_dimension']==dim
            anchor=R(trial['anchor']);assert anchor.degree()==1 and anchor(x0)==0
            assert all(a.valuation(p)>=0 for a in anchor if a)
            anchor_p=S(anchor.list())
            if anchor_p.gcd(f)!=1:
                assert trial['status']=='UNRESOLVED_BASEPOINT_MEETS_BRANCH';continue
            kept=[];omitted=[];reduced=[]
            assert len(trial['divisor_p_normalizations'])==17
            for j,(g,mult) in enumerate(zip(gs,trial['divisor_p_normalizations'])):
                mult=QQ(mult)
                assert mult>0 and mult.numerator().prime_to_m_part(p)==1 and mult.denominator().prime_to_m_part(p)==1
                normalized=g*mult
                assert all(a.valuation(p)>=0 for a in normalized if a)
                assert any(a.valuation(p)==0 for a in normalized if a)
                gp=S(normalized.list())
                if gp.gcd(f)!=1:omitted.append(j)
                else:kept.append(j);reduced.append(gp)
            assert kept==trial['generic_divisor_indices'] and omitted==trial['omitted_nonunit_divisors']
            raw=[]
            for h in factors:
                degree=h.degree()
                if degree==1:L=F;theta=-h[0]
                else:L=GF(ZZ(p)**degree,name='theta',modulus=h);theta=L.gen()
                base=anchor_p(theta);row=[]
                for j,gp in zip(kept,reduced):
                    # The rational divisor degree, not its possibly dropped
                    # reduction degree, determines its sign and basepoint power.
                    value=(-1)**gs[j].degree()*gp(theta)/base**gs[j].degree();assert value
                    norm=prod(value**(ZZ(p)**e) for e in range(degree))
                    legendre=F(norm)**((p-1)//2);assert legendre in [F(1),F(-1)]
                    row.append(int(legendre==F(-1)))
                raw.append(row)
            assert raw==trial['raw_character_rows']
            assert all(sum(row[j] for row in raw)%2==0 for j in range(len(kept)))
            parity=vector(GF(2),[d%2 for d in degrees]);k=len(factors)
            proj=matrix(GF(2),trial['scalar_quotient_projection'])
            assert proj.ncols()==k and not proj*parity
            assert binary_rank(proj.rows())==k-int(bool(parity))
            B=matrix(GF(2),k,len(kept),sum(raw,[]));image=proj*B
            rows=[list(map(int,row)) for row in image.rows()]
            assert rows==trial['image_rows']
            rank=binary_rank(rows);assert rank==trial['generic_image_rank']<=dim
            positions=[kept.index(j) for j in trial['basis_generic_divisors']]
            assert len(set(positions))==rank and binary_rank(image.matrix_from_columns(positions).rows())==rank
            if rank==dim:
                assert trial['status']=='COMPLETE_LOCAL_KUMMER_IMAGE';dimensions[dim]+=1
                if old_by_p[p]['status']!='PASS_LOCAL_KUMMER_BLOCK':new_case+=1
            else:assert trial['status']=='INCOMPLETE_GENERIC_LOCAL_SPAN'
        assert dict(counts)==d['counts']==manifest['cases'][i]['counts']
        totals.update(counts);newly_usable+=new_case
        result={'case_index':i,'counts':dict(counts),'complete_image_dimensions':dict(dimensions),
                'newly_usable_previously_skipped_pairs':new_case,
                'classification':'verified application; independent replay',
                'status':'PASS_INDEPENDENT_GOOD_LOCAL_IMAGE_CASE',
                'inputs':{str(path.relative_to(ROOT)):sha(path)},'checker_sha256':sha(Path(__file__))}
        retain(OUT/('case-%02d-replay.json'%i),result);results.append(result)
    assert dict(totals)==manifest['totals'] and sum(totals.values())==640
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_FROZEN_GOOD_LOCAL_IMAGES',
            'sage_version':sage_version,'cases':results,'totals':dict(totals),
            'newly_usable_previously_skipped_pairs':newly_usable,
            'branch_subset_proofs':[patterns[key] for key in sorted(patterns)],
            'conclusion':'At every complete pair, inherited generic section divisors span the entire true local2-Kummer image; no marked exceptional or generic point is needed.',
            'scope':'Only recorded good odd places. Full global Selmer, bad places,2,real place,and untested primes remain unresolved. Individual local spanning is not simultaneous global spanning.',
            'limits':{**protocol['limits'],'independent_replay_wall_seconds':25},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in all_paths},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();retain(OUT/'replay.json',result)
    print(result['status'],result['totals'],flush=True)
