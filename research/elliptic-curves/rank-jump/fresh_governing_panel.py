#!/usr/bin/env python3
"""Masked, bounded inherited governing/CT baseline on fresh and historic fibres."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'FRESH_GOVERNING_PANEL_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_fresh_governing_panel_inputs_v1.json'
MANIFEST=r.OUT/'rank_jump_fresh_governing_panel_manifest_v1.json'
OUTPUT=r.OUT/'rank_jump_fresh_governing_panel_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-governing-panel-v1'
LOCAL=r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,Path(r.__file__),LOCAL)}


def export():
    sources={}
    def read(p):
        raw=p.read_bytes(); sources[str(p.relative_to(r.ROOT))]=r.digest(raw);return __import__('json').loads(raw)
    inventory=read(r.OUT/r.read(PROTOCOL)['inventory'])
    old=read(r.INPUT)['rows']
    compact=read(r.OUT/'compact192_r17_results_v1.json')['curves']
    large=read(r.OUT/'full11952_64_r17_results_v1.json')['curves']
    specs=[];pairs=[];pool={}
    def shape(row):
        model,_=r.short(row.get('curve',row.get('model')),[])
        return max(max(abs(r.F(a).numerator).bit_length(),r.F(a).denominator.bit_length()) for a in model),max(abs(r.F(row['parameter']).numerator),r.F(row['parameter']).denominator).bit_length()
    def add(row,kind,source,ident=None):
        key=ident or row['family']+':'+row['parameter']
        if key in pool:return pool[key]
        n=len(specs);token=f'case-{n:02d}';pool[key]=token
        gp=row['generic_points'];model=row.get('curve',row.get('model'))
        r.short(model,gp) # membership check only on explicitly generic coordinates
        g=row.get('generic_rank',16 if row['family'].startswith('a1-') else 17)
        assert len(gp)==g
        rank=row['rank_lower_bound']
        specs.append(({'token':token,'model':model,'generic_sections':gp},
          {'token':token,'id':key,'kind':kind,'family':row['family'],'parameter':row.get('parameter'),
           'generic_rank':g,'retained_rank_lower_bound':rank,'observed_quotient_rank':rank-g,
           'full_rank':'UNKNOWN','source':source,'completed_boxes':row.get('completed_boxes',row.get('chart_count')),
           'short_coefficient_bits':shape(row)[0] if row.get('parameter') else None,
           'parameter_height_bits':shape(row)[1] if row.get('parameter') else None}))
        return token
    for high in inventory['curves']:
        if high['rank_lower_bound']<27:continue
        source=high['source_certificate'];row=read(r.OUT/source)['curves'][high['source_curve_index']]
        assert row['family']==high['family'] and row['parameter']==high['parameter'] and row['rank_lower_bound']==high['rank_lower_bound']
        h=add(row,'fresh_high',source,high['id']);family=row['family']
        if family.startswith('a1-'):
            candidates=[x for x in old if x['family']==family and x.get('generic_points') and x['rank_lower_bound']<22]
            lowgain=min(x['rank_lower_bound']-x['generic_rank'] for x in candidates)
            candidates=[x for x in candidates if x['rank_lower_bound']-x['generic_rank']==lowgain]
            low_source=r.INPUT.name
        else:
            candidates=large if family=='11952' and shape(row)[1]>13 else compact
            low_source='full11952_64_r17_results_v1.json' if candidates is large else 'compact192_r17_results_v1.json'
            candidates=[x for x in candidates if x['family']==family and x['rank_lower_bound']==17]
        hi=shape(row);low=min(candidates,key=lambda x:(abs(shape(x)[0]-hi[0]),abs(shape(x)[1]-hi[1]),x['parameter']))
        l=add(low,'matched_low',low_source)
        pairs.append({'high':h,'low':l,'coefficient_bit_gap':abs(shape(low)[0]-hi[0]),'parameter_height_bit_gap':abs(shape(low)[1]-hi[1]),
                      'exposure_caveat':'Low gain is censored; same family, not necessarily same search exposure or cohort.'})
    # The historical public first17 interpretation is licensed only by this exact lineage certificate.
    lineage=read(r.ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json')
    assert lineage['status']=='PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS'
    parents=read(r.ROOT/'elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json')
    assert parents['status']=='PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS'
    for ident in ['ICARM-356','ICARM-385','ICARM-398']:
        row=next(x for x in old if x['id']==ident);add(row,'historic_high',r.INPUT.name,ident)
    r.write_new(INPUT,{'schema':'rank-jump.fresh-governing-panel-inputs.v1','cases':[x[0] for x in specs]})
    r.write_new(MANIFEST,{'schema':'rank-jump.fresh-governing-panel-manifest.v1','sources':sources,
       'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'masked_input_sha256':r.digest(INPUT.read_bytes()),
       'rows':[x[1] for x in specs],'pairs':pairs,
       'boundary':'Labels select retrospective cases only. Arithmetic workers read neither this manifest nor its sources.'})
    print('frozen',len(specs),'cases',len(pairs),'pairs')


def case(token):
    row=next(x for x in r.read(INPUT)['cases'] if x['token']==token)
    assert set(row)=={'token','model','generic_sections'}
    return row


def model_data(token):
    from sage.all import QQ,ZZ,PolynomialRing,lcm
    row=case(token);model,points=r.short(row['model'],row['generic_sections'])
    # One rational isomorphism, with no minimization or exceptional inputs.
    scale=lcm(QQ(x).denominator() for x in model)
    A=ZZ(QQ(model[3])*scale**4);B=ZZ(QQ(model[4])*scale**6)
    f=PolynomialRing(QQ,'z')([B,A,0,1])
    pts=[(QQ(x)*scale**2,QQ(y)*scale**3) for x,y in points]
    assert all(y*y==f(x) for x,y in pts)
    return f,pts,scale


def factor_worker(token):
    from sage.all import ZZ,pari
    f,points,scale=model_data(token)
    disc=ZZ(16*f.discriminant());fac=list(disc.abs().factor(proof=True))
    assert all(p.is_prime(proof=True) for p,e in fac)
    return {'status':'PASS','integral_cubic_ascending':list(map(str,f.list())),'scale':str(scale),
            'elliptic_discriminant':str(disc),'factors':[[int(p),int(e)] for p,e in fac]}


def local_worker(token):
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari,GF,matrix,EllipticCurve
    from sage.version import version
    import local_collision as lc
    sys.path.insert(0,str(LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    f,pts,scale=model_data(token);factor=r.read(WORK/f'{token}-factor.json')
    if factor['status']!='PASS':return {'status':'UNKNOWN','reason':'factor stage incomplete'}
    primes=[p for p,e in factor['factors']]
    product=ZZ(1)
    for p,e in factor['factors']:
        assert ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert product==abs(16*f.discriminant())
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f))
    E=EllipticCurve([0,0,0,f[1],f[0]])
    S=[2]+[p for p in primes if p!=2 and E.local_data(p,proof=True).conductor_valuation()>0]
    gammas=[];initial_ideals=[];norms=[]
    for x,y in pts:
        d=ZZ(x.denominator()).sqrt();assert d in ZZ
        a,b=ZZ(x*d*d),ZZ(y*d**3);gamma=pari(a)-pari(d*d)*theta
        assert pari.nfeltnorm(nf,gamma)==b*b
        gammas.append(gamma);initial_ideals.append(pari.idealadd(nf,b,gamma));norms.append(abs(b))
    # Fresh finite-character replay proves generic mod2 independence.
    model=['0','0','0',str(f[1]),str(f[0])];blocks=[];signatures=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['generic_independence_prime_bound']):
        roots=r.roots_at(model[3],model[4],p)
        if roots:blocks.append((p,roots))
        signatures=[r.point_signature(model,list(map(str,P)),blocks) for P in pts]
        if r.rank(signatures)==len(pts):break
    assert r.rank(signatures)==len(pts)
    locals=[];columns=[]
    for p in S:
        chars=LocalSquareclasses(nf,p);sigs=[list(chars.signature(g)) for g in gammas]
        for j in range(len(sigs[0])):columns.append(r.pack(s[j] for s in sigs))
        locals.append({'place':p,'signatures':sigs,'point_kummer_dimension':chars.point_kummer_dimension,
                       'generic_image_dimension':r.rank(list(map(r.pack,sigs)))})
    roots=f.roots(AA,multiplicities=False)
    signs=[[int(x<a) for a in roots] for x,y in pts]
    for j in range(len(roots)):columns.append(r.pack(s[j] for s in signs))
    locals.append({'place':'infinity','signatures':signs,'point_kummer_dimension':int(len(roots)==3),
                   'generic_image_dimension':r.rank(list(map(r.pack,signs)))})
    masks=lc.orthogonal(columns,len(pts));k=len(masks)
    assert matrix(GF(2),[[v>>i&1 for i in range(len(pts))] for v in columns]).right_kernel().dimension()==k
    prime_ideals=[P for p in primes for P in pari.idealprimedec(nf,p)]
    betas=[];ideals=[];records=[]
    enc=lambda a:[str(pari.lift(a).polcoef(i)) for i in range(3)]
    mat=lambda H:[[str(H[i,j]) for j in range(3)] for i in range(3)]
    for mask in masks:
        beta=pari.Mod(1,pari(f));I=pari.idealhnf(nf,1)
        for i,g in enumerate(gammas):
            if mask>>i&1:beta*=g;I=pari.idealmul(nf,I,initial_ideals[i])
        for P in prime_ideals:
            v=int(pari.idealval(nf,beta,P));assert v%2==0
            correction=v//2-int(pari.idealval(nf,I,P))
            if correction:I=pari.idealmul(nf,I,pari.idealpow(nf,P,correction))
        I=pari.idealhnf(nf,I);assert pari.idealpow(nf,I,2)==pari.idealhnf(nf,beta)
        betas.append(beta);ideals.append(I);records.append({'generic_mask':mask,'beta_ascending':enc(beta),'half_ideal_hnf':mat(I)})
    A=[[None]*k for _ in range(k)];cols=[]
    for j,I in enumerate(ideals):
        reduced,alpha=pari.idealred(nf,[I,1]);assert pari.idealmul(nf,reduced,alpha)==I
        good=reduced
        for p in S:
            for P in pari.idealprimedec(nf,p):
                e=int(pari.idealval(nf,good,P));assert e>=0
                if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
        good=pari.idealhnf(nf,good);N=ZZ(pari.idealnorm(nf,good))
        assert all(N%p for p in S) and N>0
        cyclic=good[1,1]==good[2,2]==1 and good[0,0]==N
        evaluations=[]
        if cyclic:
            for i,beta in enumerate(betas):
                coords=pari.nfalgtobasis(nf,beta)
                residue=ZZ(coords[0]-good[0,1]*coords[1]-good[0,2]*coords[2])%N
                gcd=residue.gcd(N)
                if gcd==1:
                    A[i][j]=int(pari.kronecker(residue,N)==-1)
                    evaluations.append({'residue':str(residue),'artin_bit':A[i][j]})
                else:evaluations.append({'status':'UNKNOWN','nonunit_gcd':str(gcd)})
        cols.append({'coprime_ideal_hnf':mat(good),'norm':str(N),'cyclic':bool(cyclic),'evaluations':evaluations,
                     'reduced_ideal_hnf':mat(reduced),'multiplier_ascending':enc(pari.nfbasistoalg(nf,alpha))})
    complete=all(x is not None for row in A for x in row)
    M=[[A[i][j]^A[j][i] for j in range(k)] for i in range(k)] if complete else None
    galois=r.galois(model)
    return {'status':'PASS' if complete else 'PARTIAL','software':{'sage':version,'pari':str(pari.version())},
      'S_finite':S,'generic_dimension':len(pts),'strict_generic_dimension':k,'generic_strict_masks':masks,
      'independence_blocks':blocks,'independence_signatures':signatures,'local':locals,'galois':galois,
      'field_discriminant':str(nf.disc()),'maximal_order_basis':list(map(str,nf.nf_get_zk())),
      'class_records':records,'artin_columns':cols,'Artin_matrix':A,'minus_twist_CT_matrix':M,
      'Artin_rank':r.rank(list(map(r.pack,A))) if complete else None,
      'minus_twist_CT_rank':r.rank(list(map(r.pack,M))) if complete else None,
      'original_restricted_CT_rank':0,
      'additional_class_basis':'UNKNOWN','additional_quotient_CT':'UNKNOWN','full_CT':'UNKNOWN',
      'boundary':'Inherited generic strict block only. Original CT zero follows from generic rational sections. Switch rank is exact on this block; no exceptional quotient or whole-curve rank claim.'}


def independent_worker(token):
    from sage.all import pari
    f,_,_=model_data(token);factor=r.read(WORK/f'{token}-factor.json')
    if factor['status']!='PASS':return {'status':'UNKNOWN','reason':'factor stage incomplete'}
    # Equation-only path: no generic section, Kummer mask or exceptional source used by bnf.
    pari.setrand(20260906)
    nf=pari.nfinit([pari(f),[p for p,e in factor['factors']]])
    bnf=pari.bnfinit(nf,1)
    certified=bool(pari.bnfcertify(bnf)==1)
    return {'status':'PASS' if certified else 'UNKNOWN','certified':certified,
            'class_group_cyclic':list(map(str,bnf.bnf_get_cyc())),
            'class_group_generators':list(map(str,bnf.bnf_get_gen())),
            'units':list(map(str,bnf.bnf_get_fu())),
            'boundary':'Ordinary class group only. Localized strict basis, quotient and CT require further certified steps.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[];limits=r.read(PROTOCOL)['limits']
    for c in r.read(INPUT)['cases']:
        token=c['token'];row={'token':token}
        for stage,cap in [('factor',limits['factor_seconds_per_case']),('local',limits['local_and_artin_seconds_per_case']),('independent',limits['independent_class_seconds_per_case'])]:
            path=WORK/f'{token}-{stage}.json'
            initial=WORK.with_name(WORK.name+'-initial')/path.name
            if not path.exists() and stage in ('factor','independent') and initial.exists():
                prior=r.read(initial)
                assert prior['bindings'][str(INPUT.relative_to(r.ROOT))]==r.digest(INPUT.read_bytes())
                # Preserve completed equation-only work and genuine bounded failures.
                # The repaired local worker rechecks factors before using them.
                prior.pop('bindings')
                r.write_new(path,{'bindings':bindings(),'resumed_checkpoint_sha256':r.digest(initial.read_bytes()),**prior})
            if not path.exists():
                with (WORK/f'{token}-{stage}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token,'--stage',stage],stdout=log,stderr=log,timeout=cap)
                        reason=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:reason=f'{cap}-second timeout'
                if reason:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':reason})
            data=r.read(path);assert data['bindings']==bindings();row[stage]=data
            print(token,stage,data['status'],data.get('strict_generic_dimension'),data.get('minus_twist_CT_rank'),flush=True)
        rows.append(row)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-governing-panel.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--token');p.add_argument('--stage');args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        result={'factor':factor_worker,'local':local_worker,'independent':independent_worker}[args.stage](args.token)
        r.write_new(WORK/f'{args.token}-{args.stage}.json',{'bindings':bindings(),**result})
