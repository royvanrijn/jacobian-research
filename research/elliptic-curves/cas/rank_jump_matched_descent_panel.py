#!/usr/bin/env python3
"""Frozen ten-fibre retrospective comparison; complete support or UNKNOWN.

Five pairs share their literal R17 fibration. They are not height-matched or
randomized, and302 is a separate determinant1092 reference. Never label a
certified lower bound17 an exact rank17. No point search or new factorization.
"""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'r17_60_panel_results_v1.json';ROSTER=ART/'r17_60_panel_roster_v1.json'
PROTOCOL=ART/'rank_jump_matched_descent_panel_protocol_v1.json'
OUT=ART/'rank_jump_matched_descent_panel_v1.json'
WORK=ROOT/'artifacts/local/elliptic-curves/rank-jump-matched-descent-panel-v1'
PAIRS=[['103b2-low-02','103b2-high-01'],['103b2-low-05','103b2-high-05'],
       ['074d9-low-04','074d9-low-02'],['11952-low-04','11952-high-02'],
       ['08f72-low-01','08f72-high-02']]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    if p.exists():assert read(p)==d
    else:
        p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x')as stream:json.dump(d,stream,indent=2,sort_keys=True);stream.write('\n')
def freeze():
    results={r['id']:r for r in read(INPUT)['results']};roster={r['id']:r for r in read(ROSTER)['rows']}
    rows=[];sources=[Path(__file__),INPUT,ROSTER]
    for pair_index,pair in enumerate(PAIRS):
        assert results[pair[0]]['family']==results[pair[1]]['family']
        for case in pair:
            r=results[case];packet=r['packet'];hints=[];primes=set()
            for path in sorted(ART.glob('inventory291_conductors*/r17-panel-'+case+'.json')):
                d=read(path);c=d.get('certificate',{})
                if c.get('curve')!=packet['curve']:continue
                sources.append(path);hints.append(str(path.relative_to(ROOT)))
                primes.update(int(x['prime'])for x in c.get('local_data',[]))
            rows.append({'id':case,'pair':pair_index,'family':r['family'],'parameter':r['parameter'],
                         'rank_lower_bound':r['rank_lower_bound'],'generic_rank':17,
                         'curve':packet['curve'],'points':packet['points'],'proof':packet['proof'],
                         'retained_prime_hints':list(map(str,sorted(primes))),'hint_sources':hints,
                         'roster_row':roster[case]})
    d={'schema':'rank-jump.matched-descent-panel.protocol.v1','pairs':PAIRS,'rows':rows,
       'bindings':{str(p.relative_to(ROOT)):sha(p)for p in sources},
       'selection':'Five named same-fibration pairs from the sealed60-fibre panel, selected by known lower-bound strata before descent. Three17, three21, three26, one27.',
       'matching_boundary':'Within-pair fibration matches only; parameter heights and search exposure differ.302 is a separate parent/reference. This is retrospective arithmetic, not prospective predictor validation.',
       'limits':{'workers':1,'seconds_per_case':60,'good_character_prime_bound':2000,
                 'new_factorization':0,'point_searches':0,'new_parameters':0,'bnf_calls':0},
       'failure_semantics':'Unfactored support, deficient good-prime Kummer rank, local failures and timed-out stages remain UNKNOWN. Missing data do not count as low class rank.'}
    save(PROTOCOL,d);print('FROZEN ten cases',[(r['id'],r['rank_lower_bound'],len(r['retained_prime_hints']))for r in rows],flush=True)
def worker(case):
    from sage.all import AA,QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,pari,prime_range
    sys.path.insert(0,str(Path(__file__).parent))
    from r17_60_arithmetic import native_check
    from research_runtime.local_kummer import LocalSquareclasses
    p=read(PROTOCOL);row=next(r for r in p['rows']if r['id']==case)
    dest=WORK/case;dest.mkdir(parents=True,exist_ok=True)
    assert sha(Path(__file__))==p['bindings'][str(Path(__file__).relative_to(ROOT))]
    packet={k:row[k]for k in ['curve','points','proof','rank_lower_bound']}
    native_check(packet,row['roster_row'])
    E=EllipticCurve(QQ,list(map(QQ,row['curve'])));assert E.a1()==E.a2()==E.a3()==0
    R=PolynomialRing(QQ,'z');f=R([E.a6(),E.a4(),0,1]);assert all(c in ZZ for c in f)
    disc=ZZ(f.discriminant());remainder=abs(disc);S=sorted({2,*map(int,row['retained_prime_hints'])})
    factors=[]
    for q in S:
        assert ZZ(q).is_prime(proof=True)
        v=remainder.valuation(q)
        if v:remainder//=ZZ(q)**v;factors.append([str(q),int(v)])
    preliminary={'id':case,'family':row['family'],'rank_lower_bound':row['rank_lower_bound'],
                 'protocol_sha256':sha(PROTOCOL),'cubic_ascending':list(map(str,f.list())),
                 'disc_factors':factors,'remaining_cofactor':str(remainder),'stage':'support',
                 'status':'UNKNOWN_INCOMPLETE_SUPPORT'if remainder!=1 else 'PASS_COMPLETE_SUPPORT'}
    save(dest/'support.json',preliminary)
    if remainder!=1:save(dest/'result.json',preliminary);print(case,preliminary['status'],flush=True);return
    assert f.is_irreducible()
    nf=pari.nfinit([pari(f),S]);n=row['rank_lower_bound'];gs=[];polys=[]
    for P in row['points']:
        Q=E(list(map(QQ,P)));x,y=Q.xy();d=ZZ(x.denominator()).sqrt();assert d in ZZ and d*d==x.denominator()
        aa,bb=ZZ(x*d*d),ZZ(y*d**3);poly=R([aa,-d*d]);g=pari.Mod(pari(poly),pari(f))
        assert pari.nfeltnorm(nf,g)==bb*bb and aa.gcd(d)==bb.gcd(d)==1
        polys.append(poly);gs.append(g)
    goodrows=[[]for _ in gs];blocks=[]
    for q0 in prime_range(3,2001):
        q=int(q0)
        if disc%q==0:continue
        roots=f.change_ring(GF(q)).roots(multiplicities=False)
        if not roots:continue
        vals=[[int(g.change_ring(GF(q))(root))for root in roots]for g in polys]
        if any(v==0 for row0 in vals for v in row0):continue
        for destrow,vs in zip(goodrows,vals):destrow.extend(int(pow(v,(q-1)//2,q)==q-1)for v in vs)
        blocks.append({'prime':q,'roots':list(map(int,roots))})
        if matrix(GF(2),goodrows).rank()==n:break
    rank=int(matrix(GF(2),goodrows).rank())
    if rank!=n:
        save(dest/'result.json',{**preliminary,'status':'UNKNOWN_KUMMER_RANK','detected':rank,'expected':n});return
    localrows=[[]for _ in gs];local=[];valuation=[];unram=[]
    for q in S:
        L=LocalSquareclasses(nf,q);signatures=[list(map(int,L.signature(g)))for g in gs]
        for agg,sig in zip(localrows,signatures):agg.extend(sig)
        local.append({'place':q,'point_image_dimension':int(L.point_kummer_dimension),
                      'generic_image_rank':int(matrix(GF(2),signatures[:17]).rank()),
                      'known_image_rank':int(matrix(GF(2),signatures).rank()),'signatures':signatures})
        for j,P in enumerate(pari.idealprimedec(nf,q)):
            vs=[int(pari.idealval(nf,g,P))for g in gs];valuation.append([v%2 for v in vs]);unram.append([v%2 for v in vs])
            if q!=2:continue
            # Retain a complete local unit group via PARI's prime-power bid.
            ram=int(P[2]);degree=int(P[3])
            if degree!=1:
                save(dest/'result.json',{**preliminary,'status':'UNKNOWN_DYADIC_RESIDUE_DEGREE','degree':degree});return
            pi=pari(2)if ram==1 else pari.nfbasistoalg(nf,pari.idealappr(nf,P))
            assert pari.idealval(nf,pi,P)==1
            for k in range(1,2*ram+1):unram.append([int(pari.nfhilbert(nf,g,1+pi**k,P)==-1)for g in gs])
        save(dest/('local-%s.json'%q),local[-1])
    realroots=f.roots(AA,multiplicities=False)
    real=[[int(g(root)<0)for root in realroots]for g in polys]
    for agg,sig in zip(localrows,real):agg.extend(sig)
    unram.extend([list(row0)for row0 in zip(*real)])
    joint=matrix(GF(2),localrows);v=matrix(GF(2),valuation);u=matrix(GF(2),unram)
    dimensions={'generic':17,'known':n,'generic_local':int(joint[:17,:].rank()),'known_local':int(joint.rank()),
                'generic_strict':17-int(joint[:17,:].rank()),'known_strict':n-int(joint.rank()),
                'generic_even':17-int(v[:,:17].rank()),'known_even':n-int(v.rank()),
                'generic_unramified':17-int(u[:,:17].rank()),'known_unramified':n-int(u.rank()),
                'local_product':sum(x['point_image_dimension']for x in local)+(1 if len(realroots)==3 else 0)}
    unit_bound=(3+len(realroots))//2-1
    dimensions['relative_ideal_lower_bound']=max(0,dimensions['known_even']-dimensions['generic_even']-unit_bound)
    save(dest/'result.json',{**preliminary,'status':'PASS_EXACT_POINT_SUBSPACE_ANATOMY','stage':'complete',
                            'field_discriminant':str(nf.disc()),'maximal_order_basis':list(map(str,nf.nf_get_zk())),
                            'signature':[len(realroots),(3-len(realroots))//2],'dimensions':dimensions,
                            'good_character_blocks':blocks,'local':local,'real_signatures':real,
                            'valuation_matrix':valuation,'unramified_constraint_matrix':unram,
                            'unit_dimension_upper_bound':unit_bound,
                            'strict_words':[list(map(int,w))for w in joint.left_kernel().basis()],
                            'generic_strict_words':[list(map(int,w))for w in joint[:17,:].left_kernel().basis()],
                            'full_class_group':'UNKNOWN','full_Selmer':'UNKNOWN',
                            'boundary':'Exact known-point squareclass/local dimensions; ideal bounds use Dirichlet only. Complete class/unit decomposition and prospective predictor validation are not asserted.'})
    print(case,dimensions,flush=True)
def run():
    p=read(PROTOCOL);results=[]
    for row in p['rows']:
        case=row['id'];folder=WORK/case;folder.mkdir(parents=True,exist_ok=True)
        path=folder/'result.json'
        if not path.exists():
            start=time.monotonic()
            with (folder/'worker.log').open('x')as log:
                try:
                    r=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--case',case],stdout=log,stderr=log,timeout=60)
                    failure=None if r.returncode==0 else 'WORKER_FAILED'
                except subprocess.TimeoutExpired:failure='TIME_LIMIT'
            save(folder/'runtime.json',{'seconds':time.monotonic()-start,'failure':failure,'log_sha256':sha(folder/'worker.log')})
            if failure and not path.exists():save(path,{'id':case,'status':'UNKNOWN_'+failure,'protocol_sha256':sha(PROTOCOL)})
        result=read(path);assert result['protocol_sha256']==sha(PROTOCOL);results.append(result)
        print('CASE',case,result['status'],result.get('dimensions'),flush=True)
    save(OUT,{'schema':'rank-jump.matched-descent-panel.v1','protocol_sha256':sha(PROTOCOL),'cases':results,
              'complete_cases':sum(r['status']=='PASS_EXACT_POINT_SUBSPACE_ANATOMY'for r in results)})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','run','worker']);p.add_argument('--case');args=p.parse_args()
    if args.mode=='freeze':freeze()
    elif args.mode=='run':run()
    else:worker(args.case)
