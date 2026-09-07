#!/usr/bin/env python3
"""Fixed 302 generic-only control, full31 target and two certified19 siblings."""
import argparse, sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
BATCH=LOCAL/'curve302-focused-point-exposure-v1'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['curve302_focused_point_exposure.py','prepare_curve302_focused_maps.sage',
           'verify_curve302_focused_seeds.sage','verify_factor_free_rank.sage','run_curve302_focused_point_exposure.py']
    return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    assert not (BATCH/'intake.json').exists()
    recovery=LOCAL/'mixed-reduced-parent-relocation-replay-v1/ledger.json'
    assert cert.read(recovery)['status']=='PASS'
    source=ART/'curve302_recovered_mw17_parent_v1.json';parent=cert.read(source)
    # Only generic formulas go to the masked control. The public embedding and
    # point words in the parent artifact are not used for its seed or geometry.
    redacted={k:parent[k] for k in ['a_invariants','basis_weierstrass_coordinates']}
    checkpoint(BATCH/'generic-parent.json',redacted)
    value=lambda q:cert.F(q['numerator'][0])/cert.F(q['denominator'][0])
    a=tuple(value(q) for q in redacted['a_invariants']);assert a[:3]==(1,1,1)
    b2=a[0]**2+4*a[1];b4=2*a[3]+a[0]*a[2];b6=a[2]**2+4*a[4]
    c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
    model=(cert.F(0),cert.F(0),cert.F(0),-27*c4,-54*c6)
    generic_long=[tuple(value(q) for q in P) for P in redacted['basis_weierstrass_coordinates']]
    generic=[(36*x+3*b2,108*(2*y+a[0]*x+a[2])) for x,y in generic_long]
    assert len(generic)==17 and all(cert.is_on_weierstrass_curve(a,P) for P in generic_long)
    rows=[];paths=[source,recovery,ART/'blind_factor_free_28_control_v1.json',CAS/'icarm_curve302.py']
    def add(ident,model,points,old_generic,floor,role,parameter):
        from research_runtime.search_state import raw_state
        from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
        from research_runtime.memory_store import MemoryFactStore
        pts=tuple(tuple(map(cert.F,P)) for P in points);model=tuple(map(cert.F,model))
        state=raw_state(model,pts,cache=QuotientOnlyReductionCache(MemoryFactStore()),prime_bound=1000)
        assert tuple(state.basis)==pts
        proof=checked_rank(model,pts,state.reductions.primes,state.no_two_torsion_prime)
        seed=dict(family='det1092',parameter=parameter,curve=list(map(str,model)),
            points=[list(map(str,P)) for P in pts],generic_points=[list(map(str,P)) for P in old_generic],rank_certificate=proof)
        out=BATCH/ident/'seed.json';assert not out.exists();checkpoint(out,seed)
        # A copied-input complete finite-group verifier also checks every seed.
        checkpoint(BATCH/ident/'seed-cloud.json',dict(curve=seed['curve'],points=seed['points'],
            signatures=proof['signatures'],rank_certificate=proof,rank_lower_bound=len(pts),
            independent_column_indices=list(range(len(pts)))))
        rows.append(dict(id=ident,initial_rank=len(pts),generic_dimension=len(old_generic),mask_floor=floor,
            family='det1092',parameter=parameter,role=role,seed_sha256=cert.hashed(out)))
    add('curve302-generic17',model,generic,generic,0,'masked recovery control; only parent sections', '0')
    # This import occurs only in intake, never in a point or geometry worker.
    import icarm_curve302 as public
    assert model==public.short_coefficients()
    add('curve302-full31',model,public.SHORT_POINTS,[],0,'direct32 search from public31; empty generic_points marks no prefix alignment, not generic rank0','0')
    for ident in ['d1092-fibre-018','d1092-fibre-050']:
        old=LOCAL/'mixed-reduced-parent-exposure-v1'/ident
        cloudpath=ART/('mixed_reduced_parent_'+ident.replace('-','_')+'_mod2_v1.json')
        cloud=cert.read(cloudpath);proof=cert.read(old/'certification-ledger.json');seed=cert.read(old/'seed.json')
        assert proof['status']=='PASS' and proof['rank_lower_bound']==19
        assert all(v==19 for v in proof['odd_modulus_ranks'].values())
        assert cloud['independent_points'][:17]==seed['points'] and len(seed['points'])==17
        paths += [cloudpath,old/'certification-ledger.json',old/'seed.json']
        add(ident,seed['curve'],cloud['independent_points'],seed['points'],17,
            'updated own19 subgroup; centre parity outside previous17',seed['parameter'])
    assert [r['initial_rank'] for r in rows]==[17,31,19,19]
    checkpoint(BATCH/'intake.json',dict(status='PASS',rows=rows,sources=sources(),
        inputs={str(x.relative_to(ROOT)):cert.hashed(x) for x in paths},generic_parent_sha256=cert.hashed(BATCH/'generic-parent.json')))

def freeze():
    assert not (BATCH/'protocol.json').exists()
    intake=cert.read(BATCH/'intake.json');assert intake['sources']==sources()
    assert cert.read(BATCH/'seed-replay/result.json')['status']=='PASS'
    fresh=BATCH/'seed-replay';replay=cert.read(fresh/'protocol.json')
    assert all(cert.hashed(fresh/n)==h for n,h in replay['files'].items())
    assert all(cert.hashed(fresh/(r['id']+'-seed.json'))==r['seed_sha256']==cert.hashed(BATCH/r['id']/'seed.json') for r in intake['rows'])
    paths=[BATCH/'intake.json',BATCH/'generic-parent.json',fresh/'result.json',fresh/'protocol.json']
    paths += [BATCH/r['id']/'seed.json' for r in intake['rows']]
    checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.curve302-focused-point-exposure.v1',
        sources=sources(),inputs={**intake['inputs'],**{str(x.relative_to(ROOT)):cert.hashed(x) for x in paths}},
        rows=intake['rows'],sample_size=2048,sample_domain='full11952-specialized-followup-v1',charts=49,
        height=125000,seconds_per_chart=10,rank_stop=False,target_rank=32,maximum_point_boxes=196,
        geometry_wall_seconds=180,worker_wall_seconds=1200,replay_wall_seconds=1200,
        rss_bytes=2147483648,maximum_workers=1,gp_sha256=cert.hashed(Path('/usr/bin/gp')),
        gate='Full original27-only blind28 control passed. All441 mixed-parent point boxes and histories completed; independent certificates completed by relocation continuation. Two siblings each have certified19 directions. Parent generic17 formulas specialize exactly to302. Seed independence independently checked.',
        selection='User-requested302 control and direct target; precisely the two highest certified sibling bounds from the frozen six-fibre comparison, one strong and one moderate score stratum. No new parameters, validation primes, public point input to sibling selection, or control-outcome-driven replacement.',
        centre_policy='Same2048 SHA parities,384-bit rounded metric and49 largest computed norms as calibrated factor-free method. Control17 and full31 use all nonzero parity classes of their respective supplied subgroup. Siblings require nonzero coefficients above inherited17. All196 maps freeze before any point worker. No CVP optimality or coverage theorem.',
        endpoints='Control recovery relative to generic17 is recorded separately from any direct32 discovery and sibling rank gains. A null full31 search supplies no upper rank bound. Parent selection is explicitly informed by302; sibling addresses and point search remain free of its exceptional points.',
        failure_policy='Preserve any failure/censoring; no refill, rank stop, automatic larger box or subsequent wave.',following_campaign=None))
    print('FROZEN196 boxes:302 control17, target31, siblings19/19',flush=True)

def protocol():
    p=cert.read(BATCH/'protocol.json')
    assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items())
    return p

def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
        if m>>ROW['mask_floor'] and m not in result:result.append(m)
    return result

def configure(index):
    global ROW,D,SEED
    ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'
    engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def launch():
    p=protocol();out=BATCH/'ledger.json';assert not out.exists()
    ledger=dict(status='RUNNING_GEOMETRY',maps=[],rows=[]);checkpoint(out,ledger)
    for i,row in enumerate(p['rows']):
        configure(i)
        s=run([SAGE,str(CAS/'prepare_curve302_focused_maps.sage'),'--index',str(i)],
            limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['maps'].append(dict(id=row['id'],status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
        print(row['id'],'maps',s['outcome'],s['returncode'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('all196 maps required before points')
    ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
    for i,row in enumerate(p['rows']):
        configure(i);entry=dict(id=row['id'],status='RUNNING',stages=[]);ledger['rows'].append(entry);checkpoint(out,ledger)
        for name in ['worker','replay']:
            s=run([sys.executable,str(Path(__file__).resolve()),name,'--index',str(i)],
                limits=Limits(p[name+'_wall_seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            entry['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
            print(row['id'],name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
            if not ok:
                ledger['status']=entry['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed/censored stage')
        entry.update(status='PASS',rank_lower_bound=cert.read(D/'result.json')['rank_lower_bound']);checkpoint(out,ledger)
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','freeze','launch','worker','replay']);ap.add_argument('--index',type=int);a=ap.parse_args()
    if a.stage in ['worker','replay']:configure(a.index);getattr(engine,a.stage)()
    else:globals()[a.stage]()
