"""Bounded arithmetic phases for the detached foundry; no network or model client."""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path

from v3_warm_support import read, sha, require, atomic as _atomic
from high_rank_foundry_policy import parent_masks, hashed

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def save(path, value, immutable=True):
    _atomic(path,json.loads(json.dumps(value)),immutable=immutable)


def verify_packet(packet, row):
    from r17_60_arithmetic import native_check
    native_check(packet,row)
    require(packet['curve']==row['model'],'packet equation differs from intake')


def prepare_bank(row, packet_path, output, bank_index, job):
    """Independent integer/rational generic CVPs; full-space lower-shell refreshes."""
    from sage.all import ZZ, matrix, pari
    import numpy as np
    from visibility_lattice_fast import IntegerExactParity
    from visibility_lattice_v2 import ExactParity
    packet = read(packet_path); verify_packet(packet,row)
    cat = ART/'r17_exact_maximum_parity_classes_v1.json'
    family = next(f for f in read(cat)['families'] if f['family']==row['family'])
    maximum = [r['mask'] for r in family['classes']]
    masks,mode = parent_masks(row['family'],maximum,bank_index)
    g = 2*matrix(ZZ,family['gram'])
    u = matrix(ZZ,pari(g).qflllgram()).transpose(); inv = u.inverse()
    require(abs(u.det())==1,'generic LLL not unimodular')
    fast = IntegerExactParity((u*g*u.transpose()).rows())
    ref = ExactParity((u*g*u.transpose()).rows())
    candidates, checks = [],[]
    for mask in masks:
        w = matrix(ZZ,1,17,[(mask>>j)&1 for j in range(17)])
        residue = tuple(int(x)%2 for x in (w*inv).row(0))
        starts,_ = fast.babai(np.asarray([residue],dtype=np.int64)); start = tuple(map(int,starts[0]))
        proof = fast.solve(residue,start,2000000)
        require(proof==ref.solve(residue,start,2000000),'independent generic CVP differs')
        if mode=='maximum':
            require(proof['norm']==2*family['exact_maximum_parity_minimum'],'maximum norm differs')
        words = set()
        for v in proof['minima']:
            z = tuple(map(int,(matrix(ZZ,1,17,v)*u).row(0)))
            if next(x for x in z if x)<0: z=tuple(-x for x in z)
            require(sum((x%2)<<j for j,x in enumerate(z))==mask,'generic parity transport differs')
            words.add(z)
        words = sorted(words)
        # Retain all exact minimal alternatives, choose reproducibly for this bank.
        index = int(hashed(row['family'],bank_index,mask),16)%len(words)
        candidates.append({'mask':mask,'word':list(words[index]),'norm':proof['norm'],
                           'alternative_words':[list(w) for w in words]})
        checks.append({'mask':mask,'reduced_seed':start,'proof':proof})
    pivots = {}
    def insert(m):
        for k in sorted(pivots,reverse=True):
            if (m>>k)&1: m ^= pivots[k]
        if m: pivots[m.bit_length()-1]=m
        return bool(m)
    if mode != 'maximum':
        for m in maximum: insert(m)
    candidates.sort(key=lambda r:(-r['norm'],r['mask']))
    chosen = []
    for c in candidates:
        if insert(c['mask']): chosen.append(c)
        if len(chosen)==16: break
    for c in candidates:
        if len(chosen)==16: break
        if c not in chosen: chosen.append(c)
    require(bool(chosen),'empty parent bank')
    protocol = {'family':row['family'],'parameter':row['parameter'],'initial_rank':packet['rank_lower_bound'],
                'generic_rank':17,'bank_index':bank_index,'point_searches':0,
                'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (packet_path,cat,job/'request.json',Path(__file__))}}
    save(output/'protocol.json',protocol)
    save(output/f"seed-M{packet['rank_lower_bound']}.json",packet)
    save(output/'generic-cvp-proofs.json',{'gram':[list(map(int,r)) for r in g.rows()],
          'LLL':[list(map(int,r)) for r in u.rows()],'checks':checks})
    save(output/'derivation.json',{'mode':mode,'bank_index':bank_index,'selected_masks':[c['mask'] for c in chosen],
         'covered_span_rank':len(pivots),'all_candidates':candidates,
         'rule':'Disjoint exact maximum blocks; then affine-permuted full17 parity samples, extending the maximum-class span before filling by norm. All minima and both CVP proofs retained.'})
    save(output/'anchor-bank.json',{'status':'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET','dimension':17,
         'generic_gram_scale':2,'rows':[{k:c[k] for k in ('mask','word','norm')} for c in chosen],
         'shells':sorted({c['norm'] for c in chosen}), 'protocol_sha256':sha(output/'protocol.json'),
         'derivation_sha256':sha(output/'derivation.json'),'generic_cvp_sha256':sha(output/'generic-cvp-proofs.json')})
    save(output/'prepared.json',{'status':'PASS_FOUNDRY_PARENT_PREPARATION','point_searches':0,
         'files':{p.name:sha(p) for p in sorted(output.glob('*.json')) if p.name!='prepared.json'}})


def phase(job, name, index=0):
    request = read(job/'request.json'); row = request['candidate']; case = row['id']
    config = request['config']
    if name in ('search','seed-search'):
        import pari_pointed_backend as backend
        from high_rank_foundry_point_calls import install
        install(job,backend)
    import r17_60_arithmetic as panel
    panel.D = job
    folder = job/'cases'/case
    folder.mkdir(parents=True,exist_ok=True)
    if name == 'generic':
        from high_rank_foundry_seed_gate import assess
        gate=assess(row)
        save(job/'generic-seed-gate.json',gate)
        if gate['status']=='UNRESOLVED_GENERIC_SEED':
            return
        save(job/'roster.json',{'rows':[row]})
        save(job/'protocol.json',{'schema':'foundry-generic.v1','request_sha256':sha(job/'request.json')})
        panel.generic()
    elif name == 'seed-prepare':
        import run_fresh6_seed_confirmation_v2 as seed
        dest = folder/'seed-search'
        # A crash during preparation is retained. It has no completed point call.
        if dest.exists() and not (dest/'geometry.json').exists():
            dest.rename(dest.with_name('incomplete-seed-preparation-'+str(__import__('time').time_ns())))
        if not dest.exists():
            seed.prepare(job/'generic'/case/'seed-M17.json',dest)
        p = read(dest/'protocol.json')
        p['sources'][str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
        p['height']=config['height']; p['seconds_per_chart']=config['point_seconds']; p['seconds_per_map']=config['map_seconds']
        # Cooled M17 fibres receive genuinely different generic classes on revival.
        if request.get('revival'):
            family = next(f for f in read(seed.PARITY)['families'] if f['family']==row['family'])
            masks,_ = parent_masks(row['family'],[r['mask'] for r in family['classes']],request['bank_index'])
            p['masks']=masks[:16]; p['max_point_invocations']=2*len(p['masks'])
            p['rule']='Revival on a deterministic full17 parity sample; bounded seed search, no rank exclusion.'
            save(dest/'geometry.json',seed.geometry(read(dest/'seed.json'),p['masks']),immutable=False)
            p['inputs'][str((dest/'geometry.json').relative_to(ROOT))]=sha(dest/'geometry.json')
        save(dest/'protocol.json',p,immutable=False)
    elif name in ('seed-search','seed-replay'):
        import run_fresh6_seed_confirmation_v2 as seed
        seed.execute(folder/'seed-search',name=='seed-replay')
    elif name == 'seed-cloud':
        panel.reconcile_seed(case)
    elif name == 'bank':
        spec = read(job/f'segment-{index:02d}.json')
        prepare_bank(row,ROOT/spec['packet'],job/f'preparation-{index:02d}',request['bank_index'],job)
    elif name in ('search','replay'):
        import run_complement_seed_v3 as runner
        spec = read(job/f'segment-{index:02d}.json'); dest=job/f'search-{index:02d}'
        if dest.exists() and not (dest/'protocol.json').exists():
            dest.rename(dest.with_name(dest.name+'-incomplete-freeze-'+str(__import__('time').time_ns())))
        if spec.get('parent'):
            import run_complement_cached_seed_v3 as cached
            cached.SOURCE=ROOT/spec['parent']; runner.obtain_map=cached.obtain_map
            if not dest.exists():
                p=cached.freeze(dest)
                p['max_charts']=p['inherited_charts']+spec['allowance']
                p['additional_chart_budget']=spec['allowance']
                save(dest/'protocol.json',p,immutable=False)
        else:
            runner.PREP=runner.BANK=job/f'preparation-{index:02d}'
            if not dest.exists():
                p=runner.freeze(dest);p['max_charts']=spec['allowance']
                p['height']=config['height'];p['seconds_per_chart']=config['point_seconds'];p['seconds_per_map']=config['map_seconds']
                p['sources'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
                p['inputs'][str((job/'request.json').relative_to(ROOT))]=sha(job/'request.json')
                save(dest/'protocol.json',p,immutable=False)
        (runner.search if name=='search' else runner.replay)(dest)
    elif name == 'cloud':
        from reconcile_verified_v3_cloud import run
        run(job/f'search-{index:02d}',job/f'reconciled-{index:02d}.json')
    elif name == 'verify-packet':
        packet=read(job/'packet.json');verify_packet(packet,row)
        # Separate older implementation: recomputes finite quotients without the
        # producer's ephemeral ReductionCache. It certifies the same lower bound.
        import certify_compact_r17_candidates as cert
        points=tuple(tuple(map(F,p)) for p in packet['points']);model=tuple(map(F,packet['curve']))
        proof=packet['proof']
        other=cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        require(other['rank_lower_bound']==packet['rank_lower_bound'],'independent portable certificate failed')
        save(job/'packet-verified.json',{'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','packet_sha256':sha(job/'packet.json'),
             'rank_lower_bound':packet['rank_lower_bound'],'point_searches':0})
    elif name == 'conductor-build' or name=='conductor-replay':
        from inventory_conductor_worker_v2 import run
        run(job/'conductor-input.json',job/'conductor.json',name=='conductor-replay')
    else:
        raise ValueError(name)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase');p.add_argument('--job',type=Path,required=True);p.add_argument('--index',type=int,default=0)
    a=p.parse_args();phase(a.job,a.phase,a.index)
