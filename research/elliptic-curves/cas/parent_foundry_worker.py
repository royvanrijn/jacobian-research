"""Generic-rank-aware point evaluator using the replayed V3 search machinery."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from v3_warm_support import read, sha, require, atomic

CAS=Path(__file__).resolve().parent; ROOT=CAS.parents[1]


def save(path, data):
    data=json.loads(json.dumps(data))
    if path.exists():require(read(path)==data,'immutable checkpoint changed: '+str(path))
    else:atomic(path,data,immutable=True)


def specialize(parent,parameter):
    from compact_atlas_specialization import polynomial,rational
    from certify_compact_r17_candidates import is_on_weierstrass_curve
    t=F(parameter);d=t.denominator
    a=polynomial(parent['A_coefficients_low_to_high'],t)*d**8
    b=polynomial(parent['B_coefficients_low_to_high'],t)*d**12
    if 4*a**3+27*b**2==0:raise ValueError('SINGULAR_FIBRE')
    model=(F(0),F(0),F(0),a,b)
    points=tuple((rational(s['X'],t)*d**4,rational(s['Y'],t)*d**6) for s in parent['sections'])
    require(len(points)==parent['generic_rank_lower_bound'] and all(is_on_weierstrass_curve(model,p) for p in points),
            'native generic specialization failed')
    return model,points


def verify(packet,parent,parameter):
    from v3_warm_engine import certified_state
    import certify_compact_r17_candidates as independent
    model,points=specialize(parent,parameter)
    require(packet['curve']==list(map(str,model)) and packet['points'][:len(points)]==[list(map(str,p)) for p in points],
            'generic point prefix or equation changed')
    certified_state(packet['curve'],packet['points'],packet['proof'])
    proof=packet['proof']
    other=independent.checked_rank(model,tuple(tuple(map(F,p)) for p in packet['points']),
                                   [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    require(other['rank_lower_bound']==packet['rank_lower_bound'],'second finite implementation disagrees')


def generic_seed(parent,parameter,job):
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import short_curve_has_no_rational_2_torsion_modular_certificate
    try:model,points=specialize(parent,parameter)
    except (ValueError,ZeroDivisionError) as e:
        save(job/'seed-gate.json',{'status':'UNRESOLVED_GENERIC_SPECIALIZATION','reason':str(e)})
        return None
    admission=FinitePointAdmission(model,(),prime_bound=1000)
    outcomes=[admission.consider(p) for p in points]
    torsion=next((p for p in admission.primes if p<=200 and
        short_curve_has_no_rational_2_torsion_modular_certificate(model,p)),None)
    full=len(admission.points)==len(points) and torsion is not None
    save(job/'seed-gate.json',{'status':'PASS_GENERIC_SEED_GATE' if full else 'UNRESOLVED_GENERIC_SPECIALIZATION',
        'required':len(points),'admitted':len(admission.points),'torsion_prime':torsion,
        'curve':list(map(str,model)),'points':[list(map(str,p)) for p in points],
        'outcomes':outcomes,'boundary':'A bounded finite witness failure is UNKNOWN, not dependence or an upper bound.'})
    if not full:return None
    proof=checked_rank(model,points,admission.primes,torsion)
    packet={'family':parent['family'],'parameter':parameter,'curve':list(map(str,model)),
        'points':[list(map(str,p)) for p in points],'proof':proof,
        'generic_rank':len(points),'rank_lower_bound':len(points)}
    verify(packet,parent,parameter);save(job/'generic.json',packet)
    return packet


def bank(parent,packet,folder,index):
    from sage.all import matrix,ZZ,QQ,pari
    import numpy as np
    from visibility_lattice_fast import IntegerExactParity
    from visibility_lattice_v2 import ExactParity
    rank=parent['generic_rank_lower_bound'];g=2*matrix(QQ,parent['generic_height_gram'])
    g=matrix(ZZ,g);u=matrix(ZZ,pari(g).qflllgram()).transpose();inv=u.inverse()
    require(abs(u.det())==1,'nonunimodular generic lattice reduction')
    fast=IntegerExactParity((u*g*u.transpose()).rows());reference=ExactParity((u*g*u.transpose()).rows())
    modulus=1<<rank
    seed=int(hashlib.sha256(('parent-banks-v1/'+parent['family']).encode()).hexdigest(),16)
    multiplier=(seed%modulus)|1;offset=(seed>>rank)%modulus
    rows=[];proofs=[]
    for k in range(index*64,(index+1)*64):
        mask=(multiplier*k+offset)%modulus
        if not mask:continue
        w=matrix(ZZ,1,rank,[(mask>>i)&1 for i in range(rank)])
        residue=tuple(int(x)%2 for x in (w*inv).row(0))
        starts,_=fast.babai(np.asarray([residue],dtype=np.int64));start=tuple(map(int,starts[0]))
        p=fast.solve(residue,start,2000000)
        require(p==reference.solve(residue,start,2000000),'generic exact CVP implementations disagree')
        words=set()
        for v in p['minima']:
            word=tuple(map(int,(matrix(ZZ,1,rank,v)*u).row(0)))
            if next(x for x in word if x)<0:word=tuple(-x for x in word)
            require(sum((x%2)<<i for i,x in enumerate(word))==mask,'generic parity transport differs')
            words.add(word)
        words=sorted(words);word=words[(index+mask)%len(words)]
        rows.append({'mask':mask,'word':list(word),'norm':p['norm'],'alternatives':[list(w) for w in words]})
        proofs.append({'mask':mask,'start':start,'proof':p})
    rows.sort(key=lambda r:(-r['norm'],r['mask']))
    # Twelve deep classes plus four full-space exploration classes. Never restrict
    # later banks to the linear span of previously productive generic classes.
    chosen=rows[:12]
    rest=sorted(rows[12:],key=lambda r:hashlib.sha256(f'{index}/{r["mask"]}'.encode()).hexdigest())
    chosen+=rest[:4]
    save(folder/'protocol.json',{'family':parent['family'],'parameter':packet['parameter'],
        'generic_rank':rank,'bank_index':index,'rule':'64 full-space parity proposals; 12 deep and 4 exploratory, all minimal alternatives retained.'})
    save(folder/'derivation.json',{'candidates':rows,'selected':[r['mask'] for r in chosen]})
    save(folder/'generic-cvp-proofs.json',{'gram':[list(map(int,r)) for r in g.rows()],
        'LLL':[list(map(int,r)) for r in u.rows()],'proofs':proofs})
    save(folder/'anchor-bank.json',{'status':'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET','dimension':rank,
        'generic_gram_scale':2,'rows':[{k:r[k] for k in ('mask','word','norm')} for r in chosen],
        'shells':sorted({r['norm'] for r in chosen}),
        'protocol_sha256':sha(folder/'protocol.json'),'derivation_sha256':sha(folder/'derivation.json')})
    save(folder/f'seed-M{packet["rank_lower_bound"]}.json',packet)
    save(folder/'prepared.json',{'status':'PASS_PARENT_FOUNDRY_BANK','files':{
        p.name:sha(p) for p in sorted(folder.glob('*.json')) if p.name!='prepared.json'}})


def evaluate(job):
    request=read(job/'request.json');parent=read(ROOT/request['parent'])
    require(sha(ROOT/request['parent'])==request['parent_sha256'],'parent input changed')
    parameter=request['parameter'];generic=parent['generic_rank_lower_bound']
    if request.get('packet'):
        packet=read(ROOT/request['packet']);require(sha(ROOT/request['packet'])==request['packet_sha256'],'continuation changed')
        packet.update(family=parent['family'],parameter=parameter)
        verify(packet,parent,parameter)
    else:packet=generic_seed(parent,parameter,job)
    if packet is None:
        save(job/'result.json',{'status':'UNRESOLVED_GENERIC_SPECIALIZATION','rank_lower_bound':None,
            'calls':0,'request_sha256':sha(job/'request.json')});return
    initial=packet['rank_lower_bound'];used=0;segments=[];events=[]
    import run_complement_seed_v3 as runner
    from parent_foundry_seed import freeze as freeze_search
    import pari_pointed_backend as backend
    from high_rank_foundry_point_calls import install
    from reconcile_verified_v3_cloud import run as reconcile
    from high_rank_foundry_certificate import compact_packet
    install(job,backend)
    for segment in range(8):
        if used>=request['allowance'] or packet['rank_lower_bound']>=32:break
        prep=job/f'bank-{segment:02d}';dest=job/f'search-{segment:02d}'
        if not (prep/'prepared.json').exists():bank(parent,packet,prep,request['bank_index']*8+segment)
        runner.PREP=runner.BANK=prep
        if not dest.exists():
            p=freeze_search(dest,prep)
            p.update(max_charts=request['allowance']-used,height=request['height'],
                     seconds_per_chart=10,seconds_per_map=5)
            p['inputs'][str((job/'request.json').relative_to(ROOT))]=sha(job/'request.json')
            p['sources'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
            atomic(dest/'protocol.json',p)
        runner.search(dest);runner.replay(dest)
        reconciled=job/f'reconciled-{segment:02d}.json'
        reconcile(dest,reconciled);reconcile(dest,reconciled)
        terminal=read(dest/'terminal.json');next_packet=read(reconciled)
        n=terminal['charts'];require(0<=n<=request['allowance']-used,'point allowance exceeded')
        offsets={};at=used
        for stage in terminal['stages']:
            offsets[stage['epoch']]=at;at+=stage['charts']
            if stage['after']>stage['before']:
                events.append({'before':stage['before'],'after':stage['after'],'call':at,'phase':'adaptive'})
        for gain in next_packet['gains']:
            chart=Path(gain['chart']);epoch=int(chart.parts[0].split('-')[1]);i=int(chart.stem.split('-')[1])
            events.append({'before':gain['after']-1,'after':gain['after'],
                           'call':offsets[epoch]+i+1,'phase':'full-cloud'})
        used+=n
        segments.append({'search':str(dest.relative_to(ROOT)),'calls':n,
            'before':packet['rank_lower_bound'],'after':next_packet['rank_lower_bound'],
            'stop':terminal['stop_reason'],'point_timeouts':sum(s['censored'] for s in terminal['stages']),
            'map_timeouts':sum(s['censored_maps'] for s in terminal['stages'])})
        packet=next_packet;packet.update(family=parent['family'],parameter=parameter)
        require(packet['rank_lower_bound']>=max([terminal['rank_lower_bound']]+[
            r for s in terminal['stages'] for r in s['odd_ranks'].values()]),'unresolved stronger odd-prime cloud')
        if n==0:break
    packet=compact_packet(packet);verify(packet,parent,parameter)
    packet['status']='PASS_TWO_FINITE_IMPLEMENTATIONS'
    save(job/'packet.json',packet)
    save(job/'verified.json',{'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','packet_sha256':sha(job/'packet.json'),
        'generic_prefix_sha256':sha(job/'generic.json') if (job/'generic.json').exists() else request['packet_sha256']})
    save(job/'result.json',{'status':'PASS_CERTIFIED_PARENT_EVALUATION','family':parent['family'],
        'parameter':parameter,'generic_rank':generic,'initial_rank':initial,'rank_lower_bound':packet['rank_lower_bound'],
        'certified_jump_lower_bound':packet['rank_lower_bound']-generic,'calls':used,
        'allowance':request['allowance'],'call_allowance_consumed':used==request['allowance'],
        'exposure_complete':packet['rank_lower_bound']>=32 or (used==request['allowance'] and
            not any(s['point_timeouts'] or s['map_timeouts'] for s in segments)),
        'segments':segments,'gain_timeline':events,'packet':str((job/'packet.json').relative_to(ROOT)),
        'packet_sha256':sha(job/'packet.json'),'request_sha256':sha(job/'request.json'),
        'claim_boundary':'Certified rank lower bound, not exact rank. Jump is a certified lower bound. Empirical tail rates describe the fixed detection protocol; censored exposures stay explicit.'})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--job',type=Path,required=True)
    parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
    args.job=args.job.resolve()
    if args.verify_only:
        r=read(args.job/'request.json');verify(read(args.job/'packet.json'),read(ROOT/r['parent']),r['parameter'])
    else:evaluate(args.job)
