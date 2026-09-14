#!/usr/bin/env sage -python
"""Four cold, generic-only fresh fibres after the independent transfer gate."""
import argparse
import ast
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import signal
import tempfile
import time

import cancellation_scheduler_round4 as round4
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, OUT as CORPUS, canonical, digest, write

CAS=Path(__file__).resolve().parent
OUT=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_fresh_v1'
RAW=LOCAL/'cancellation-scheduler-fresh-v1'
PRIOR=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_transfer_v1'
FOUNDRY=ROOT/'artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1'
RUNTIME=LOCAL/'x948-seed-foundry-v1/runtime/research'


def load(path):return json.loads(path.read_text())
def sha(path):return digest(path.read_bytes())
def require(test,message):
    if not test:raise ArithmeticError(message)


def source_closure(start):
    """Conservatively bind local imports and literal script dependencies."""
    pending=list(start);seen=set()
    while pending:
        p=pending.pop().resolve()
        if p in seen or not p.is_file():continue
        seen.add(p)
        if p.suffix not in ('.py','.sage'):continue
        tree=ast.parse(p.read_text())
        for node in ast.walk(tree):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):
                names=[node.module or '']
                if node.module:names += [node.module+'.'+a.name for a in node.names]
            for name in names:
                for base in (CAS,p.parent):
                    stem=base.joinpath(*name.split('.'))
                    pending += [stem.with_suffix('.py'),stem/'__init__.py']
            if isinstance(node,ast.Constant) and isinstance(node.value,str):
                value=Path(node.value)
                if value.suffix in ('.py','.sage') and len(node.value)<200:
                    pending += [CAS/value,p.parent/value,ROOT/value]
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(seen)}


def prepare():
    require(load(PRIOR/'summary.json')['fresh_fibre_gate'] is True,'high-rank transfer gate failed')
    prior=load(PRIOR/'protocol.json');roster=load(FOUNDRY/'roster.json')
    oldplan=load(FOUNDRY/'plan.json');exposure=load(FOUNDRY/'exposure-state.json')
    # Fix the address before reading specialized equations, ranks or j values.
    ordinal=64;parameter=oldplan['parameters'][ordinal]
    chosen=sorted((r for r in roster['parents'] if not r['baseline']),
        key=lambda r:digest(('cancellation-fresh-v1/'+r['family']).encode()))[:4]
    require(len(chosen)==4,'four retained parents required; no replacement')
    inputs=[];bindings={};audit=[]
    for row in chosen:
        parent=RUNTIME/row['path'];require(sha(parent)==row['sha256'],'parent bytes changed')
        old=[e for e in exposure['exposures'] if e['family']==row['family'] and e['parameter']==parameter]
        require(not old,'declared fresh address has prior exposure; no replacement')
        # Reject a checkpoint absent from the stopped exposure snapshot too.
        live=[]
        for request in (RUNTIME/'jobs').rglob('request.json'):
            r=load(request)
            if r.get('parameter')==parameter and (r.get('family')==row['family'] or
                Path(r.get('parent','')).name==parent.name):live.append(str(request.relative_to(ROOT)))
        require(not live,'unlisted prior address checkpoint exists')
        ident=digest(canonical(['cancellation-fresh-v1',row['family'],parameter]))[:20]
        inputs.append({'id':ident,'family':row['family'],'stratum':'fresh_generic_M16',
            'parameter':parameter,'address_index':ordinal,'parent':str(parent.relative_to(ROOT)),
            'parent_sha256':row['sha256']})
        bindings[str(parent.relative_to(ROOT))]=row['sha256']
        audit.append({'family':row['family'],'parameter':parameter,'prior_exposures':0,'prior_checkpoints':live})
    design={'status':'DESIGN_BEFORE_FRESH_SPECIALIZATION','arms':prior['arms'],
        'candidate':prior['candidate'],'cases':4,'heights':[8000,32000,125000],
        'maximum_centres':48,'maximum_calls':432,'search_cpu_seconds':40,
        'point_wall_seconds':5,'cold_wall_seconds':90,'arm_wall_seconds':300,
        'hard_process_cpu_seconds':240,'fit':prior['fit'],'gp_sha256':prior['gp_sha256'],
        'selection':'Four nonbaseline admitted X948 parents in SHA256(cancellation-fresh-v1/ + family) order, each at the old plan zero-based address 64. No specialized equation, j, rank, point or jump chooses the roster. No replacement. Original foundry remains stopped.',
        'cold_recipe':'Per arm, specialize only the parent generic sections, certify the initial subgroup, and independently solve the existing generic exact CVP proposals for banks 0, 1 and 2 (48 selected anchors from 192 proposals). Rebuild and charge all work. No old specialized point cloud or cached bank is loaded.',
        'search':'Unmodified V4 candidate, fitted before these addresses are specialized. Same 40 CPU seconds, H125000 ceiling and 432-call allowance. Stop search at first independently certified direction; retain complete failed and censored calls.',
        'cloud':'Without more point calls, reconcile every returned cloud with fixed finite places through 1000. Certify the complete resulting subgroup using portable arithmetic and independent Sage finite-group cosets; do not claim maximal rank or saturation.',
        'cost':'Outer isolated process CPU covers cold specialization and bank construction, first-direction search and its independent replay, full-cloud admission and its independent rank proof, interpreter and discarded work. Both arms reconstruct identical inputs. Offline fit and parent construction remain retained one-time development costs, not repeated per fibre.',
        'gate':'Descriptive paired pilot on four predetermined fibres, after passing the 24-j validation and all three high-rank controls. No population inference, automatic expansion, basis reuse after enlargement, or rank32 claim.',
        'freshness':'Before execution check parent/address exposure only. Compare exact j against the retained corpus only after all arms finish. Novelty is relative to retained corpus and exposure, not all mathematics.'}
    new_write(OUT/'design.json',design);new_write(OUT/'inputs.json',inputs)
    source_inputs={str(FOUNDRY/n):sha(FOUNDRY/n) for n in ('roster.json','plan.json','exposure-state.json')}
    new_write(OUT/'preflight.json',{'status':'PASS_UNSEARCHED_RECIPE','address_audit':audit,
        'parent_sha256':bindings,'historical_input_sha256':{str(Path(k).relative_to(ROOT)):v for k,v in source_inputs.items()},
        'prerequisite_summary_sha256':sha(PRIOR/'summary.json')})
    print(json.dumps({'status':'FROZEN_FRESH_RECIPES','inputs':inputs}),flush=True)


def seal():
    require(load(PRIOR/'summary.json')['fresh_fibre_gate'] is True,'transfer gate changed')
    prior=load(PRIOR/'protocol.json');plan=load(OUT/'design.json');pre=load(OUT/'preflight.json')
    sources=source_closure([ROOT/p for p in prior['source_sha256']]+[Path(__file__),CAS/'parent_foundry_worker.py'])
    for name,h in prior['source_sha256'].items():require(sources[name]==h,'frozen ancestor source changed')
    sources.update(pre['parent_sha256']);sources.update(pre['historical_input_sha256'])
    sources[str(PRIOR.relative_to(ROOT)/'summary.json')]=sha(PRIOR/'summary.json')
    plan.update(status='FROZEN_COLD_FRESH_FIBRES',source_sha256=sources,
        prerequisite_summary_sha256=sha(PRIOR/'summary.json'),
        input_sha256={n:sha(OUT/n) for n in ('design.json','inputs.json','preflight.json')})
    new_write(OUT/'protocol.json',plan)
    print(json.dumps({'status':'SEALED_FRESH','sources':len(sources),'protocol_sha256':sha(OUT/'protocol.json')}),flush=True)


def sage_rank(packet):
    checker=SourceFileLoader('fresh_independent_rank',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module()
    return checker.sage_rank(packet)


def reconcile(case,dest):
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from pari_pointed_backend import replay
    from pointed_quartic_search import PointedQuarticSearch
    from v3_warm_engine import certified_state
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof']);admission=FinitePointAdmission(curve,basis,prime_bound=1000)
    result=load(dest/'result.json');gains=[]
    for record in result['calls']:
        path=dest/record['file'];require(sha(path)==record['sha256'],'cloud call hash mismatch')
        row=load(path);mapping=row['mapping']
        search=PointedQuarticSearch(state=state,centre={'coefficients':row['centre']},coordinate_policy=mapping['coordinate_policy'])
        for point in replay(search,mapping,row['search']):
            if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':
                gains.append({'call':record['file'],'point':list(map(str,point)),'after':len(admission.points)})
    proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
    packet={'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof}
    independent=sage_rank(packet)
    require(independent['rank']==len(basis)+len(gains),'full-cloud rank differs')
    return packet,{'status':'PASS_FULL_CLOUD_TWO_IMPLEMENTATIONS','rank':independent,
        'initial_rank':len(basis),'gains':gains,'point_search_calls_added':0,
        'boundary':'Certified lower bound from all returned points. Finite-column failure remains UNKNOWN.'}


def worker(case_id,arm):
    import cancellation_scheduler_cpu as cpu_worker
    from finite_cancellation_features import alarm
    from parent_foundry_worker import generic_seed,bank
    from future_point_admission import FinitePointAdmission
    plan=load(OUT/'protocol.json');request=next(c for c in load(OUT/'inputs.json') if c['id']==case_id)
    require(sha(OUT/'inputs.json')==plan['input_sha256']['inputs.json'],'recipe inputs changed')
    for name,h in plan['source_sha256'].items():require(sha(ROOT/name)==h,'sealed source/input changed: '+name)
    require(arm in plan['arms'],'unsealed arm')
    resource.setrlimit(resource.RLIMIT_CPU,(plan['hard_process_cpu_seconds']-5,plan['hard_process_cpu_seconds']))
    signal.signal(signal.SIGALRM,alarm)
    dest=RAW/'arms'/case_id/arm;cold=dest/'cold';parent=load(ROOT/request['parent'])
    before=cpu_worker.cpu();signal.alarm(plan['cold_wall_seconds']);failure=None;packet=None
    try:
        packet=generic_seed(parent,request['parameter'],cold)
        if packet is None:failure='UNRESOLVED_GENERIC_SPECIALIZATION'
        else:
            # The worker's fixed finite places must independently span this seed.
            FinitePointAdmission(tuple(map(F,packet['curve'])),[tuple(map(F,p)) for p in packet['points']],prime_bound=500)
            sage_rank(packet)
            centres=[]
            for i in range(3):
                folder=cold/f'bank-{i:02d}';bank(parent,packet,folder,i)
                centres.extend(r['word'] for r in load(folder/'anchor-bank.json')['rows'])
            require(len(centres)==48,'incomplete declared bank')
    except TimeoutError:failure='COLD_PREPARATION_TIMEOUT'
    except ArithmeticError as e:failure='COLD_PREPARATION_UNKNOWN: '+str(e)
    finally:signal.alarm(0)
    cold_receipt={'status':failure or 'PASS_COLD_PREPARATION','cpu_seconds':cpu_worker.cpu()-before,
        'parent_sha256':request['parent_sha256'],'global_protocol_sha256':sha(OUT/'protocol.json'),
        'files_sha256':{str(p.relative_to(dest)):sha(p) for p in sorted(cold.rglob('*')) if p.is_file()}}
    if failure:
        write(dest/'cold-receipt.json',cold_receipt)
        write(dest/'result.json',{'status':failure,'success':False,'initial_rank':None,
            'rank_lower_bound':None,'calls':[],'protocol_sha256':sha(OUT/'protocol.json'),
            'boundary':'Unresolved preparation; no negative point-search inference.'});return
    case={k:request[k] for k in ('id','family','stratum')}
    case.update(seed={k:packet[k] for k in ('curve','points','proof')},centres=centres)
    derived=dest/'derived';new_write(derived/'inputs.json',[case])
    localplan={**plan,'input_sha256':{'inputs.json':sha(derived/'inputs.json')},
        'global_protocol_sha256':sha(OUT/'protocol.json'),'parent_recipe_sha256':digest(canonical(request))}
    new_write(derived/'protocol.json',localplan)
    cold_receipt.update(derived_inputs_sha256=sha(derived/'inputs.json'),derived_protocol_sha256=sha(derived/'protocol.json'))
    write(dest/'cold-receipt.json',cold_receipt)
    cpu_worker.OUT=derived;cpu_worker.RAW=RAW;cpu_worker.run(case_id,arm)
    before=cpu_worker.cpu();cloud,verification=reconcile(case,dest)
    write(dest/'cloud-certificate.json',cloud)
    verification.update(certificate_sha256=sha(dest/'cloud-certificate.json'),cpu_seconds=cpu_worker.cpu()-before,
        result_sha256=sha(dest/'result.json'),derived_inputs_sha256=sha(derived/'inputs.json'))
    write(dest/'cloud-verification.json',verification)
    print(json.dumps({'status':'FRESH_COLD_ARM_COMPLETE','case':case_id,'arm':arm,
        'cloud_rank_lower_bound':verification['rank']['rank'],'cold_cpu_seconds':cold_receipt['cpu_seconds']}),flush=True)


def configure():
    round4.OUT=OUT;round4.RAW=RAW;round4.ENTRY=Path(__file__).resolve()


def run():configure();round4.run()
def pack():configure();round4.pack()


def report():
    from finite_cancellation_validation_audit import j_group
    from collections import Counter
    plan=load(OUT/'protocol.json');supervision=load(OUT/'supervision.json')
    require(supervision['status']=='COMPLETE' and len(supervision['records'])==8,'incomplete fresh pilot')
    # This is the first permitted read of target-j corpus data in this stage.
    corpus=load_corpus=json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))
    oldj={r['j_group'] for r in corpus};rows=[];paired={};status_counts=Counter()
    for receipt in supervision['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm'];result=load(dest/'result.json');cold=load(dest/'cold-receipt.json')
        require(receipt['status']=='COMPLETE' and receipt['result_sha256']==sha(dest/'result.json'),'invalid process receipt')
        row={'case':receipt['case'],'family':receipt['family'],'arm':receipt['arm'],'status':result['status'],
            'success':result['success'],'calls':len(result['calls']),'cpu_seconds':receipt['charged_cpu_seconds'],
            'cold_cpu_seconds':cold['cpu_seconds'],'rank_lower_bound':result['rank_lower_bound']}
        if cold['status']=='PASS_COLD_PREPARATION':
            case=load(dest/'derived/inputs.json')[0];verification=load(dest/'independent-verification.json')
            cloud=load(dest/'cloud-verification.json');j=j_group(case['seed']['curve'])
            require(verification['status']=='PASS' and result['independent_verification_sha256']==sha(dest/'independent-verification.json'),'invalid independent proof')
            require(cloud['status']=='PASS_FULL_CLOUD_TWO_IMPLEMENTATIONS' and cloud['certificate_sha256']==sha(dest/'cloud-certificate.json') and cloud['result_sha256']==sha(dest/'result.json'),'invalid cloud proof')
            require(cold['derived_inputs_sha256']==sha(dest/'derived/inputs.json') and cold['derived_protocol_sha256']==sha(dest/'derived/protocol.json'),'cold input seal differs')
            require(result['protocol_sha256']==sha(dest/'derived/protocol.json'),'worker used other input seal')
            for name,h in cold['files_sha256'].items():require(sha(dest/name)==h,'cold evidence changed')
            old=paired.setdefault(receipt['case'],cold['derived_inputs_sha256']);require(old==cold['derived_inputs_sha256'],'cold arm inputs differ')
            for record in result['calls']:status_counts[record['status']]+=1
            row.update(initial_rank=len(case['seed']['points']),cloud_rank_lower_bound=cloud['rank']['rank'],
                j_group=j,absent_retained_corpus=j not in oldj,preparation_unknowns=verification['preparation_unknowns'],
                cloud_certificate=str((dest/'cloud-certificate.json').relative_to(ROOT)))
        rows.append(row)
    totals={a:{'gains':sum(r['success'] for r in rows if r['arm']==a),
        'cpu_seconds':sum(r['cpu_seconds'] for r in rows if r['arm']==a),
        'cold_cpu_seconds':sum(r['cold_cpu_seconds'] for r in rows if r['arm']==a),
        'calls':sum(r['calls'] for r in rows if r['arm']==a)} for a in plan['arms']}
    ratio=totals[plan['candidate']]['cpu_seconds']/totals['factor_free']['cpu_seconds']
    write(OUT/'summary.json',{'status':'COMPLETE_COLD_FRESH_PILOT','rows':rows,'totals':totals,
        'candidate_cpu_ratio':ratio,'point_call_statuses':dict(status_counts),
        'protocol_sha256':sha(OUT/'protocol.json'),'novelty_corpus_sha256':sha(CORPUS/'corpus.json.gz'),
        'rank32':'UNKNOWN','boundary':'Four predetermined generic-only fibres; descriptive paired cold costs, no population speed claim. Freshness means absent retained corpus and prior parent/address exposure. All ranks are certified lower bounds. First-direction and full-cloud outcomes are separate. No additional point calls or enlarged-basis search.'})
    print(json.dumps({'rows':rows,'totals':totals,'candidate_cpu_ratio':ratio,'point_call_statuses':dict(status_counts)},indent=2),flush=True)


def replay():
    from parent_foundry_worker import verify,bank
    from verify_cancellation_scheduler import verify_arm
    plan=load(OUT/'protocol.json');requests={r['id']:r for r in load(OUT/'inputs.json')};rows=[];start=time.process_time()
    for receipt in load(OUT/'supervision.json')['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm'];cold=load(dest/'cold-receipt.json')
        if cold['status']!='PASS_COLD_PREPARATION':
            rows.append({'case':receipt['case'],'arm':receipt['arm'],'status':'RETAINED_UNKNOWN'});continue
        request=requests[receipt['case']];parent=load(ROOT/request['parent']);packet=load(dest/'cold/generic.json')
        require(sha(ROOT/request['parent'])==request['parent_sha256'],'parent hash changed')
        verify(packet,parent,request['parameter']);sage_rank(packet)
        case=load(dest/'derived/inputs.json')[0];centres=[]
        with tempfile.TemporaryDirectory(prefix='cancellation-fresh-replay-') as temp:
            for i in range(3):
                folder=Path(temp)/f'bank-{i:02d}';bank(parent,packet,folder,i)
                for p in folder.iterdir():require(p.read_bytes()==(dest/'cold'/folder.name/p.name).read_bytes(),'generic bank replay differs')
                centres.extend(r['word'] for r in load(folder/'anchor-bank.json')['rows'])
        require(case['centres']==centres and case['seed']=={k:packet[k] for k in ('curve','points','proof')},'cold recipe reconstruction differs')
        checked=verify_arm(case,receipt['arm'],dest,load(dest/'derived/protocol.json'),load(dest/'result.json'))
        cloud,v=reconcile(case,dest);require(cloud==load(dest/'cloud-certificate.json'),'cloud certificate replay differs')
        rows.append({'case':receipt['case'],'arm':receipt['arm'],'status':'PASS','first_direction':checked['rank'],'cloud':v['rank']})
    write(OUT/'replay.json',{'status':'PASS_WITH_EXPLICIT_UNKNOWNS','rows':rows,'cpu_seconds':time.process_time()-start,
        'protocol_sha256':sha(OUT/'protocol.json')})
    print(json.dumps({'status':'PASS_WITH_EXPLICIT_UNKNOWNS','arms':len(rows)}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','seal','run','worker','report','pack','replay'])
    parser.add_argument('--case');parser.add_argument('--arm');args=parser.parse_args()
    worker(args.case,args.arm) if args.command=='worker' else globals()[args.command]()
