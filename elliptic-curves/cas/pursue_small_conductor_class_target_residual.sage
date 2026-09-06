#!/usr/bin/env sage-python
"""Search inexpensive prime representatives of unresolved quotient directions, with exact replay."""
import argparse
import hashlib
import json
from math import gcd, prod
from pathlib import Path
import runpy
import time
from sage.all import AA, QQ, PolynomialRing, pari
import extend_small_conductor_norm_batch as batch
import capped_primorial_remainders as capped
batch.residues = capped.residues
from audit_bnf_free_s_class_quotient import packed_rank
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT, ART, cert = batch.ROOT, batch.ART, batch.cert
SOURCE = Path(__file__).resolve()
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-class-target-residual-v1'
BASE = ART/'small_conductor_norm_batch_relations_v1.json'
EXTRA = ART/'small_conductor_special_prime_relations_v1.json'
SMALL = ART/'small_conductor_smaller_base_v2.json'
PREVIOUS = ART/'small_conductor_small_base_relations_v1.json'
OLD_SOURCE = ROOT/'elliptic-curves/cas/target_small_conductor_small_base.sage'
old = runpy.run_path(str(OLD_SOURCE))
SEED_SOURCE = ROOT/'elliptic-curves/cas/pursue_small_conductor_class_target_capped.sage'
seed = runpy.run_path(str(SEED_SOURCE))
SEED_WAVES = 2
SLOPE_SCALE = 2**96
CHAR_SOURCE = ROOT/'elliptic-curves/cas/certify_small_conductor_class_characters.sage'
CHARACTERS = ART/'small_conductor_class_characters_v1.json'
ANCHORS = {r['column'] for r in cert.read(CHARACTERS)['anchors']}
ANCHOR_MASK = sum(1<<c for c in ANCHORS)
BOUND = 37638


def directory(wave):
    return D/('wave_%03d'%wave)


def audit_path(wave):
    return ART/('small_conductor_class_target_residual_wave_%03d_v1.json'%wave)


def bits(row):
    return sum(1<<c for c in row)


def add_basis(basis, row):
    while row:
        unprotected = row & ~ANCHOR_MASK
        if not unprotected: raise ArithmeticError('relation contradicts independent ideal anchors')
        p = unprotected.bit_length()-1
        if p not in basis:
            basis[p] = row
            return True
        row ^= basis[p]
    return False


class Matrix:
    """Eliminate outside coordinates exactly, then supported coordinates over F2."""
    def __init__(self, base):
        self.base = base
        self.canonical = {max(r['parity_columns']):bits(r['parity_columns'])
                          for r in base['canonical_rational_relations']}
        self.inside_mask = bits(i for i,c in enumerate(base['columns']) if c['p']<=BOUND)
        self.free = [i for i,c in enumerate(base['columns']) if c['p']<=BOUND and i not in self.canonical]
        self.inside, self.outside, self.rows = {}, {}, []

    def add(self, factors):
        row = bits(c for c,e in factors if e%2)
        for c,e in factors:
            if e%2 and c in self.canonical:
                row ^= self.canonical[c]
        self.rows.append(row)
        inside, outside = row & self.inside_mask, row & ~self.inside_mask
        while outside:
            p = outside.bit_length()-1
            if p not in self.outside:
                self.outside[p] = (outside, inside)
                return False
            o,i = self.outside[p]
            outside ^= o
            inside ^= i
        return add_basis(self.inside, inside)

    def dimension(self):
        return len(self.free)-len(self.inside)

    def report(self):
        # Independent rank identity checks the supported-intersection dimension.
        total = packed_rank(self.rows)
        outside = packed_rank([r & ~self.inside_mask for r in self.rows])
        if total-outside != len(self.inside) or outside != len(self.outside):
            raise ArithmeticError('intersection rank identity failed')
        return {'free_columns':len(self.free), 'all_relation_rank':total,
                'outside_projection_rank':outside, 'supported_rank':len(self.inside),
                'quotient_dimension':self.dimension(), 'further_rows_to_16':max(0,self.dimension()-16)}


def baseline():
    matrix = seed['prior_state'](SEED_WAVES+1)
    result = Matrix(matrix.base)
    # Transfer the seed's already normalized relation rows. Canonical rational
    # relations are principal, so normalization does not change validity mod2.
    for row in matrix.rows:
        factors = []
        while row:
            bit = row & -row; factors.append([bit.bit_length()-1,1]); row ^= bit
        result.add(factors)
    if result.report()!=matrix.report():
        raise ArithmeticError('seed matrix differs')
    # The lower-bound proof also gives exact principal parity relations for
    # each known point: every omitted valuation is proved even.
    lower=cert.read(ART/'small_conductor_class_lower16_v1.json')
    lookup={(str(c['p']),c['hnf']):i for i,c in enumerate(result.base['columns'])}
    for j in range(22):
        factors=[[lookup[(r['p'],r['hnf'])],1] for r in lower['bad_prime_valuations'] if r['valuations'][j]%2]
        result.add(factors)
    if result.dimension()<16:raise ArithmeticError('point parity rows contradict proved lower bound')
    return result


def protocol(wave):
    p = cert.read(directory(wave)/'protocol.json')
    for name,h in p['sources'].items():
        if cert.hashed(ROOT/name)!=h:
            raise ArithmeticError('frozen source differs: '+name)
    return p


def prior_state(wave):
    matrix = baseline()
    for k in range(1,wave):
        a = cert.read(audit_path(k))
        if a['status']!='PASS':
            raise ArithmeticError('previous wave not audited')
        p = protocol(k)
        for item in a['chunks']:
            path = ROOT/item['path']
            if cert.hashed(path)!=item['sha256']:
                raise ArithmeticError('prior witness chunk differs')
            for r in cert.read(path)['relations']:
                matrix.add(r['ideal_factorization'])
        if matrix.dimension()!=a['matrix']['quotient_dimension']:
            raise ArithmeticError('prior matrix replay differs')
    return matrix


def history(wave, smooth_bound):
    seen = seed['history'](SEED_WAVES+1,smooth_bound)
    for k in range(1,wave):
        p = protocol(k)
        if p['smooth_bound']<smooth_bound: continue
        # Only non-skipped completed target chunks count as searched.
        for item in cert.read(audit_path(k))['chunks']:
            chunk = cert.read(ROOT/item['path'])
            if chunk.get('skipped'): continue
            t = chunk['target']
            seen[t['column']] = max(seen.get(t['column'],0),p['vmax'])
    return seen


def residual_projection(matrix,column):
    # Normal form in the presented quotient modulo the16 anchor directions.
    # This is only a target-selection computation; no relation is inferred.
    row=(1<<column)&matrix.inside_mask
    outside=(1<<column)&~matrix.inside_mask
    while outside:
        pivot=outside.bit_length()-1
        if pivot not in matrix.outside:return 0
        o,i=matrix.outside[pivot];outside^=o;row^=i
    row&=~ANCHOR_MASK;result=0
    while row:
        pivot=row.bit_length()-1
        if pivot in matrix.inside:row^=matrix.inside[pivot]
        else:result^=1<<pivot;row^=1<<pivot
        row&=~ANCHOR_MASK
    return result


def selection(matrix, wave, box, count, smooth_bound):
    data,poly,nf,x,w = old['setup']()
    a = int(data['integral_norm_generator']['fixed_a'])
    M = data['sl2_matrix']
    seen = history(wave,smooth_bound)
    ring = PolynomialRing(QQ,'t')
    roots = ring(list(reversed(list(map(int,data['reduced_binary_cubic_descending']))))).roots(AA,multiplicities=False)
    if len(roots)!=3: raise ArithmeticError('three real norm roots required')
    eligible = []
    for i,c in enumerate(matrix.base['columns']):
        q = c['p']
        if i in matrix.canonical or q>smooth_bound:continue
        if i in ANCHORS or c['f']!=1 or c['e']!=1 or seen.get(i,0)>=box:
            continue
        if int(nf.disc())%q==0 or int(data['defining_order_index'])%q==0 or a%q==0:
            continue
        if residual_projection(matrix,i):eligible.append(i)
    span={};selected=[]
    for i in eligible:
        if add_basis(span,residual_projection(matrix,i)):selected.append(i)
        if len(selected)==count:break
    chosen=set(selected)
    selected+= [i for i in eligible if i not in chosen][:count-len(selected)]
    selected.sort()
    targets = []
    for column in selected:
        col = matrix.base['columns'][column]; q = col['p']
        root = next(int(pari.lift(r)) for r in pari.polrootsmod(poly,q)
                    if str(pari.idealhnf(nf,q,x-int(pari.lift(r))))==col['hnf'])
        wbar = 0
        for i,c in enumerate(data['integral_norm_generator']['w_power_basis']):
            c = QQ(c)
            wbar += int(c.numerator())*pow(int(c.denominator()),-1,q)*pow(root,i,q)
        bm,bn = (a*M[0][0]+wbar*M[1][0])%q,(a*M[0][1]+wbar*M[1][1])%q
        slope = -bn*pow(bm,-1,q)%q
        lattice = old['reduce_lattice'](q,slope,data['reduced_hessian'])
        v1,v2 = lattice
        slopes = [int((((r*v2[1]-v2[0])/(v1[0]-r*v1[1]))*SLOPE_SCALE).floor()) for r in roots]
        targets.append({'column':column,'p':q,'hnf':col['hnf'],'theta_root':root,
                        'norm_slope':slope,'lattice_basis':lattice,'root_slopes_scaled':slopes,
                        'unresolved_projection_hex':hex(residual_projection(matrix,column)),
                        'old_vmax':seen.get(column,0)})
    return {'eligible_columns':len(eligible),'selected_projection_rank':packed_rank([residual_projection(matrix,i) for i in selected]),'targets':targets}


def prepare(args):
    path = directory(args.wave)/'protocol.json'
    if path.exists():
        raise FileExistsError('preserve protocol')
    if not 2<=args.targets<=4096 or not 1<=args.box<=1048576 or not BOUND<=args.smooth_bound<=400000:
        raise ValueError('parameters outside declared bounds')
    matrix = prior_state(args.wave)
    choose = selection(matrix,args.wave,args.box,args.targets,args.smooth_bound)
    if not choose['targets']:
        raise ArithmeticError('no eligible targets at this box; choose a larger bounded box')
    paths = [SOURCE,Path(capped.__file__).resolve(),ART/'small_conductor_class_lower16_v1.json',SEED_SOURCE,CHAR_SOURCE,CHARACTERS,OLD_SOURCE,BASE,EXTRA,SMALL,PREVIOUS,Path(batch.__file__).resolve(),batch.forms.OUT,
             ROOT/'elliptic-curves/cas/audit_bnf_free_s_class_quotient.py']
    paths += [seed['audit_path'](k) for k in range(1,SEED_WAVES+1)]
    paths += [audit_path(k) for k in range(1,args.wave)]
    checkpoint(path,{'schema':'elliptic-curves.small-conductor-class-target-residual-protocol.v1',
        'wave':args.wave,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'selection':choose,'starting_dimension':matrix.dimension(),'vmax':args.box,
        'smooth_bound':args.smooth_bound,'workers':1,'rss_bytes':1610612736,
        'worker_seconds':args.seconds,'audit_seconds':600,
        'point_parity_relations_seeded':True,'remainder_method':'Exact product tree capped at primorial+1; every leaf remainder is unchanged.',
        'protected_anchor_columns':sorted(ANCHORS),'seed_waves':SEED_WAVES,'slope_scale':str(SLOPE_SCALE),
        'gate':'User goal: reduce the GRH-certified small-base quotient to16. Protect the16 independently certified ideal classes and skip a target only when its current normal form lies in the certified anchor span. Select prime representatives with nonzero unresolved normal form, including pivot columns and outside-base ideals whose outside normal form is zero; first choose a greedy spanning set in ascending prime order and then fill with the smallest eligible primes. Use thin strips near real zeros of the transformed norm form, to test lower norms per enumerated candidate. Measure exact independent supported rank gain.',
        'population':'For each exact index-q lattice and each of three real norm roots, center u at floor(v*floor(root_in_lattice*2^96)/2^96). For each new integer v<=vmax, try center-1, center, center+1; deduplicate pairs and require primitive(u,v) and mapped(m,n). Skip a target only if its normal form has no unresolved coordinate at its turn. Stop after a completed target if dimension<=16.',
        'verification':'Every retained principal ideal and exact norm identity is independently replayed. No exhaustive miss or smoothness-completeness claim; rejected candidates need not enter the rank certificate.',
        'outside_rule':'If the smoothness bound exceeds37638, keep all extra coordinates through exact elimination. Admit only combinations with zero outside component; verify rank(all)-rank(outside).',
        'claim_boundary':'Generation and resulting class/rank upper bounds are GRH-conditional. Known curve rank is at least22. A residual<=16 would certify rank22 under GRH, not unconditionally.'})
    print('PREPARED',args.wave,'TARGETS',len(choose['targets']),'DIMENSION',matrix.dimension(),'VMAX',args.box,flush=True)


def context():
    data,poly,nf,x,w = old['setup']()
    base = cert.read(BASE)
    lookup = {(r['p'],r['hnf']):i for i,r in enumerate(base['columns'])}
    return data,nf,w,lookup


def exact_relation(m,n,value,target,ctx,bound):
    data,nf,w,lookup = ctx
    a = int(data['integral_norm_generator']['fixed_a']); M = data['sl2_matrix']
    beta = a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
    if pari.nfeltnorm(nf,beta)!=a*a*value:
        raise ArithmeticError('principal norm identity failed')
    factor = pari.idealfactor(nf,beta)
    if pari.idealhnf(nf,pari.idealfactorback(nf,factor))!=pari.idealhnf(nf,beta):
        raise ArithmeticError('principal ideal identity failed')
    row = []
    for j in range(factor.nrows()):
        ideal,e = factor[j,0],int(factor[j,1]); q = int(ideal[0])
        if e<0 or q>bound:
            raise ArithmeticError('factor outside declared support')
        row.append([lookup[(q,str(pari.idealhnf(nf,ideal)))],e])
    if not any(c==target['column'] and e>0 for c,e in row):
        raise ArithmeticError('intended prime ideal absent')
    return {'m':m,'n':n,'value':str(value),
            'beta_power_basis':[str(pari.polcoef(pari.lift(beta),i)) for i in range(3)],
            'ideal_factorization':row}


def worker(wave):
    p = protocol(wave); matrix = prior_state(wave); ctx = context()
    primes = batch.prior.primes_to(p['smooth_bound']); primorial = prod(primes)
    c0,c1,c2,c3 = map(int,ctx[0]['reduced_binary_cubic_descending'])
    b = p['vmax']; chunks = []; started = time.monotonic(); count = 0
    for index,target in enumerate(p['selection']['targets']):
        path = directory(wave)/('target_%04d.json'%index)
        if path.exists():
            chunk = cert.read(path)
            if chunk['target']!=target:
                raise ArithmeticError('resume target differs')
            for r in chunk['relations']:
                matrix.add(r['ideal_factorization'])
        elif not residual_projection(matrix,target['column']):
            chunk = {'target':target,'candidate_occurrences':0,'relations':[],
                     'skipped':'normal form supported only on certified anchors','ending_dimension':matrix.dimension(),
                     'wall_seconds':0.0,'population_digest_sha256':hashlib.sha256(b'').hexdigest()}
            checkpoint(path,chunk)
        else:
            begin = time.monotonic(); v1,v2 = target['lattice_basis']; pairs = []; values = []
            digest = hashlib.sha256()
            seen_pairs = set()
            for v in range(target['old_vmax']+1,b+1):
                for slope in target['root_slopes_scaled']:
                    center = slope*v//SLOPE_SCALE
                    for u in range(center-1,center+2):
                        if (u,v) in seen_pairs or gcd(u,v)!=1:
                            continue
                        seen_pairs.add((u,v))
                        m,n = u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                        if gcd(m,n)!=1:
                            continue
                        value = ((c0*m+c1*n)*m+c2*n*n)*m+c3*n*n*n
                        if not value or value%target['p'] or (m-target['norm_slope']*n)%target['p']:
                            raise ArithmeticError('norm lattice identity failed')
                        pairs.append((u,v,m,n)); values.append(value)
            relations = []
            if values:
                for pair,value,r in zip(pairs,values,batch.residues(list(map(abs,values)),primorial)):
                    remainder = batch.strip(abs(value),gcd(abs(value),r))
                    digest.update(('%s,%s,%s,%s\n'%(pair[0],pair[1],value,remainder)).encode())
                    if remainder==1:
                        u,v,m,n = pair
                        relation = exact_relation(m,n,value,target,ctx,p['smooth_bound'])
                        relation.update({'u':u,'v':v})
                        relations.append(relation); matrix.add(relation['ideal_factorization'])
            chunk = {'target':target,'candidate_occurrences':len(values),'population_digest_sha256':digest.hexdigest(),
                     'relations':relations,'ending_dimension':matrix.dimension(),'wall_seconds':time.monotonic()-begin}
            checkpoint(path,chunk)
        count += chunk['candidate_occurrences']
        chunks.append({'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)})
        checkpoint(directory(wave)/'progress.json',{'completed_targets':index+1,'candidates':count,
                   'dimension':matrix.dimension(),'rank_gain':p['starting_dimension']-matrix.dimension(),
                   'wall_seconds':time.monotonic()-started})
        if (index+1)%16==0 or matrix.dimension()<=16:
            print('WAVE',wave,'TARGETS',index+1,'DIMENSION',matrix.dimension(),'GAIN',p['starting_dimension']-matrix.dimension(),flush=True)
        if matrix.dimension()<=16:
            break
    path = directory(wave)/'result.json'
    if path.exists():
        raise FileExistsError('preserve wave result')
    checkpoint(path,{'status':'PASS','chunks':chunks,'candidate_occurrences':count,'matrix':matrix.report(),
                     'protocol_sha256':cert.hashed(directory(wave)/'protocol.json'),
                     'claim_boundary':p['claim_boundary']})
    print('WORKER PASS',wave,'DIMENSION',matrix.dimension(),flush=True)


def audit_all(last,write=False):
    # Replay every inherited proof once, then every new witness and rank transition.
    characters=runpy.run_path(str(CHAR_SOURCE))['expected']()
    if characters!=cert.read(CHARACTERS): raise ArithmeticError('anchor character proof differs')
    seed['audit_all'](SEED_WAVES,write=False)
    matrix = baseline(); ctx = context()
    c = list(map(int,ctx[0]['reduced_binary_cubic_descending']))
    for wave in range(1,last+1):
        p = protocol(wave); summary = cert.read(directory(wave)/'result.json')
        if summary['status']!='PASS' or summary['protocol_sha256']!=cert.hashed(directory(wave)/'protocol.json'):
            raise ArithmeticError('wave summary not bound to protocol')
        if matrix.dimension()!=p['starting_dimension'] or selection(matrix,wave,p['vmax'],len(p['selection']['targets']),p['smooth_bound'])!=p['selection']:
            raise ArithmeticError('target-selection replay failed')
        start = matrix.dimension(); count = 0; relations = 0; distinct = set()
        for index,item in enumerate(summary['chunks']):
            path = directory(wave)/('target_%04d.json'%index)
            if str(path.relative_to(ROOT))!=item['path'] or cert.hashed(path)!=item['sha256']:
                raise ArithmeticError('chunk binding failed')
            chunk = cert.read(path); target = p['selection']['targets'][index]
            if chunk['target']!=target:
                raise ArithmeticError('target differs')
            should_skip = not residual_projection(matrix,target['column'])
            if bool(chunk.get('skipped'))!=should_skip:
                raise ArithmeticError('adaptive skip differs')
            if should_skip and (chunk['relations'] or chunk['candidate_occurrences']):
                raise ArithmeticError('skipped target has candidates')
            for r in chunk['relations']:
                u,v,m,n = r['u'],r['v'],r['m'],r['n']; v1,v2 = target['lattice_basis']
                if not (target['old_vmax']<v<=p['vmax'] and any(abs(u-slope*v//SLOPE_SCALE)<=1 for slope in target['root_slopes_scaled']) and gcd(u,v)==gcd(m,n)==1):
                    raise ArithmeticError('witness outside primitive population')
                if (m,n)!=(u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]):
                    raise ArithmeticError('witness lattice map failed')
                value = sum(a*m**(3-i)*n**i for i,a in enumerate(c))
                expected = exact_relation(m,n,value,target,ctx,p['smooth_bound'])
                expected.update({'u':u,'v':v})
                if expected!=r:
                    raise ArithmeticError('principal relation replay failed')
                matrix.add(r['ideal_factorization']); relations += 1
                distinct.add((m,n) if n>0 or n==0 and m>0 else (-m,-n))
            if matrix.dimension()!=chunk['ending_dimension']:
                raise ArithmeticError('target checkpoint matrix differs')
            count += chunk['candidate_occurrences']
        if len(summary['chunks'])!=len(p['selection']['targets']) and matrix.dimension()>16:
            raise ArithmeticError('incomplete wave without stopping target')
        result = matrix.report()
        if result!=summary['matrix'] or count!=summary['candidate_occurrences']:
            raise ArithmeticError('summary matrix differs')
        record = {'schema':'elliptic-curves.small-conductor-class-target-residual-wave.v1','status':'PASS','wave':wave,
                  'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in [SOURCE,directory(wave)/'protocol.json',directory(wave)/'result.json']},
                  'chunks':summary['chunks'],'candidate_occurrences':count,'relation_occurrences':relations,
                  'distinct_pairs_up_to_sign':len(distinct),'starting_dimension':start,'independent_supported_gain':start-matrix.dimension(),
                  'matrix':result,'conditional_on_grh_class_two_rank_upper_bound':matrix.dimension(),
                  'conditional_on_grh_curve_rank_upper_bound':2*((matrix.dimension()+7)//2),
                  'point_parity_relations_seeded':True,'protected_anchor_columns':sorted(ANCHORS),'unconditional_class_two_rank_lower_bound':16,
                  'unconditional_rank_lower_bound':22,'unconditional_rank_upper_bound':None,
                  'conditional_on_grh_exact_rank':22 if matrix.dimension()<=16 else None,
                  'claim_boundary':p['claim_boundary']}
        if matrix.dimension()<16 or record['conditional_on_grh_curve_rank_upper_bound']<22:
            raise ArithmeticError('contradiction with known lower bound; fail closed')
        if write and wave==last:
            if audit_path(wave).exists():
                raise FileExistsError('preserve audit')
            checkpoint(audit_path(wave),record)
        elif cert.read(audit_path(wave))!=record:
            raise ArithmeticError('audit record differs')
        print('AUDIT PASS',wave,'GAIN',start-matrix.dimension(),'DIMENSION',matrix.dimension(),flush=True)


def launch(args):
    p = protocol(args.wave); stage = args.stage[7:]
    previous = list(directory(args.wave).glob(stage+'_*.supervisor.json'))
    name = stage+'_%02d'%len(previous)
    command = ['/home/royvanrijn/.local/bin/sage','-python',str(SOURCE),stage,'--wave',str(args.wave)]
    outcome = run(command,limits=Limits(p['worker_seconds'] if stage=='worker' else p['audit_seconds'],p['rss_bytes']),
                  cwd=ROOT,log_path=directory(args.wave)/(name+'.log'),checkpoint_path=directory(args.wave)/(name+'.supervisor.json'))
    print(stage,args.wave,outcome['outcome'],outcome['returncode'],flush=True)
    if outcome['outcome']!='completed' or outcome['returncode']!=0:
        raise SystemExit(1)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','worker','audit','check','launch-worker','launch-audit'])
    parser.add_argument('--wave',type=int,required=True)
    parser.add_argument('--vmax',dest='box',type=int,default=1024)
    parser.add_argument('--targets',type=int,default=512)
    parser.add_argument('--smooth-bound',type=int,default=50000)
    parser.add_argument('--seconds',type=int,default=600)
    args = parser.parse_args()
    if args.stage=='prepare': prepare(args)
    elif args.stage=='worker': worker(args.wave)
    elif args.stage in ['audit','check']: audit_all(args.wave,write=args.stage=='audit')
    else: launch(args)
