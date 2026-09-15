#!/usr/bin/env python3
"""Independent factor-incidence/contact replay and complete norm-census validation."""
import argparse
from collections import defaultdict
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path
import resource
import subprocess
import tempfile
import time
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
K1='artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1'
PRODUCER='elkies-k3/scripts/certify_q80_genus_one_k0_gate.sage'
CPP='elkies-k3/scripts/certify_q80_cubic_rational_norms.cpp'
REPLAY='elkies-k3/scripts/replay_q80_cubic_rational_norms.cpp'
HELPER='elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
QUARTIC='elkies-k3/scripts/run_q80_quartic_norm_census.py'
P=131
def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
h=load(ROOT/HELPER,'q80_k0_polynomials')
require=h.require;read=h.read;digest=h.digest
add=h.add;sub=h.sub;mul=h.mul;scale=h.scale;divrem=h.divrem;gcd=h.gcd;exact=h.exact;ev=h.ev;rank=h.rank;trim=h.trim
def low_factors(y):
    """All degree <=2 factors, independent of Sage, with exact multiplicities."""
    y=trim(y);original=y;out=[]
    for b in range(P):
        f=[-b%P,1];n=0
        while y and not ev(y,b):y=exact(y,f);n+=1
        if n:out.append((f,n))
    if len(y)<3:return out
    g=gcd(y,sub(h.powmod([0,1],P*P,y),[0,1]));quadratics=[]
    while len(g)>3:
        factor=None
        for b in range(P):
            for a in range(P):
                if pow((b*b-4*a)%P,65,P)!=P-1:continue
                f=[a,b,1]
                if not divrem(g,f)[1]:factor=f;break
            if factor is not None:break
        require(factor is not None,'split every quadratic-factor product');quadratics.append(factor);g=exact(g,factor)
    if len(g)==3:quadratics.append(g)
    else:require(len(g)==1,'degree-two distinct-factor remainder')
    for f in quadratics:
        yy=original;n=0
        while not divrem(yy,f)[1]:yy=exact(yy,f);n+=1
        out.append((f,n))
    return out
def quadratic_divisors(y):
    fac=low_factors(y);hh=[];linears=[f for f,n in fac if len(f)==2]
    for f,n in fac:
        if len(f)==3:hh.append(f)
        elif n>=2:hh.append(mul(f,f))
    hh += [mul(a,b) for a,b in itertools.combinations(linears,2)]
    infinity=7-len(trim(y))
    if infinity:hh+=linears
    if infinity>=2:hh.append([1])
    return [tuple(f+[0]*(3-len(f))) for f in hh]
def check_contacts(record,points):
    groups=defaultdict(list);nodal=[]
    for i,point in enumerate(points):
        require(point['index']==i,'literal norm4 section indices')
        for H in quadratic_divisors(point['y']):groups[H].append(i)
        yy=point['y'];n=0
        while not divrem(yy,[88,62,1])[1]:yy=exact(yy,[88,62,1]);n+=1
        if n:nodal.append({'index':i,'q0_order':n,'x_mod_q0':divrem(point['x'],[88,62,1])[1]})
    shared=set();found=defaultdict(list)
    for H,indices in groups.items():
        for i,j in itertools.combinations(indices,2):
            shared.add((i,j));diff=sub(points[i]['x'],points[j]['x'])
            if H[2]==0 and ((points[i]['x']+[0]*5)[4]==(points[j]['x']+[0]*5)[4]):continue
            if len(gcd(diff,trim(H)))>1:continue
            found[(i,j)].append(list(H))
    rows=[{'pair':list(pair),'H':sorted(found[pair])} for pair in sorted(found,key=lambda a:(a[1],a[0]))]
    require(record['pairs_examined']==len(points)*(len(points)-1)//2==861328,'complete unordered section range')
    require(record['shared_degree_at_least2']==len(shared)==69,'complete shared quadratic incidence')
    require(record['distinct_root_contact_pairs']==rows and len(rows)==30,'complete independent distinct-root pair roster')
    require(record['sections_with_q0_ordinate_factor']==nodal==[{'index':732,'q0_order':1,'x_mod_q0':[73,6]}],'unique ordinary norm4 nodal contact')
    return rows
def first_lift(pair,H,points,model):
    D=mul(H,H);A=h.fraction_poly(model['A_coefficients_low_to_high'],P)
    A2=h.fraction_poly(model['A_coefficients_low_to_high'],P*P);B2=h.fraction_poly(model['B_coefficients_low_to_high'],P*P)
    require(len(H)==3 and H[2]==1,'monic quadratic contact chart')
    J=[[0]*24 for _ in range(26)];target=[]
    for n,index in enumerate(pair):
        x=points[index]['x'];r=exact(points[index]['y'],H)
        require(len(x)<=5 and len(r)<=5,'quartic coordinate budgets')
        # Put common D columns first; this is a separate Gaussian elimination.
        cols=[[0]*k+scale(mul(r,r),-1) for k in range(4)]
        cols += [[0]*k+scale(mul(D,r),-2) for k in range(5)]
        cols += [[0]*k+add(scale(mul(x,x),3),A) for k in range(5)]
        for j,f in enumerate(cols):
            dest=j if j<4 else 4+10*n+j-4
            for k,v in enumerate(f):J[13*n+k][dest]=v
        error=sub(add(add(mul(mul(x,x,P*P),x,P*P),mul(A2,x,P*P),P*P),B2,P*P),mul(D,mul(r,r,P*P),P*P),P*P)
        error += [0]*(13-len(error));require(all(c%P==0 for c in error),'literal first-lift divisibility')
        target += [(-c//P)%P for c in error]
    return rank(J),rank([row+[c] for row,c in zip(J,target)]),target
def check_composite(record,rat,quad):
    rs=[r for r in rat if len(r['rational_roots'])==3];qs=[r for r in quad if len(r['roots'])==3]
    codes=[tuple(r['rational_codes']) for r in rs];rp=defaultdict(set)
    # Fix the first fibre's labelling, quotienting the simultaneous S3 action.
    for i,j in itertools.combinations(range(len(rs)),2):
        for b in set(itertools.permutations(codes[j])):
            rp[tuple(sorted(a^c for a,c in zip(codes[i],b)))].add((i,j))
    four=sum(set(a).isdisjoint(b) for rows in rp.values() for a,b in itertools.combinations(rows,2))
    mixed=sum(bool(rp.get(tuple(sorted(r['norm_codes'])))) for r in qs)
    qkeys=[tuple(sorted(r['norm_codes'])) for r in qs]
    require(len(qkeys)==len(set(qkeys)),'no distinct quadratic signature collision')
    require(four==mixed==0,'no four-rational or mixed squarefree norm tuple')
    expected={'rational_pair_choices':4320,'four_rational_matches':0,'two_rational_quadratic_matches':0,'two_quadratic_matches':0,'rational_signatures':[sorted(c) for c in codes]}
    require(record==expected and len(codes)==16,'complete composite gate and signature input')
def check_quartic(out,recompute=False):
    q=load(ROOT/QUARTIC,'q80_k0_quartic');packet,text=q.validate_packet(out)
    require(packet['pair_ranges']==q.RANGES and len(q.RANGES)==10,'ten exact quartic checkpoint ranges')
    aggregates=[]
    for kind,source in [('producer',q.PRODUCER),('replay',q.REPLAY)]:
        rows=[];record=read(out/('quartic-'+kind+'.json'))
        for lo,hi in q.RANGES:
            path=out/('quartic-'+kind)/f'{lo:04d}-{hi:04d}.json';rr=read(path)
            require(rr['source_sha256']==digest(ROOT/source) and rr['input_sha256']==digest(out/'quartic-input.json'),'quartic source/input binding')
            require(record['checkpoints'][str(path.relative_to(out))]==digest(path),'quartic checkpoint hash')
            r=rr['result'];require(r['status']=='COMPLETE_RANGE' and [r['pair_begin'],r['pair_end']]==[lo,hi] and r['orbits']==(hi-lo)*P*P,'gap-free genuine-degree-four range')
            require(len(r['first_failed_character'])==18 and sum(r['first_failed_character'])==r['smooth_split'],'full split partition')
            require(r['first_failed_character'][-1]==len(r['candidate_rows'])==0,'no all-zero character tuple');rows.append(r)
        require(record['orbits']==sum(r['orbits'] for r in rows)==73620690,'independent quartic coverage total')
        require(record['smooth_split']==sum(r['smooth_split'] for r in rows)==12267623,'complete split-fibre total')
        require(record['first_failed_character']==[sum(r['first_failed_character'][i] for r in rows) for i in range(18)],'aggregate character partition')
        aggregates.append(rows)
    for a,b in zip(*aggregates):
        for k in ['status','pair_begin','pair_end','orbits','smooth_split','first_failed_character','split_character_checksum','candidate_rows']:require(a[k]==b[k],'independent quartic arithmetic '+k)
    if recompute:
        with tempfile.TemporaryDirectory(prefix='q80-k0-full-quartic-') as tmp:
            exe=Path(tmp)/'replay';subprocess.run(['g++','-O3','-std=c++17',str(ROOT/q.REPLAY),'-o',str(exe)],check=True,timeout=30)
            for (lo,hi),expected in zip(q.RANGES,aggregates[1]):
                proc=subprocess.run([str(exe),str(lo),str(hi)],input=text,text=True,capture_output=True,check=True,timeout=95);r=json.loads(proc.stdout)
                for k in expected:
                    if k!='cpu_seconds':require(r[k]==expected[k],'fresh full quartic replay '+k)
    return read(out/'quartic-comparison.json')
def verify(out,recompute=False):
    started=time.process_time();wall=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,PRODUCER,CPP,REPLAY,HELPER,QUARTIC,str(Path(__file__).relative_to(ROOT)),'tests/test_q80_genus_one_k0_gate.py',
       *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
       *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']],
       *[K1+'/'+n for n in ['input.json','result.json','fibre-hypotheses.json','independent-replay.json','cubic-input.txt']],
       *[str((out/n).relative_to(ROOT)) for n in ['quartic-input.json','quartic-input.txt','quartic-producer.json','quartic-replay.json','quartic-comparison.json']]}
    require(set(packet['bindings'])==expected,'complete frozen source set')
    for path,hsh in packet['bindings'].items():require(digest(ROOT/path)==hsh,'source binding '+path)
    require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
    for path,hsh in result['records'].items():require(digest(out/path)==hsh,'finite record binding '+path)
    source=read(ROOT/SOURCE);model=source['weierstrass_model'];points=read(ROOT/NORM4/'norm4-sections.json')['records'];require(len(points)==1313,'complete retained norm4 roster')
    rat=read(ROOT/DEGREE/'rational-fibres.json')['rows'];quad=read(ROOT/DEGREE/'quadratic-fibres.json')['rows'];comp=read(out/'composite-norms.json');check_composite(comp,rat,quad)
    contacts=check_contacts(read(out/'square-contacts.json'),points);lifts=read(out/'square-lifts.json')['rows']
    require({(tuple(r['pair']),tuple(r['H'])) for r in lifts}=={(tuple(r['pair']),tuple(H)) for r in contacts for H in r['H']} and len(lifts)==30,'every contact first lift')
    for row in lifts:
        a,b,target=first_lift(row['pair'],row['H'],points,model)
        require(a==row['rank']==24 and b==row['augmented_rank']==25 and target==row['error_target'],'independently inconsistent first lift')
    quartic=check_quartic(out,recompute)
    text=(ROOT/K1/'cubic-input.txt').read_text()+'16\n'+'\n'.join(' '.join(map(str,row)) for row in comp['rational_signatures'])+'\n'
    require(text==(out/'cubic-input.txt').read_text(),'cubic rational-signature coefficient input')
    with tempfile.TemporaryDirectory(prefix='q80-k0-cubic-independent-') as tmp:
        exe=Path(tmp)/'replay';subprocess.run(['g++','-O3','-std=c++17',str(ROOT/REPLAY),'-o',str(exe)],check=True,timeout=30)
        proc=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True,timeout=95)
    replay=json.loads(proc.stdout);cubic=read(out/'cubic-census.json')
    for k in ['status','orbits','smooth_split','nodal','signature_rows']:require(cubic[k]==replay[k],'independent cubic plus rational census '+k)
    require(cubic['status']=='COMPLETE' and cubic['orbits']==749320 and cubic['smooth_split']==125595 and cubic['nodal']==1 and not cubic['signature_rows'],'empty complete cubic plus rational gate')
    require(result['integral_pair_candidates']==0 and result['necessary_D_mod131']==mul([88,62,1],[88,62,1]),'necessary branch reduction')
    require(result['at_least_one_negative_abscissa_required'] and result['remaining_nodal_denominator_boundary']=='UNKNOWN' and result['genus1_k0_allocations_remaining']==5,'remaining unbounded boundary')
    require(result['written_integrality_and_norm_specialization_required'] and not result['formal_verification'] and not result['positive_mw17_target_complete'] and not result['old_norm8_replay_upgraded'],'assurance and positive target limits')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),'checker_sha256':digest(Path(__file__)),
            'quartic_orbits':73620690,'quartic_smooth_split':12267623,'quartic_candidates':0,'quartic_full_recomputation_this_run':recompute,
            'quartic_completed_independent_replay_bound':True,'cubic_orbits':749320,'cubic_signature_matches':0,'square_contact_pairs':30,'inconsistent_first_lifts':30,
            'python_cpu_seconds':time.process_time()-started,'cubic_replay_cpu_seconds':replay['cpu_seconds'],'wall_seconds':time.monotonic()-wall,
            'remaining_boundary':'UNKNOWN','positive_mw17_target_complete':False}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,default=DEFAULT);p.add_argument('--receipt',type=Path);p.add_argument('--replay-quartic',action='store_true');args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,95));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    record=verify(args.directory,args.replay_quartic)
    if args.receipt:
        with args.receipt.open('x') as f:json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(record),flush=True)
