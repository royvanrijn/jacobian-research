#!/usr/bin/env sage
"""Finite certificates for the Q80 genus-one k0 integral-coefficient exclusion."""
from sage.all import GF, Integers, PolynomialRing, QQ, matrix, vector
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
import argparse, importlib.util, itertools, json, resource, subprocess, tempfile, time
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
K1='artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1'
CPP='elkies-k3/scripts/certify_q80_cubic_rational_norms.cpp'
REPLAY='elkies-k3/scripts/replay_q80_cubic_rational_norms.cpp'
CHECKER='elkies-k3/scripts/verify_q80_genus_one_k0_gate.py'
TEST='tests/test_q80_genus_one_k0_gate.py'
HELPER='elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
P=131
def read(p):return json.loads(p.read_text())
def digest(p):return sha256(p.read_bytes()).hexdigest()
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def bindings():
    files=[SOURCE,CPP,REPLAY,CHECKER,TEST,HELPER,str(Path(__file__).relative_to(ROOT)),
           *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
           *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']],
           *[K1+'/'+n for n in ['input.json','result.json','fibre-hypotheses.json','independent-replay.json','cubic-input.txt']],
           'elkies-k3/scripts/run_q80_quartic_norm_census.py']
    files += [str((OUT/n).relative_to(ROOT)) for n in ['quartic-input.json','quartic-input.txt','quartic-producer.json','quartic-replay.json','quartic-comparison.json']]
    return files
def freeze():
    write('input.json',{'schema':1,'prime':P,'scope':'Complete finite quartic norm and square-contact first-lift gates; the necessary nodal denominator boundary remains UNKNOWN.',
          'driver_cpu_seconds':30,'driver_memory_bytes':2*1024**3,'cubic_cpu_seconds':90,'cubic_memory_bytes':512*1024**2,
          'bindings':{p:digest(ROOT/p) for p in bindings()}})
def composite(rat,quad):
    rs=[r for r in rat if len(r['rational_roots'])==3];qs=[r for r in quad if len(r['roots'])==3]
    def options(r,k):return set(itertools.permutations(zip(r['roots'],r[k])))
    index=defaultdict(list);count=0
    for a,b in itertools.combinations(rs,2):
        for aa in options(a,'rational_codes'):
            for bb in options(b,'rational_codes'):
                codes=tuple(sorted(aa[j][1]^bb[j][1] for j in range(3)))
                index[codes].append([a['t'],b['t']]);count+=1
    four=sum(set(a).isdisjoint(b) for rows in index.values() for a,b in itertools.combinations(rows,2))
    mixed=sum(len(index.get(tuple(sorted(q['norm_codes'])),[])) for q in qs)
    qindex=defaultdict(list)
    for q in qs:qindex[tuple(sorted(q['norm_codes']))].append(q['t'])
    two=sum(len(rows)*(len(rows)-1)//2 for rows in qindex.values())
    assert count==4320 and four==mixed==two==0
    return {'rational_pair_choices':count,'four_rational_matches':four,'two_rational_quadratic_matches':mixed,'two_quadratic_matches':two,
            'rational_signatures':[sorted(r['rational_codes']) for r in rs]}
def contacts(points,R):
    t=R.gen();xs=[R(r['x']) for r in points];ys=[R(r['y']) for r in points];found=[];shared=0;nodal=[];q0=t*t+62*t+88
    for i in range(len(points)):
        if ys[i]%q0==0:nodal.append({'index':i,'q0_order':int(ys[i].valuation(q0)),'x_mod_q0':[int(c) for c in (xs[i]%q0).list()]})
        for j in range(i):
            g=ys[i].gcd(ys[j]);inf=min(6-ys[i].degree(),6-ys[j].degree())
            if g.degree()+inf<2:continue
            shared+=1;diff=xs[i]-xs[j]
            while g.gcd(diff).degree()>0:g//=g.gcd(diff)
            if xs[i][4]==xs[j][4]:inf=0
            if g.degree()+inf<2:continue
            fac=[(f,e) for f,e in g.factor() if f.degree()<=2];hh=[]
            for f,e in fac:
                if f.degree()==2:hh.append(f)
                elif e>=2:hh.append(f*f)
            linear=[f for f,e in fac if f.degree()==1]
            hh += [a*b for a,b in itertools.combinations(linear,2)]
            if inf:hh+=linear
            if inf>=2:hh.append(R(1))
            if hh:found.append({'pair':[j,i],'H':[[int(f[k]) for k in range(3)] for f in hh]})
    assert len(found)==30 and all(len(r['H'])==1 and r['H'][0][2]==1 for r in found)
    return {'pairs_examined':len(points)*(len(points)-1)//2,'shared_degree_at_least2':shared,'distinct_root_contact_pairs':found,'sections_with_q0_ordinate_factor':nodal}
def lifts(contact,points,model,R):
    F=R.base_ring();t=R.gen();Z=Integers(P*P);R2=PolynomialRing(Z,'t')
    A=R([F(QQ(v)) for v in model['A_coefficients_low_to_high']]);A2=R2([Z(QQ(v).numerator())/Z(QQ(v).denominator()) for v in model['A_coefficients_low_to_high']])
    B2=R2([Z(QQ(v).numerator())/Z(QQ(v).denominator()) for v in model['B_coefficients_low_to_high']]);lift=lambda f:R2([int(c) for c in f.list()]);rows=[]
    for row in contact['distinct_root_contact_pairs']:
        for hh in row['H']:
            H=R(hh);D=H*H;assert H.degree()==2
            J=matrix(F,26,24);error=[]
            for i,n in enumerate(row['pair']):
                X,Y=R(points[n]['x']),R(points[n]['y']);r,rem=Y.quo_rem(H);assert not rem
                local=[(3*X*X+A)*t**k for k in range(5)]+[-2*D*r*t**k for k in range(5)]
                for j,f in enumerate(local):
                    for k in range(13):J[13*i+k,10*i+j]=f[k]
                for j in range(4):
                    f=-r*r*t**j
                    for k in range(13):J[13*i+k,20+j]=f[k]
                err=lift(X)**3+A2*lift(X)+B2-lift(D)*lift(r)**2
                assert all(int(err[k])%P==0 for k in range(13));error.extend([int(F(-int(err[k])//P)) for k in range(13)])
            rank=int(J.rank());aug=int(J.augment(vector(F,error).column()).rank());assert (rank,aug)==(24,25)
            rows.append({'pair':row['pair'],'H':hh,'rank':rank,'augmented_rank':aug,'error_target':error})
    return {'rows':rows}
def run():
    started=time.process_time();wall=time.monotonic();packet=read(OUT/'input.json');assert set(packet['bindings'])==set(bindings())
    for path,h in packet['bindings'].items():assert digest(ROOT/path)==h
    source=read(ROOT/SOURCE);model=source['weierstrass_model'];rat=read(ROOT/DEGREE/'rational-fibres.json')['rows'];quad=read(ROOT/DEGREE/'quadratic-fibres.json')['rows']
    gate=composite(rat,quad);write('composite-norms.json',gate)
    points=read(ROOT/NORM4/'norm4-sections.json')['records'];R=PolynomialRing(GF(P),'t')
    contact=contacts(points,R);write('square-contacts.json',contact);write('square-lifts.json',lifts(contact,points,model,R))
    assert R([2,0,-2,0,1]).is_irreducible()
    text=(ROOT/K1/'cubic-input.txt').read_text()+str(len(gate['rational_signatures']))+'\n'+'\n'.join(' '.join(map(str,row)) for row in gate['rational_signatures'])+'\n'
    with (OUT/'cubic-input.txt').open('x') as f:f.write(text)
    with tempfile.TemporaryDirectory(prefix='q80-k0-cubic-') as tmp:
        exe=Path(tmp)/'census';subprocess.run(['g++','-O3','-std=c++17',str(ROOT/CPP),'-o',str(exe)],check=True,timeout=30)
        proc=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True,timeout=95)
    cubic=json.loads(proc.stdout);assert cubic['status']=='COMPLETE' and cubic['orbits']==749320 and not cubic['signature_rows'];write('cubic-census.json',cubic)
    write('result.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),
          'records':{n:digest(OUT/n) for n in ['composite-norms.json','square-contacts.json','square-lifts.json','cubic-input.txt','cubic-census.json']},
          'quartic_modulus_irreducible':True,'integral_pair_candidates':0,'necessary_D_mod131':[int(c) for c in R([88,62,1])**2],
          'at_least_one_negative_abscissa_required':True,'remaining_nodal_denominator_boundary':'UNKNOWN','genus1_k0_allocations_remaining':5,
          'written_integrality_and_norm_specialization_required':True,'positive_mw17_target_complete':False,'formal_verification':False,'old_norm8_replay_upgraded':False})
    write('execution.json',{'status':'PASS','sage_cpu_seconds':time.process_time()-started,'cubic_cpu_seconds':cubic['cpu_seconds'],'wall_seconds':time.monotonic()-wall})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','run']);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,95));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    freeze() if a.mode=='freeze' else run()
