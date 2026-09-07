#!/usr/bin/env python3
"""Record-scale million-parameter intake after the302 half-gap calibration."""
import argparse,gzip,hashlib,heapq,json,math,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import certify_compact_r17_candidates as cert
import score_retained_native19 as scalar
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/det1092-record-scale-selection-v1'
SOURCE=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
GATE=ART/'curve302_recovery_calibration_v1.json'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),CAS/'prepare_det1092_record_scale_tables.sage',Path(scalar.__file__),CAS/'research_runtime/finite_fields.py']}

def freeze():
    assert not (D/'protocol.json').exists();g=cert.read(GATE)
    assert g['status']=='PASS' and g['significant_fraction_gate_met'] and g['recovered_exceptional_directions']>=7
    checkpoint(D/'protocol.json',dict(schema='elliptic-curves.det1092-record-scale-selection.v1',sources=sources(),
        inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,GATE]},
        population=1048576,maximum_draws=2097152,parameter_height_min=65536,parameter_height_max=1048576,
        sample_domain='det1092-record-scale-population-v1',block_size=16384,
        height_bins=[10,11],height_bin_width=64,retain_per_bin=1024,select_per_bin=24,
        first_primes=[p for p in scalar.PRIMES if p<=997],extension_primes=scalar.PRIMES,
        validation_primes='65537..131071 are neither computed nor read',
        gp_sha256=cert.hashed(scalar.GP),maximum_workers=1,rss_bytes=2147483648,
        table_seconds=180,population_seconds=600,extension_seconds=900,seconds_per_score=30,
        gate='Completed original17-only and recovered-only302 calibration reaches24 with seven exactly identified directions in the14-dimensional public exceptional quotient. Public points enter only its terminal identity audit, never this prospective selector or geometry.',
        population_policy='SHA256(domain:counter) supplies m from bytes0:8 modulo2H+1 minusH and n from bytes8:16 moduloH plus1. Reduce by gcd; discard zero, height below65536 and repeated fractions. Take1048576 distinct addresses, with at most2097152 draws. No known302 coordinate, exceptional point, rank label or catalogue lookup chooses addresses.',
        selection='Evaluate every sampled equation and exact reduced j height. Score only640..767-bit j numerators, the user-requested302-scale region. Independently retain1024 largest first-stage scores in each64-bit band; ties smaller denominator,numerator,id. Extend all retained rows through32749 with corrected good-prime scores. Select16 strongest, four consecutive rows starting at floor(N/3), and four SHA-selected rows from the bottom third in each band. Exactly48 distinct selected addresses; no point-based replacement.',
        score='Sum round(1e12*(2-ap)*log(p)/(p+1-ap)) for good minimal short reductions. Projective tables handle unambiguous residues; every A=B=0 residue is recomputed on exact integer homogeneous coefficients with repeated p4/p6 division before good/bad classification. No omitted5/13 good contributions. All table entries independently replay against PARI; every retained first-stage score cross-checks against the corrected scalar extension.',
        point_policy='Immediately expose all48 selected fibres from their exact specialized generic17 seed using the successful adaptive policy: initial49 generic-only charts, then at most four additional49-chart waves, continuing only after an independently certified gain. Require centre parity outside the preceding input subgroup. Sample2048 fixed SHA parities, use384-bit rounded metric and49 largest computed norms, factor-free maps, height125000 and ten seconds per chart. No partial-wave rank stop. Maximum11760 point boxes. One worker;180-second map stages,1200-second point/history stages,180-second individual certificate stages. Failures are retained and never cause refill.',
        scope='Primary search family is the recovered determinant1092 MW17 parent. Comparable arithmetic height is an explicit user-requested scheduling preference, not record-coordinate seeding or a rank guarantee. No automatic second population.'))
    print('FROZEN1048576 independent addresses;48 adaptive point targets',flush=True)

def protocol():
    p=cert.read(D/'protocol.json');assert p['sources']==sources() and cert.hashed(scalar.GP)==p['gp_sha256']
    assert all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p

def coefficients():
    parent=cert.read(SOURCE);assert all(v['denominator']==['1'] for v in parent['a_invariants'])
    assert all(not v['numerator'] for v in parent['a_invariants'][:3])
    a,b=[list(map(int,v['numerator'])) for v in parent['a_invariants'][3:]]
    assert len(a)==9 and len(b)==13;return a,b

def homogeneous(co,m,n):
    out=co[-1];power=n
    for v in reversed(co[:-1]):out=out*m+v*power;power*=n
    return out

def population():
    p=protocol();tab=cert.read(D/'tables.json');assert tab['status']=='PASS' and tab['primes']==p['first_primes']
    assert tab['protocol_sha256']==cert.hashed(D/'protocol.json')
    assert not (D/'population-result.json').exists() and not (D/'population-progress.json').exists()
    tables=[]
    for t in tab['tables']:
        q=t['prime'];inv=np.array([0]+[pow(i,-1,q) for i in range(1,q)],dtype=np.int64)
        weights=np.array([0 if a is None else round((2-a)*math.log(q)/(q+1-a)*10**12) for a in t['traces']],dtype=np.int64)
        tables.append((q,inv,weights,np.array(t['ambiguous'],dtype=np.bool_)))
    ac,bc=coefficients();seen=set();counter=count=0;heaps={k:[] for k in p['height_bins']};counts={};blocks=[];fallback={};fallback_counts={}
    while count<p['population']:
        rows=[]
        while len(rows)<min(p['block_size'],p['population']-count):
            assert counter<p['maximum_draws'];draw=counter
            h=hashlib.sha256((p['sample_domain']+':'+str(counter)).encode()).digest();counter+=1;H=p['parameter_height_max']
            m=int.from_bytes(h[:8],'big')%(2*H+1)-H;n=int.from_bytes(h[8:16],'big')%H+1
            g=math.gcd(m,n);m//=g;n//=g
            if not m or max(abs(m),n)<p['parameter_height_min'] or (m,n) in seen:continue
            seen.add((m,n));a=homogeneous(ac,m,n);b=homogeneous(bc,m,n);a3=a*a*a;disc=4*a3+27*b*b
            if not disc:raise ArithmeticError('singular sampled address; preserve intake')
            jnum=6912*a3;jg=math.gcd(jnum,disc);jb=abs(jnum//jg).bit_length();band=jb//64
            rows.append(dict(index=count+len(rows),draw=draw,m=m,n=n,A=a,B=b,j_numerator_bits=jb,j_denominator_bits=abs(disc//jg).bit_length(),height_bin=band))
            counts[str(band)]=counts.get(str(band),0)+1
        admitted=[r for r in rows if r['height_bin'] in heaps];nums=np.array([r['m'] for r in admitted],dtype=np.int64);dens=np.array([r['n'] for r in admitted],dtype=np.int64);scores=np.zeros(len(admitted),dtype=np.int64)
        for q,inv,weights,amb in tables:
            den=dens%q;ix=(nums%q)*inv[den]%q;ix[den==0]=q;scores+=weights[ix]
            for j in np.flatnonzero(amb[ix]):
                a,b=admitted[j]['A'],admitted[j]['B'];k=0;q4=q**4;q6=q**6
                while a%q4==0 and b%q6==0:a//=q4;b//=q6;k+=1
                if not k:continue
                av,bv=a%q,b%q;bad=(4*av**3+27*bv*bv)%q==0
                key=(q,av,bv)
                if key not in fallback:
                    trace=None if bad else scalar.direct(av,bv,q)
                    fallback[key]=0 if trace is None else round((2-trace)*math.log(q)/(q+1-trace)*10**12)
                scores[j]+=fallback[key]
                label=str(q)+('-bad' if bad else '-good');fallback_counts[label]=fallback_counts.get(label,0)+1
        score_by_index={r['index']:int(v) for r,v in zip(admitted,scores)}
        for r,s in zip(admitted,scores):
            score=int(s);key=(score,-r['n'],-r['m'],-r['index']);heap=heaps[r['height_bin']]
            if len(heap)<p['retain_per_bin'] or key>heap[0][0]:
                item=(key,dict(id='scale-'+str(r['index']).zfill(7),family='det1092-reduced',parameter=str(F(r['m'],r['n'])),
                    model=['0','0','0',str(r['A']),str(r['B'])],model_coefficient_bits=max(abs(r['A']).bit_length(),abs(r['B']).bit_length()),
                    j_numerator_bits=r['j_numerator_bits'],j_denominator_bits=r['j_denominator_bits'],height_bin=r['height_bin'],first_score_units=score,draw=r['draw']))
                if len(heap)<p['retain_per_bin']:heapq.heappush(heap,item)
                else:heapq.heapreplace(heap,item)
        record=[[r['index'],r['draw'],r['m'],r['n'],r['j_numerator_bits'],r['j_denominator_bits'],score_by_index.get(r['index'])] for r in rows]
        path=D/'population-blocks'/('block-'+str(len(blocks)).zfill(3)+'.json.gz');path.parent.mkdir(parents=True,exist_ok=True)
        assert not path.exists();path.write_bytes(gzip.compress(json.dumps(record,separators=(',',':')).encode(),mtime=0))
        count+=len(rows);blocks.append(dict(path=str(path.relative_to(ROOT)),sha256=cert.hashed(path),rows=len(rows)))
        checkpoint(D/'population-progress.json',dict(status='RUNNING',parameters=count,draws=counter,height_bin_counts=counts,blocks=blocks,scaled_reduction_counts=fallback_counts))
        print('POPULATION',count,'/',p['population'],'admitted',sum(counts.get(str(k),0) for k in heaps),flush=True)
    retained=[]
    for band,heap in heaps.items():
        assert len(heap)==p['retain_per_bin']
        retained += [dict(r,first_stage_position=i) for i,(key,r) in enumerate(sorted(heap,reverse=True))]
    checkpoint(D/'population-result.json',dict(status='PASS',parameters=count,draws=counter,height_bin_counts=counts,blocks=blocks,
        scaled_reduction_counts=fallback_counts,retained=retained,protocol_sha256=cert.hashed(D/'protocol.json'),table_sha256=cert.hashed(D/'tables.json')))

def extend():
    p=protocol();pop=cert.read(D/'population-result.json');assert pop['status']=='PASS' and pop['parameters']==p['population']
    assert not (D/'selection-result.json').exists();scalar.D=D/'extended-scores';rows=[]
    for r in pop['retained']:
        v=scalar.evaluate(r,True)
        score=sum(round((2-a)*math.log(q)/(q+1-a)*10**12) for q,a in v['traces'] if q<=997 and a is not None)
        assert score==r['first_score_units'];v.pop('traces');rows.append(v)
        checkpoint(D/'extension-progress.json',dict(status='RUNNING',completed=len(rows),rows=rows))
        if len(rows)%64==0:print('EXTENDED',len(rows),'/',len(pop['retained']),flush=True)
    selected=[];orders={}
    for band in p['height_bins']:
        order=sorted([r for r in rows if r['height_bin']==band],key=lambda r:(-r['score_units'],F(r['parameter']).denominator,F(r['parameter']).numerator,r['id']));N=len(order)
        assert N==1024;orders[str(band)]=[r['id'] for r in order]
        selected += [dict(r,stratum='strong') for r in order[:16]]
        selected += [dict(r,stratum='moderate') for r in order[N//3:N//3+4]]
        lower=sorted(order[2*N//3:],key=lambda r:hashlib.sha256(('det1092-record-scale-lower-v1:'+r['id']).encode()).hexdigest())
        selected += [dict(r,stratum='lower_fixed') for r in lower[:4]]
    assert len(selected)==48 and len({r['id'] for r in selected})==48
    checkpoint(D/'selection-result.json',dict(status='PASS',protocol_sha256=cert.hashed(D/'protocol.json'),
        population_sha256=cert.hashed(D/'population-result.json'),rows=rows,orders=orders,selected=selected,
        scope='Completed million-address record-scale selection. First-stage scores replay against exact corrected full scalar traces for all2048 retained equations.48 fixed prospective fibres; no point, validation-prime or catalogue input. No rank inference.'))
    print('SELECTED48 record-scale fibres for calibrated adaptive exposure',flush=True)

def launch():
    p=protocol();out=D/'ledger.json';assert not out.exists();ledger=dict(status='RUNNING',stages=[]);checkpoint(out,ledger)
    jobs=[('tables',[SAGE,str(CAS/'prepare_det1092_record_scale_tables.sage')],p['table_seconds']),
          ('population',[sys.executable,str(Path(__file__).resolve()),'population'],p['population_seconds']),
          ('extend',[sys.executable,str(Path(__file__).resolve()),'extend'],p['extension_seconds'])]
    for name,cmd,seconds in jobs:
        s=run(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
        print(name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
        if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed record-scale selection stage')
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','launch','population','extend']);a=ap.parse_args();globals()[a.stage]()
