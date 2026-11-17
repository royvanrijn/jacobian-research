#!/usr/bin/env python3
"""Frozen projective13 subgroup homology and exceptional-characteristic test."""

from collections import Counter
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import signal
import time

from audit_core_subgroup_homology import covering_matrix
from integer_lifts import smith_certificate
from ternary_transvections import image,inverse
from verify_integer_lifts import nullspace

HERE=Path(__file__).resolve().parent


def canonical_line(v):
    first=next(x for x in v if x);scale=1 if first==1 else 2
    return tuple(scale*x%3 for x in v)


LINES=sorted({canonical_line(v) for v in product(range(3),repeat=3) if any(v)})


def projective_action(matrices):
    return [tuple(LINES.index(canonical_line(image(a,v))) for v in LINES) for a in matrices]


def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def invert(p):return tuple(p.index(i) for i in range(len(p)))


def intertwiner(source,target):
    for target_zero in range(13):
        mapping={0:target_zero};queue=[0];valid=True
        for v in queue:
            for a,b in zip(source,target):
                w=a[v];z=b[mapping[v]]
                if w in mapping:
                    if mapping[w]!=z:valid=False;break
                else:mapping[w]=z;queue.append(w)
            if not valid:break
        if valid and len(mapping)==len(set(mapping.values()))==13:return [mapping[i] for i in range(13)]
    return None


def exponent_matrix(words,n):
    out=[]
    for word in words:
        row=[0]*n
        for letter in word:row[abs(letter)-1]+=1 if letter>0 else -1
        out.append(row)
    return out


def strip13(value):
    while value%13==0:value//=13
    return value


def maximum_gain_cyclic(D,free,torsion):
    # p first divides torsion[j] => its cover multiplicity is len(torsion)-j.
    # Separate p|D from p not dividingD, without factoring large integers.
    candidates=[{'gain':free,'kind':'generic prime'}]
    for j,d in enumerate(torsion):
        count=len(torsion)-j;common=strip13(gcd(d,D))
        if common>1:candidates.append({'gain':free+count-1,'kind':'prime dividing baseD','divisor':common,'first_candidate_index':j})
        remainder=strip13(d)
        while (common:=gcd(remainder,D))>1:remainder//=common
        if remainder>1:candidates.append({'gain':free+count,'kind':'prime coprime to baseD','divisor':remainder,'first_candidate_index':j})
    return max(c['gain'] for c in candidates),candidates


def augmentation(p):
    return [[int(p[i]==j)-int(p[12]==j) for i in range(12)] for j in range(12)]


def augmentation_cocycles(images,words,prime):
    n=len(images);rows=[];ident=tuple(range(13))
    for word in words:
        prefix=ident;block=[[0]*(12*n) for _ in range(12)]
        for letter in word:
            i=abs(letter)-1
            if letter<0:prefix=compose(prefix,invert(images[i]))
            coefficient=augmentation(prefix);sign=1 if letter>0 else -1
            for a in range(12):
                for b in range(12):block[a][12*i+b]+=sign*coefficient[a][b]
            if letter>0:prefix=compose(prefix,images[i])
        assert prefix==ident;rows.extend(block)
    basis=nullspace(rows,12*n,prime)
    principal=[[(v-int(a==b))%prime for b,v in enumerate(row)] for p in images for a,row in enumerate(augmentation(p))]
    fixed=len(nullspace(principal,12,prime))
    return {'prime':prime,'rows':len(rows),'columns':12*n,'Z1_dimension':len(basis),
        'H1_dimension':len(basis)-12+fixed,'fixed_dimension':fixed,'kernel_basis':basis,
        'matrix_sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
        'rank3_gate_passed':len(basis)>=25}


def homology(images,words):
    matrix,tree=covering_matrix(images,words);columns=len(matrix[0]);smith=smith_certificate(matrix,columns)
    return {'generator_images':images,'presentation_words':words,'tree_edges':tree,'rows':len(matrix),
        'columns':columns,'smith':smith,'H1_free_rank':columns-len(smith['diagonal']),
        'H1_torsion_factors':[d for d in smith['diagonal'] if d>1]}


def timeout_handler(*_):raise TimeoutError('frozen per-cover time cap')


def main():
    started=time.monotonic();plan=json.loads((HERE/'projective13_cover_plan.json').read_text())
    source=json.loads((HERE/'ternary_cover_module_results.json').read_text())
    verified=json.loads((HERE/'ternary_cover_module_verification.json').read_text())
    assert verified['status']=='PASS independent natural-dual ternary cover modules'
    for p,digest in verified['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    files=['projective13_cover_plan.json','audit_projective13_covers.py','audit_core_subgroup_homology.py',
        'integer_lifts.py','ternary_transvections.py','verify_integer_lifts.py',
        'ternary_cover_module_results.json','ternary_cover_module_verification.json','ternary_transvection_control.json']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}
    path=HERE/'projective13_cover_results.json';report={'status':'RUNNING','input_sha256':hashes,'controls':[],'base_homology':{},'results':[]}
    if path.exists():report=json.loads(path.read_text());assert report['input_sha256']==hashes
    signal.signal(signal.SIGALRM,timeout_handler)
    if not report['controls']:
        control=json.loads((HERE/'ternary_transvection_control.json').read_text())
        matrices=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in control['chain_matrices'][:3]]
        dual=[tuple(inverse(a)[3*j+i] for i in range(3) for j in range(3)) for a in matrices]
        for label,model in [('point',matrices),('dual_point',dual)]:
            images=projective_action(model);words=[[1,1,1],[2,2,2],[3,3,3]]
            h=homology(images,words);assert h['H1_free_rank']==6 and h['H1_torsion_factors']==[3]*12
            at3=augmentation_cocycles(images,words,3);assert at3['Z1_dimension']==27
            report['controls'].append({'action':label,**h,'base_F3_betti':3,'cover_F3_betti':18,'gain_F3':15,
                'augmentation_F3':at3,'augmentation_F13':augmentation_cocycles(images,words,13)})
        path.write_text(json.dumps(report,indent=2)+'\n')
    done={(r['name'],r['representation_index'],r['action']) for r in report['results']}
    for row in source['results']:
        name=row['name'];words=row['core_words'];n=len(row['modules'][0]['matrices'])
        if name not in report['base_homology']:
            base=smith_certificate(exponent_matrix(words,n),n);torsion=[d for d in base['diagonal'] if d>1]
            assert n==len(base['diagonal']) and len(torsion)==1
            report['base_homology'][name]={'smith':base,'H1_free_rank':0,'H1_torsion_factors':torsion}
        D=report['base_homology'][name]['H1_torsion_factors'][0]
        for module in row['modules']:
            action='point' if module['module']=='natural3' else 'dual_point'
            key=(name,row['representation_index'],action)
            if key in done:continue
            if time.monotonic()-started>plan['limits']['seconds_total']:break
            images=projective_action(module['matrices']);alias=None
            for i,previous in enumerate(report['results']):
                if previous['name']!=name or previous['status']!='COMPLETE' or 'equivalent_to' in previous:continue
                mapping=intertwiner(images,previous['generator_images'])
                if mapping is not None:alias=(i,mapping);break
            result={'name':name,'representation_index':row['representation_index'],'action':action,
                'source_group':f'pi1(Sigma2({name}));core-root0','generator_images':images,'presentation_words':words}
            if alias is not None:
                previous=report['results'][alias[0]]
                result.update({'status':'COMPLETE','equivalent_to':alias[0],'sheet_intertwiner':alias[1],
                    'maximum_gain_all_primes_except13':previous['maximum_gain_all_primes_except13'],
                    'char13_Z1_dimension':previous['augmentation_F13']['Z1_dimension'],
                    'sharp_four':previous['sharp_four']})
            else:
                signal.alarm(plan['limits']['seconds_per_distinct_cover'])
                try:
                    result.update(homology(images,words));gain,candidates=maximum_gain_cyclic(D,result['H1_free_rank'],result['H1_torsion_factors'])
                    at13=augmentation_cocycles(images,words,13)
                    result.update({'status':'COMPLETE','maximum_gain_all_primes_except13':gain,'gain_candidates':candidates,
                        'augmentation_F13':at13,'base_D':D,
                        'sharp_four':gain>=13 or at13['rank3_gate_passed']})
                except TimeoutError as error:result.update({'status':'CAPPED_INCOMPLETE','reason':str(error)})
                finally:signal.alarm(0)
            report['results'].append(result);path.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({k:result[k] for k in ['name','representation_index','action','status','H1_free_rank','H1_torsion_factors','equivalent_to','maximum_gain_all_primes_except13','char13_Z1_dimension','sharp_four'] if k in result}),flush=True)
    report['status']='COMPLETE' if len(report['results'])==32 and all(r['status']=='COMPLETE' for r in report['results']) else 'CAPPED_INCOMPLETE'
    report['distinct_source_actions']=sum('equivalent_to' not in r for r in report['results'])
    report['sharp_names']=sorted({r['name'] for r in report['results'] if r.get('sharp_four')})
    report['seconds']=round(time.monotonic()-started,6);path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','distinct_source_actions','sharp_names','seconds']},indent=2))


if __name__=='__main__':main()
