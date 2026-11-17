#!/usr/bin/env python3
"""Exact abelian-cover homology, deriving affine quotient primes from torsion."""
import hashlib
import json
from pathlib import Path
import time
from sympy import primerange,factorint,isprime
from integer_lifts import smith_certificate
from verify_integer_lifts import nullspace

HERE=Path(__file__).resolve().parent


def homology_matrix(words,N,r):
    rows=[]
    for word in words:
        for start in range(N):
            sheet=start;row=[0]*(N+1)
            for a in word:
                shift=r if abs(a)==1 else 1
                if a>0:
                    index=sheet if abs(a)==1 else N if sheet==N-1 else None
                    if index is not None:row[index]+=1
                    sheet=(sheet+shift)%N
                else:
                    sheet=(sheet-shift)%N
                    index=sheet if abs(a)==1 else N if sheet==N-1 else None
                    if index is not None:row[index]-=1
            assert sheet==start;rows.append(row)
    return rows


def affine_matrix(words,N,r):
    rows=[]
    for word in words:
        block=[[0]*N for _ in range(N)];prefix=0
        for a in word:
            shift=r if abs(a)==1 else 1
            if a<0:prefix=(prefix-shift)%N
            if abs(a)==1:
                for j in range(N):block[(j+prefix)%N][j]+=1 if a>0 else -1
            if a>0:prefix=(prefix+shift)%N
        assert prefix==0;rows.extend(block)
    rows.append([1]*N)
    return rows


def shift(v,s):return [v[(j-s)%len(v)] for j in range(len(v))]
def multiply(a,b,p):return ((a[0]+b[0])%len(a[1]),[(x+y)%p for x,y in zip(a[1],shift(b[1],a[0]))])
def inverse(a,p):return (-a[0]%len(a[1]),[-x%p for x in shift(a[1],-a[0])])
def evaluate(word,assignment,p):
    N=len(next(iter(assignment.values()))[1]);value=(0,[0]*N)
    for a in word:value=multiply(value,assignment[a] if a>0 else inverse(assignment[-a],p),p)
    return value


def main():
    plan=json.loads((HERE/'three_sheet_metabelian_plan.json').read_text())
    structure=json.loads((HERE/'three_sheet_structural_results.json').read_text())
    verified=json.loads((HERE/'three_sheet_structure_verification.json').read_text());assert verified['status']=='PASS'
    summaries={r['name']:r for r in verified['results']};old={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    rows=[];started=time.monotonic()
    for source in structure['results']:
        if len(source['active_generators'])!=2:continue
        assert time.monotonic()-started<plan['limits']['seconds_total']
        name=source['name'];N=summaries[name]['abelianization_order'];words=source['final_relators']
        possible=[r for r in range(N) if all(sum((1 if a>0 else -1)*(r if abs(a)==1 else 1) for a in w)%N==0 for w in words)]
        assert len(possible)==1;r=possible[0]
        matrix=homology_matrix(words,N,r);smith=smith_certificate(matrix,N+1)
        torsion=[x for x in smith['diagonal'] if x>1];free=N+1-len(smith['diagonal'])
        record={'name':name,'abelianization_order':N,'generator_shifts':[r,1],'cyclic_cover_matrix':matrix,'smith':smith,
                'H1_free_rank':free,'H1_torsion':torsion,'affine_tests':[]}
        primes=[]
        if free:primes=[2]
        elif torsion:
            primes=sorted({p for p in primerange(2,1000) if any(x%p==0 for x in torsion) and N%p})
            if not primes:
                factors=factorint(torsion[-1],limit=10000)
                primes=sorted(p for p in factors if isprime(p) and N%p)
        for p in primes:
            if N%p==0:continue
            cocycle_rows=affine_matrix(words,N,r);basis=nullspace(cocycle_rows,N,int(p))
            record['affine_tests'].append({'prime':int(p),'gauge':'y_translation0, x_augmentation0','matrix':cocycle_rows,'kernel_basis':basis})
            if not basis:continue
            u=basis[0];assignment={1:(r,u),2:(1,[0]*N)}
            assert multiply(assignment[1],assignment[2],p)!=multiply(assignment[2],assignment[1],p)
            assert all(evaluate(w,assignment,p)==(0,[0]*N) for w in words)
            row=old[name];assignment={row['remaining_generators'][j]:evaluate(source['original_generator_words'][str(j+1)],assignment,p) for j in range(2)}
            for step in reversed(row['eliminations']):assignment[step['generator']]=evaluate(step['replacement'],assignment,p)
            assert all(evaluate(w,assignment,p)==(0,[0]*N) for w in row['relators'])
            pair=next((a,b) for a in assignment for b in assignment if multiply(assignment[a],assignment[b],p)!=multiply(assignment[b],assignment[a],p))
            record['nonabelian_witness']={'prime':int(p),'regular_module_dimension':N,'final_generator_translation':u,
                                          'original_generator_images':assignment,'noncommuting_original_generators':pair}
            record['sharp_four_candidate']=True;break
        rows.append(record)
        report={'status':'IN_PROGRESS','results':rows,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['three_sheet_metabelian_plan.json','audit_three_sheet_metabelian.py','three_sheet_structural_results.json','three_sheet_structure_verification.json','three_sheet_nonabelian_results.json','integer_lifts.py','verify_integer_lifts.py']}}
        (HERE/'three_sheet_metabelian_results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'name':name,'cyclic_degree':N,'H1_free':free,'H1_torsion':torsion,'prime':record.get('nonabelian_witness',{}).get('prime'),'sharp_candidate':record.get('sharp_four_candidate',False)}),flush=True)
    report['status']='COMPLETE';(HERE/'three_sheet_metabelian_results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
