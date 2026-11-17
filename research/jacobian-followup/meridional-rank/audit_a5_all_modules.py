#!/usr/bin/env python3
"""Frozen full cover cohomology for the untested A5 modules."""
import hashlib
import json
from pathlib import Path
import time
from integer_lifts import smith_certificate
from parabolic_lifts import oriented_diagram
from quadratic_integer_lifts import ZERO,ONE,IDENTITY,add,neg,mul,matmul,h3_model,integral_expansion,factor_exact,residue_fields,kernel

HERE=Path(__file__).resolve().parent


def perm_inverse(p):return tuple(p.index(i) for i in range(len(p)))


def cocycle_matrix(images,inverses,words):
    n=len(images);out=[]
    for word in words:
        prefix=IDENTITY;block=[[ZERO]*(3*n) for _ in range(3)]
        for letter in word:
            i=abs(letter)-1
            if letter<0:prefix=matmul(prefix,inverses[i])
            for j in range(3):
                for k in range(3):
                    entry=prefix[3*j+k] if letter>0 else neg(prefix[3*j+k])
                    block[j][3*i+k]=add(block[j][3*i+k],entry)
            if letter>0:prefix=matmul(prefix,images[i])
        assert prefix==IDENTITY
        out.extend(block)
    return out


def principal_matrix(images):
    return [[add(m[3*i+j],neg(ONE if i==j else ZERO)) for j in range(3)] for m in images for i in range(3)]


def dimensions(matrix,principal,field):
    z=len(kernel(matrix,len(matrix[0]),field))
    inv=len(kernel(principal,3,field))
    return {'prime':field.p,'tau_root':field.root,'field_degree':2 if field.root is None else 1,
            'dimension':3,'Z1_dimension':z,'invariants':inv,'H1_dimension':z-3+inv,
            'rank_lower_bound':(z+2)//3,'sharp_four':z>6}


def quotient_two(images,field):
    assert field.p==2
    # W3 has a one-dimensional invariant line in characteristic2.
    fixed=kernel(principal_matrix(images),3,field)
    assert len(fixed)==1
    line=fixed[0];pivot=next(i for i,x in enumerate(line) if x!=ZERO)
    free=[i for i in range(3) if i!=pivot]
    scale=field.inv(line[pivot]);line=[field.mul(x,scale) for x in line]
    matrices=[]
    for m in images:
        a=[]
        for i in free:
            for j in free:
                a.append(field.add(field.reduce(m[3*i+j]),field.neg(field.mul(line[i],field.reduce(m[3*pivot+j])))))
        matrices.append(a)
    return matrices,line,pivot


def small_mm(a,b,n,field):
    out=[]
    for i in range(n):
        for j in range(n):
            x=ZERO
            for k in range(n):x=field.add(x,field.mul(a[n*i+k],b[n*k+j]))
            out.append(x)
    return out


def small_inverse(a,n,field):
    rows=[[field.reduce(a[n*i+j]) for j in range(n)]+[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
    for j in range(n):
        i=next(i for i in range(j,n) if rows[i][j]!=ZERO)
        rows[j],rows[i]=rows[i],rows[j];s=field.inv(rows[j][j])
        rows[j]=[field.mul(s,x) for x in rows[j]]
        for i in range(n):
            if i!=j:
                c=rows[i][j];rows[i]=[field.add(x,field.neg(field.mul(c,y))) for x,y in zip(rows[i],rows[j])]
    return [x for row in rows for x in row[n:]]


def small_cocycles(images,words,n,field):
    inverses=[small_inverse(a,n,field) for a in images]
    unit=[ONE if i==j else ZERO for i in range(n) for j in range(n)]
    rows=[]
    for word in words:
        prefix=unit;block=[[ZERO]*(n*len(images)) for _ in range(n)]
        for letter in word:
            i=abs(letter)-1
            if letter<0:prefix=small_mm(prefix,inverses[i],n,field)
            for j in range(n):
                for k in range(n):block[j][n*i+k]=field.add(block[j][n*i+k],prefix[n*j+k] if letter>0 else field.neg(prefix[n*j+k]))
            if letter>0:prefix=small_mm(prefix,images[i],n,field)
        assert prefix==unit
        rows.extend(block)
    principal=[[field.add(a[n*i+j],field.neg(ONE if i==j else ZERO)) for j in range(n)] for a in images for i in range(n)]
    z=len(kernel(rows,n*len(images),field));invariants=len(kernel(principal,n,field))
    return {'dimension':n,'images':images,'cocycle_matrix':rows,'principal_matrix':principal,
            'Z1_dimension':z,'invariants':invariants,'H1_dimension':z-n+invariants,
            'rank_lower_bound':(z+n-1)//n,'sharp_four':z>2*n}


def main():
    plan=json.loads((HERE/'a5_all_modules_plan.json').read_text());assert plan['frozen_before_execution']
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases'][:16]
    maps=json.loads((HERE/'joint_a5_results.json').read_text())['results'][:16]
    catalogue=json.loads((HERE/'joint_a5_catalogue.json').read_text())
    model,_,group_certificate=h3_model()
    paths=['a5_all_modules_plan.json','audit_a5_all_modules.py','joint_a5_panel.json','joint_a5_results.json',
           'joint_a5_catalogue.json','quadratic_integer_lifts.py','integer_lifts.py','parabolic_lifts.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    output=HERE/'a5_all_modules_results.json'
    report=json.loads(output.read_text()) if output.exists() else {'status':'RUNNING','results':[],
        'input_sha256':hashes,'H3_group_certificate':group_certificate}
    assert report['input_sha256']==hashes
    done={(r['name'],r['orbit_index']) for r in report['results']};start=time.monotonic()
    for case,source in zip(panel,maps):
        assert case['name']==source['name']
        words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in oriented_diagram(case['pd'])]+[[source['normalized_arc']+1]]
        for mapping in source['epimorphism_representatives']:
            key=(case['name'],mapping['orbit_index'])
            if key in done:continue
            assert time.monotonic()-start<plan['limits']['seconds_total'];stamp=time.monotonic()
            elements=[tuple(catalogue['elements'][i]) for i in mapping['arc_images']]
            images=[model[(p,0)] for p in elements]
            inverses=[model[(perm_inverse(p),0)] for p in elements]
            matrix=cocycle_matrix(images,inverses,words);principal=principal_matrix(images)
            assert all(sum_entries([mul(x,principal[j][k]) for j,x in enumerate(r)])==ZERO for r in matrix for k in range(3))
            smith=smith_certificate(integral_expansion(matrix),39*2)
            invsmith=smith_certificate(integral_expansion(principal),3*2)
            assert len(smith['diagonal'])%2==len(invsmith['diagonal'])%2==0
            generic_z=39-len(smith['diagonal'])//2;generic_inv=3-len(invsmith['diagonal'])//2
            candidates=sorted(set([2,3,5]+factor_exact(smith['diagonal'][-1],1000000)+factor_exact(invsmith['diagonal'][-1],1000000)))
            fields=[dimensions(matrix,principal,f) for p in candidates for f in residue_fields(p)]
            binary=[]
            for f in residue_fields(2):
                two,line,pivot=quotient_two(images,f)
                item=small_cocycles(two,words,2,f)
                item.update({'prime':2,'tau_root':f.root,'field_degree':2 if f.root is None else 1,
                             'invariant_line':line,'pivot':pivot})
                binary.append(item)
                # The other simple is obtained by Frobenius field conjugation.
                conjugate=[[f.mul(x,x) for x in a] for a in two]
                item=small_cocycles(conjugate,words,2,f)
                item.update({'prime':2,'frobenius_conjugate':True,'tau_root':f.root,
                             'field_degree':2 if f.root is None else 1})
                binary.append(item)
            result={'name':key[0],'orbit_index':key[1],'source_group':mapping['source_group'],'images':images,
                    'words':words,'cocycle_matrix':matrix,'principal_matrix':principal,
                    'integer_smith':smith,'invariant_smith':invsmith,
                    'generic_Z1':generic_z,'generic_invariants':generic_inv,
                    'candidate_primes':candidates,'field_dimensions':fields,'binary_simple_modules':binary,
                    'sharp_four':generic_z>6 or any(r['sharp_four'] for r in fields+binary),
                    'seconds':time.monotonic()-stamp}
            report['results'].append(result);output.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'name':key[0],'orbit_index':key[1],'generic_Z1':generic_z,
                'primes':candidates,'field_dims':[(r['prime'],r['tau_root'],r['Z1_dimension'],r['invariants']) for r in fields],
                'binary_Z1':[r['Z1_dimension'] for r in binary],'sharp_four':result['sharp_four'],'seconds':result['seconds']}),flush=True)
            assert result['seconds']<plan['limits']['seconds_per_map']
    assert len(report['results'])==16;report['status']='COMPLETE'
    output.write_text(json.dumps(report,indent=2)+'\n')


def sum_entries(values):
    result=ZERO
    for v in values:result=add(result,v)
    return result


if __name__=='__main__':main()
