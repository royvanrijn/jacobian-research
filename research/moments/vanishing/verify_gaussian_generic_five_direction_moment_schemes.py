"""Verify generic Gaussian five-direction determinant witnesses.

Verification reconstructs literal integer Macaulay minors. It calls no
Groebner generator, RREF, rational reconstruction or ideal-membership
oracle. Generic Hilbert conclusions additionally use the written
complete-intersection dimension argument.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial, gcd, isqrt, comb
from pathlib import Path
import platform
import time
from flint import nmod_mat
import flint

if not __debug__:raise RuntimeError('Exact assertions must be enabled')

PRIME=32003
HILBERT=[1,5,14,29,49,71,90,100,95,70,22,0]
ZERO=(0,0,0,0,0)

@lru_cache(None)
def monomials(degree,variables=5):
    if variables==1:
        return ((degree,),)
    return tuple((a,)+tail for a in range(degree,-1,-1)
                 for tail in monomials(degree-a,variables-1))

@lru_cache(None)
def centered_exponential(n):
    if n==0:return 1
    if n==1:return 0
    return (n-1)*(centered_exponential(n-1)+centered_exponential(n-2))

@lru_cache(None)
def gaussian(n):
    return 0 if n%2 else factorial(n)//(2**(n//2)*factorial(n//2))

def multiply_univariate(a,b):
    result={}
    for i,c in a.items():
        for j,d in b.items():result[i+j]=result.get(i+j,0)+c*d
    return {i:c for i,c in result.items() if c}

@lru_cache(None)
def hermite_product(i,k):
    result={0:1}
    for _ in range(i):result=multiply_univariate(result,{3:1,1:-3})
    for _ in range(k):result=multiply_univariate(result,{2:1,0:-1})
    return result

@lru_cache(None)
def moment(n):
    result={}
    for powers in monomials(n):
        i,j,k,d,l=powers
        radial=centered_exponential(j+d)
        if not radial:continue
        ordinary=sum(c*gaussian(power+j+l) for power,c in hermite_product(i,k).items())
        if not ordinary:continue
        coefficient=factorial(n)*radial*ordinary
        coefficient//=factorial(i)*factorial(j)*factorial(k)*factorial(d)*factorial(l)
        result[powers]=coefficient
    content=0
    for coefficient in result.values():content=gcd(content,coefficient)
    assert content>0
    return {p:c//content for p,c in result.items()},content

def rows_in_degree(degree):
    return [(n,multiplier) for n in range(2,min(10,degree)+1)
            for multiplier in monomials(degree-n)]

def row_entries(label,column_lookup):
    n,multiplier=label
    result={}
    for powers,coefficient in moment(n)[0].items():
        shifted=tuple(a+b for a,b in zip(powers,multiplier))
        column=column_lookup.get(shifted)
        if column is not None:result[column]=coefficient%PRIME
    return result

def literal_matrix(degree,selected_rows=None,selected_columns=None):
    columns=monomials(degree)
    labels=rows_in_degree(degree)
    if selected_rows is None:selected_rows=range(len(labels))
    if selected_columns is None:selected_columns=range(len(columns))
    selected_rows=list(selected_rows)
    selected_columns=list(selected_columns)
    lookup={columns[index]:j for j,index in enumerate(selected_columns)}
    width=len(selected_columns)
    flat=[]
    for index in selected_rows:
        row=[0]*width
        for column,coefficient in row_entries(labels[index],lookup).items():row[column]=coefficient
        flat.extend(row)
    return nmod_mat(len(selected_rows),width,flat,PRIME)

def verify_record(record):
    degree=record['degree']
    assert 2<=degree<=11
    columns=monomials(degree)
    labels=rows_in_degree(degree)
    rank=len(columns)-HILBERT[degree]
    rows,cols=record['selected_rows'],record['selected_columns']
    assert record['ambient_monomials']==len(columns)
    assert record['macaulay_rows']==len(labels)
    assert record['rank_lower_bound']==rank
    assert len(rows)==len(set(rows))==rank
    assert len(cols)==len(set(cols))==rank
    assert all(isinstance(i,int) and 0<=i<len(labels) for i in rows)
    assert all(isinstance(i,int) and 0<=i<len(columns) for i in cols)
    determinant=int(literal_matrix(degree,rows,cols).det())
    assert determinant!=0
    assert determinant==record['nonzero_determinant_mod_prime']
    return dict(degree=degree,minor_size=rank,determinant_mod_prime=determinant,
                quotient_dimension_upper_bound=HILBERT[degree])

def generators_description():
    return {str(n):dict(content=moment(n)[1],terms=[[list(p),c] for p,c in sorted(moment(n)[0].items())])
            for n in range(2,11)}


def multiply_spatial(a,b):
    result={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():
            key=(i+k,j+l)
            result[key]=result.get(key,0)+c*d
    return {key:value for key,value in result.items() if value}

def original_source_controls():
    fixtures=[(1,2,-1,3,2),(-1,0,2,-2,1),(0,1,-2,-1,3)]
    checked=0
    for a,b,c,d,l in fixtures:
        phase={(0,3):a,(1,1):b,(0,2):c,(1,0):d,
               (0,1):l-3*a-b,(0,0):-c-d}
        phase={key:value for key,value in phase.items() if value}
        power={(0,0):1}
        for n in range(1,8):
            power=multiply_spatial(power,phase)
            original=sum(coefficient*factorial(i)*gaussian(j)
                         for (i,j),coefficient in power.items())
            if n==1:assert original==0
            else:
                polynomial,content=moment(n)
                expected=content*sum(coefficient*a**p[0]*b**p[1]*c**p[2]*d**p[3]*l**p[4]
                                     for p,coefficient in polynomial.items())
                assert original==expected,(fixtures,n)
            checked+=1
    return dict(fixtures=len(fixtures),orders=[1,7],exact_comparisons=checked)

def matrix_for_labels(degree,labels,selected_rows,selected_columns):
    columns=monomials(degree)
    lookup={columns[index]:j for j,index in enumerate(selected_columns)}
    result=nmod_mat(len(selected_rows),len(selected_columns),PRIME)
    for i,index in enumerate(selected_rows):
        for j,c in row_entries(labels[index],lookup).items():result[i,j]=c
    return result

def unscaled_data(degree,labels,record,primitive_determinant):
    counts=Counter(labels[i][0] for i in record['selected_rows'])
    scale=1
    for n,count in counts.items():scale=scale*pow(moment(n)[1],count,PRIME)%PRIME
    determinant=primitive_determinant*scale%PRIME
    assert determinant
    return dict(degree=degree,minor_size=len(record['selected_rows']),
                moment_row_counts=dict(sorted(counts.items())),
                universal_determinant_degree=sum(n*c for n,c in counts.items()),
                unscaled_determinant_mod_prime=determinant)

def verify_ci(record):
    assert record['schema']=='generic-centered-five-direction-CI-minor-v1'
    assert record['prefix']==6 and record['degree']==16
    assert record['prime']==PRIME and record['coefficient_field']=='QQ'
    assert record['variables']==['a','b','c','d','l']
    degree=16
    columns=monomials(degree)
    labels=[(n,u) for n in range(2,7) for u in monomials(degree-n)]
    selected_rows=record['selected_rows']
    selected_columns=record['selected_columns']
    assert len(columns)==4845 and len(labels)==9626
    assert record['ambient_monomials']==len(columns)
    assert record['macaulay_rows']==len(labels)
    assert record['rank_lower_bound']==len(columns)
    assert len(selected_rows)==len(set(selected_rows))==len(columns)
    assert all(type(i) is int and 0<=i<len(labels) for i in selected_rows)
    assert selected_columns==list(range(len(columns)))
    description={str(n):dict(content=moment(n)[1],terms=[[list(p),c] for p,c in sorted(moment(n)[0].items())])
                 for n in range(2,7)}
    digest=sha256(json.dumps(description,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert record['generators_sha256']==digest
    determinant=int(matrix_for_labels(degree,labels,selected_rows,selected_columns).det())
    assert determinant==record['nonzero_determinant_mod_prime']==6771
    return unscaled_data(degree,labels,record,determinant)

def hilbert_controls():
    profile=[comb(r+4,4) for r in range(17)]
    for degree in range(2,7):
        profile=[value-(profile[r-degree] if r>=degree else 0)
                 for r,value in enumerate(profile)]
    expected=[1,5,14,29,49,71,90,101,101,90,71,49,29,14,5,1,0]
    assert profile==expected and sum(profile)==720
    lower=[profile[r]-sum(profile[r-j] for j in range(7,min(10,r)+1))
           for r in range(11)]
    assert lower==HILBERT[:11] and sum(lower)==546
    prefix9_lower=profile[10]-sum(profile[10-j] for j in range(7,10))
    assert prefix9_lower==23 and lower[10]==22
    return dict(complete_intersection_degrees=[2,3,4,5,6],CI_Hilbert=profile,
                CI_length=720,full_Hilbert=HILBERT,full_length=546,
                prefix9_degree10_lower_bound=prefix9_lower,
                prefix10_degree10_dimension=22)

def centering_controls():
    directions=[{(0,3):1,(0,1):-3},{(1,1):1,(0,1):-1},
                {(0,2):1,(0,0):-1},{(1,0):1,(0,0):-1},{(0,1):1}]
    for phase in directions:
        assert sum(c*factorial(i)*gaussian(j) for (i,j),c in phase.items())==0
        assert max(2*i+j for i,j in phase)<=3
    return dict(directions=5,spatial_degree_max=3,all_original_Gaussian_means_zero=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve()
    parser.add_argument('--axial-minors',type=Path,default=here.with_name('gaussian_axial_cubic_macaulay_minors.json'))
    parser.add_argument('--ci-minor',type=Path,default=here.with_name('gaussian_generic_five_direction_ci_minor.json'))
    parser.add_argument('--proof',type=Path,default=here.with_name('GAUSSIAN_GENERIC_FIVE_DIRECTION_MOMENT_SCHEMES.md'))
    parser.add_argument('--primary-manifest',type=Path,default=here.parent/'sources/generic-five-direction/manifest.json')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    started=time.monotonic()
    assert all(PRIME%d for d in range(2,isqrt(PRIME)+1))
    axial=json.loads(args.axial_minors.read_bytes())
    assert axial['coefficient_field']=='QQ' and axial['prime']==PRIME
    description=generators_description()
    digest=sha256(json.dumps(description,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert axial['generators_sha256']==digest
    assert axial['variables']==['a','b','c','d','l']
    assert axial['target_hilbert']==HILBERT
    assert [v['degree'] for v in axial['records']]==list(range(2,12))
    primitive=[verify_record(v) for v in axial['records']]
    unscaled=[unscaled_data(v['degree'],rows_in_degree(v['degree']),v,w['determinant_mod_prime'])
              for v,w in zip(axial['records'],primitive)]
    assert max(unscaled[-1]['moment_row_counts'])==8
    assert unscaled[-1]['universal_determinant_degree']==4133
    assert unscaled[-1]['unscaled_determinant_mod_prime']==5815
    ci=json.loads(args.ci_minor.read_bytes())
    ci_result=verify_ci(ci)
    corrupt=dict(ci)
    corrupt['selected_rows']=list(corrupt['selected_rows'])
    corrupt['selected_rows'][0]=corrupt['selected_rows'][1]
    try:verify_ci(corrupt)
    except AssertionError:pass
    else:raise AssertionError('A duplicate-row CI minor was accepted')
    corrupt=dict(ci)
    corrupt['prefix']=5
    try:verify_ci(corrupt)
    except AssertionError:pass
    else:raise AssertionError('Wrong-prefix CI metadata was accepted')
    corrupt=dict(axial['records'][-1])
    corrupt['selected_rows']=list(corrupt['selected_rows'])
    corrupt['selected_rows'][0]=corrupt['selected_rows'][1]
    try:verify_record(corrupt)
    except AssertionError:pass
    else:raise AssertionError('A duplicate-row moment minor was accepted')
    manifest=json.loads(args.primary_manifest.read_bytes())
    assert {v['tag'] for v in manifest['sources']}=={'00NQ','02JN'}
    primary=[]
    for record in manifest['sources']:
        path=args.primary_manifest.parent/record['file']
        assert sha256(path.read_bytes()).hexdigest()==record['sha256']
        primary.append(path)
    files=[here,args.axial_minors,args.ci_minor,args.proof,args.primary_manifest]+primary
    result=dict(status='passed',coefficient_field='QQ',prime=PRIME,
                checked_on=datetime.now(timezone.utc).date().isoformat(),
                python=platform.python_version(),python_flint=flint.__version__,
                original_centering_controls=centering_controls(),
                original_moment_controls=original_source_controls(),
                primitive_axial_minors=primitive,unscaled_universal_minors=unscaled,
                complete_intersection_minor=ci_result,Hilbert_dimension_controls=hilbert_controls(),
                universal_ideal_containments=['Delta11*m^11 subset I8','epsilon*m^16 subset I6'],
                generic_ordinary_prefix=10,generic_geometric_prefix=6,
                generic_sharp_coefficient_nilpotence=11,generic_moment_length=546,
                generic_original_marked_tail_under_first8=11,
                rejection_controls=['duplicate_CI_minor_row','wrong_CI_prefix_metadata','duplicate_degree11_minor_row'],
                primary_manifest_consistency=dict(status='passed',sources=len(primary)),
                sha256={p.resolve().relative_to(here.parents[3]).as_posix():sha256(p.read_bytes()).hexdigest() for p in files},
                elapsed_seconds=round(time.monotonic()-started,3),
                scope='Literal nonzero integer minors and finite Hilbert arithmetic; unrestricted generic/family/field claims use the pinned adjugate and complete-intersection proof. No analytic comparison, Groebner generation, RREF, reconstruction or runtime scratch dependency.',
                independent_replay=False,formal_verification=False,external_review=False)
    output=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(output,encoding='utf-8')
    print(output,end='')

if __name__=='__main__':main()
