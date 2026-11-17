"""Replay exact QQ/Fp certificates for the full axial cubic moment ideal.

The default replay uses literal polynomial division, S-pairs, nonzero
Macaulay minors, and a rational separator. It calls no Groebner generator,
ideal-membership oracle, Singular, RREF or rational reconstruction.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial, gcd, isqrt
from pathlib import Path
import time
import gzip
import platform
import copy

from flint import nmod_mat, fmpq, fmpq_mpoly, fmpq_mpoly_ctx
import flint

if not __debug__:
    raise RuntimeError('Exact assertions must be enabled')

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


FAMILY="a*(Y^3-3*Y)+b*(S-1)*Y+c*(Y^2-1)+d*(S-1)+l*Y"

def divides(smaller,larger):
    return all(a<=b for a,b in zip(smaller,larger))

def decode_basis(certificate):
    assert certificate['schema']=='axial-cubic-moment-ideal-certificate-v1'
    assert certificate['coefficient_field']=='QQ'
    assert certificate['parameters']==['a','b','c','d','l']
    assert certificate['phase']==FAMILY
    assert certificate['term_order']=='degrevlex'
    assert certificate['prefix']==10 and certificate['nilpotence_index']==11
    assert certificate['hilbert']==HILBERT
    assert len(certificate['groebner_basis'])==128
    context=fmpq_mpoly_ctx.get(('a','b','c','d','l'),'degrevlex')
    basis=[fmpq_mpoly(value,ctx=context) for value in certificate['groebner_basis']]
    assert all(basis)
    basis=[value/value.leading_coefficient() for value in basis]
    leading=[next(value.terms())[0] for value in basis]
    assert len(set(leading))==len(leading)
    for value,lm in zip(basis,leading):
        assert all(sum(powers)==sum(lm) for powers,coefficient in value.terms())
    return context,basis,leading

def divide_by_basis(value,basis):
    # Each operation is ordinary single-polynomial division. A zero final
    # remainder supplies an actual standard representation, without
    # presuming that this list is already a Groebner basis.
    passes=0
    while value:
        changed=False
        for divisor in basis:
            quotient,remainder=divmod(value,divisor)
            if quotient:
                value=remainder
                changed=True
                if not value:break
        if not changed:return value
        passes+=1
        if passes>10000:raise RuntimeError('Division guard exceeded')
    return value

def verify_basis(certificate):
    context,basis,leading=decode_basis(certificate)
    for n in range(2,11):
        assert not divide_by_basis(context.from_dict(moment(n)[0]),basis),n
    divided=coprime=0
    started=time.monotonic()
    for i,lmi in enumerate(leading):
        for j in range(i):
            lmj=leading[j]
            if all(min(a,b)==0 for a,b in zip(lmi,lmj)):
                coprime+=1
                continue
            lcm=tuple(max(a,b) for a,b in zip(lmi,lmj))
            mi=context.from_dict({tuple(a-b for a,b in zip(lcm,lmi)):1})
            mj=context.from_dict({tuple(a-b for a,b in zip(lcm,lmj)):1})
            assert not divide_by_basis(mi*basis[i]-mj*basis[j],basis),(i,j)
            divided+=1
        if i%16==15:
            print(json.dumps(dict(QQ_basis_rows_checked=i+1,elapsed_seconds=round(time.monotonic()-started,2))),flush=True)
    assert divided+coprime==128*127//2
    profile=[]
    standard10=None
    for degree in range(12):
        standard=[powers for powers in monomials(degree)
                  if not any(divides(lm,powers) for lm in leading)]
        profile.append(len(standard))
        if degree==10:
            assert (0,0,0,0,10) in standard
            standard10=[list(powers) for powers in standard]
    assert profile==HILBERT
    return dict(rows=len(basis),S_pairs_divided=divided,coprime_S_pairs=coprime,
                hilbert=profile,quotient_dimension=sum(profile),
                degree10_standard_monomials=standard10,
                nonzero_nilpotence_witness='l^10')

def decode_separator(certificate):
    separator=certificate['least_prefix_separator']
    assert separator['degree']==10 and separator['prefix']==9
    assert separator['target_moment']==10 and separator['target_pairing']=='1'
    result={}
    for term in separator['terms']:
        powers=tuple(term['powers'])
        assert len(powers)==5 and all(type(v) is int and v>=0 for v in powers)
        assert sum(powers)==10 and powers not in result
        coefficient=fmpq(term['coefficient'])
        assert coefficient
        result[powers]=coefficient
    assert len(result)==496
    return result

def separator_pairing(separator,polynomial,multiplier=ZERO):
    return sum(coefficient*separator.get(tuple(a+b for a,b in zip(powers,multiplier)),0)
               for powers,coefficient in polynomial.items())

def verify_separator(certificate):
    separator=decode_separator(certificate)
    assert separator_pairing(separator,moment(10)[0])==1
    checked=0
    for n in range(2,10):
        for multiplier in monomials(10-n):
            assert separator_pairing(separator,moment(n)[0],multiplier)==0,(n,multiplier)
            checked+=1
    return dict(nonzero_entries=len(separator),annihilated_degree10_products=checked,
                moment10_pairing='1',least_generating_prefix=10)

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

def rejection_controls(certificate,minors):
    bad=copy.deepcopy(certificate)
    bad['phase']='wrong original phase'
    try:decode_basis(bad)
    except AssertionError:pass
    else:raise AssertionError('Wrong original source was accepted')
    bad=copy.deepcopy(certificate)
    bad['groebner_basis'][0]+='+1'
    try:decode_basis(bad)
    except AssertionError:pass
    else:raise AssertionError('Inhomogeneous corrupt basis was accepted')
    bad=copy.deepcopy(certificate)
    separator=decode_separator(bad)
    target=next(p for p in separator if p in moment(10)[0])
    separator[target]+=1
    assert separator_pairing(separator,moment(10)[0])!=1
    corrupt=dict(minors['records'][-1])
    corrupt['selected_rows']=list(corrupt['selected_rows'])
    corrupt['selected_rows'][0]=corrupt['selected_rows'][1]
    try:verify_record(corrupt)
    except AssertionError:pass
    else:raise AssertionError('Duplicate-row false minor was accepted')
    return ['wrong_source','corrupt_basis','altered_dual','duplicate_minor_row']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    directory=Path(__file__).resolve().parent
    parser.add_argument('--certificate',type=Path,default=directory/'gaussian_axial_cubic_moment_ideal_certificate.json.gz')
    parser.add_argument('--minors',type=Path,default=directory/'gaussian_axial_cubic_macaulay_minors.json')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    started=time.monotonic()
    assert PRIME>1 and all(PRIME%d for d in range(2,isqrt(PRIME)+1))
    description=generators_description()
    generator_hash=sha256(json.dumps(description,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    with gzip.open(args.certificate,'rt',encoding='utf-8') as stream:certificate=json.load(stream)
    minors=json.loads(args.minors.read_bytes())
    assert certificate['generators_sha256']==generator_hash
    assert minors['coefficient_field']=='QQ' and minors['prime']==PRIME
    assert minors['variables']==['a','b','c','d','l']
    assert minors['generators_sha256']==generator_hash
    assert minors['target_hilbert']==HILBERT
    assert [record['degree'] for record in minors['records']]==list(range(2,12))
    controls=original_source_controls()
    rejections=rejection_controls(certificate,minors)
    basis=verify_basis(certificate)
    minor_results=[verify_record(record) for record in minors['records']]
    separator=verify_separator(certificate)
    pins={}
    imports=['GAUSSIAN_AXIAL_CUBIC_MOMENT_IDEAL.md',
             'GAUSSIAN_W_AFFINE_ARBITRARY_FIBRE_FINITE_IMAGES.md',
             'GAUSSIAN_PRIMITIVE_CLOSURE_UNDER_SUMS.md',
             'GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md',
             'MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md',
             '../../../unpublished/gaussian-boundary/polynomial_stationary_balls.tex',
             '../../../unpublished/gaussian-boundary/manuscript.tex']
    for name in imports:pins[name]=sha256((directory/name).read_bytes()).hexdigest()
    pins[Path(__file__).name]=sha256(Path(__file__).read_bytes()).hexdigest()
    pins[args.certificate.name]=sha256(args.certificate.read_bytes()).hexdigest()
    pins[args.minors.name]=sha256(args.minors.read_bytes()).hexdigest()
    result=dict(status='passed',coefficient_field='QQ',python=platform.python_version(),
                python_flint=flint.__version__,original_source_controls=controls,
                generators_sha256=generator_hash,groebner_superset=basis,
                nonzero_macaulay_minors=minor_results,least_prefix_separator=separator,
                proved_infinite_tail='moment degree n>=11 lies in m^11 subset I10',
                ideal_equality='I_infinity=I10 globally',nilpotence_index=11,
                quotient_length=546,least_generating_prefix=10,
                rejection_controls=rejections,sha256=pins,
                elapsed_seconds=round(time.monotonic()-started,3),
                scope='Complete finite certificate replay for the stated homogeneous axial slice; no ideal-prefix extrapolation, GB generation, RREF, or runtime scratch dependency.',
                independent_replay=False,formal_verification=False,external_review=False)
    output=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(output,encoding='utf-8',newline='\n')
    print(output,end='')

if __name__=='__main__':main()
