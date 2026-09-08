#!/usr/bin/env python3
"""Exact horizontal variants of the preserved v1 pointed sieve.

All finite searches use the same GMP modular worker. A variant is an explicit
PGL2(Q) substitution plus an exactly verified square ordinate scaling.
No exceptional points, public complements or rank labels enter construction.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt, gcd, lcm
from pathlib import Path
import json
import os
import subprocess
import time

import half_lattice_pointed_sieve as base


def normalize(values):
    values=tuple(map(Q,values))
    denominator=lcm(*(v.denominator for v in values))
    integers=tuple(int(v*denominator**2) for v in values)
    content=gcd(*integers)
    if not content:
        raise ArithmeticError('zero quartic')
    root=1
    for prime in base.SMALL_PRIMES:
        exponent=0
        while content%prime==0:
            content//=prime; exponent+=1
        root*=prime**(exponent//2)
    residual=isqrt(content)
    if residual**2==content:
        root*=residual
    scale=Q(denominator,root)
    return tuple(base.divide(v,root**2) for v in integers),scale


def multiply(m,n):
    a,b,c,d=m; e,f,g,h=n
    return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)


def inverse(m):
    a,b,c,d=m; det=a*d-b*c
    if not det:
        raise ValueError('singular horizontal change')
    return tuple(Q(x,det) for x in (d,-b,-c,a))


def reduction_matrix(f, seconds=15):
    polynomial='+'.join(f'({v})*x^{i}' for i,v in enumerate(f))
    program=f'C=hyperellred([{polynomial},0],&m); print("M|",m[2][1,1],"|",m[2][1,2],"|",m[2][2,1],"|",m[2][2,2]); quit\n'
    process=subprocess.run(['gp','-q','-s','256000000'],input=program,text=True,
                           capture_output=True,timeout=seconds,check=True)
    rows=[r for r in process.stdout.splitlines() if r.startswith('M|')]
    if len(rows)!=1 or '***' in process.stderr:
        raise ArithmeticError('PARI horizontal reduction failed: '+process.stderr[-500:])
    m=tuple(Q(x) for x in rows[0].split('|')[1:])
    if not m[0]*m[3]-m[1]*m[2]:
        raise ArithmeticError('PARI supplied singular horizontal matrix')
    return m


def horizontal(chart, specification):
    parts=specification.split(':')
    name=parts[0]
    if name=='gauss':
        m=(1,0,0,1)
    elif name=='red':
        m=reduction_matrix(chart.coefficients)
    elif name=='metric':
        factor=Q(parts[1])
        if factor<=0:
            raise ValueError('metric multiplier must be positive')
        den,k=chart.denominator,chart.shift
        aa=chart.model[3]*chart.curve_scale**4
        a=chart.base_point[0]*chart.curve_scale**2*den**2
        w=(abs(a)+den**2*(isqrt(abs(int(aa)))+1))*factor
        m=multiply(inverse(chart.matrix),base.gauss_matrix(den,k,w))
    elif name=='raw':
        den,k,u=chart.denominator,chart.shift,chart.curve_scale
        original=multiply((Q(den)/u,Q(k,den)/u,0,1),chart.matrix)
        m=inverse(original)
    else:
        raise ValueError('unknown horizontal specification: '+specification)
    # Rational slope boxes: z=(p*s+r)/q after the selected horizontal chart.
    # Positive p,q and arbitrary r are explicit, invertible PGL2 changes.
    if name in ('gauss','red') and len(parts)>1:
        p,q,r=map(int,parts[1].split(','))
        if p<=0 or q<=0:
            raise ValueError('invalid rational slope parametrization')
        m=multiply(m,(p,r,0,q))
    f,scale=normalize(base.binary_transform(chart.coefficients,m))
    if tuple(Q(v)*scale**2 for v in base.binary_transform(chart.coefficients,m))!=tuple(f):
        raise ArithmeticError('horizontal five-coefficient identity failed')
    return tuple(m),tuple(map(int,f)),Q(scale)


def sources():
    return {**base.provenance(),**{str(p.relative_to(base.ROOT)):sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__).resolve(),)}}


def run(*, model, points, representative, mask, specification, height, seconds):
    started=time.monotonic()
    centre=base.linear_combination(model,points,representative)
    if centre is None:
        raise ArithmeticError('nonzero centre at infinity')
    chart=base.make_chart(model,centre)
    m,f,scale=horizontal(chart,specification)
    record,hits=base.search_box(f,height,seconds)
    if f[4]>=0 and isqrt(f[4])**2==f[4]:
        hits+=((1,0,isqrt(f[4])),)
    found=set()
    a,b,c,d=m
    for n,den,root in hits:
        upper,lower=a*n+b*den,c*n+d*den
        for signed in {root,-root}:
            point=chart.map_point(upper,lower,Q(signed)/scale)
            if point is not None:
                if not base.point_on_short_curve(model,point):
                    raise ArithmeticError('horizontal point left E')
                found.add(point)
    record.update({'backend':'mw16_sensitivity_pointed_sieve_v1','mask':int(mask),
        'representative':list(map(int,representative)),
        'specification':specification,'base_point':{'x':str(centre[0]),'y':str(centre[1])},
        'pointed_chart':chart.record(),'horizontal_matrix':[str(x) for x in m],
        'ordinate_scale':str(scale),'coefficients':[str(v) for v in f],
        'primitive_square_hits':[list(map(str,p)) for p in hits],
        'infinity_checked':True,'wall_seconds':time.monotonic()-started,
        'hyperellminimalmodel_called':False,'hyperellred_called':specification.startswith('red'),
        'finite_curve_points':[{'x':str(x),'y':str(y)} for x,y in sorted(found)]})
    return record


def checkpoint(directory, **kwargs):
    key={'sources':sources(),'arguments':kwargs}
    canonical=json.dumps(key,sort_keys=True,default=str,separators=(',',':'))
    path=Path(directory)/(sha256(canonical.encode()).hexdigest()+'.json')
    if path.exists():
        saved=json.loads(path.read_text()); record=saved['record']
        if saved['input']!=json.loads(canonical) or saved['record_sha256']!=sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest():
            raise ArithmeticError('sensitivity checkpoint integrity failed')
        if record['status']=='bounded_search_complete':
            return record
    record=run(**kwargs)
    document={'input':json.loads(canonical),'record':record,
        'record_sha256':sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()}
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(document,sort_keys=True)+'\n'); temporary.replace(path)
    return record
