#!/usr/bin/env python3
"""Exact terminal-nine retained rows above the reciprocal cap 1/8."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv
from verify_multivariate_terminal_nine_private_linear import paired_box
from verify_multivariate_terminal_seven_paired_light import radial as paired_radial
from verify_multivariate_terminal_six_intermediate_cayley import rational

from flint import acb_poly, arb, ctx

from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_retained_rows import core_powers, ball, RAY_FRACTIONS
from verify_multivariate_terminal_six_private_angles import multiply_interval, require

from verify_multivariate_terminal_nine_small_outer import (
    exponent as anchor_exponent, exponents, verify_convexity,
)

ORDERS=(8,9)
CURVATURE_CAP = Q(53, 10)
CAP=Q(1,8)
POSITIONS=tuple(dict.fromkeys((*RAY_FRACTIONS,*(Q(k,20) for k in range(1,20)))))


def rows_above(cap,budget):
    for d in range(1,1+int(budget*(1+1/(8*cap)))):
        for a in range(1,1+int(Q(d,8)/cap-Q(1,8))):
            s=Q(d*d,a+d)
            if s<=budget and Q(d,8*a+1)>cap:
                yield a,d,s


def cores():
    core_set={}
    for a,d,s in rows_above(CAP,CURVATURE_CAP):
        rest=CURVATURE_CAP-s
        candidates=tuple((aa,dd,o,ss) for aa,dd,ss in rows_above(CAP,rest) for o in (-1,1))
        def collect(i,budget,selected):
            if i==len(candidates):
                key=min(tuple(sorted(selected)),tuple(sorted((aa,dd,-o) for aa,dd,o in selected)))
                core_set[key]=budget
                return
            aa,dd,o,ss=candidates[i]
            for n in range(1+int(budget/ss)):
                collect(i+1,budget-n*ss,selected+((aa,dd,o),)*n)
        collect(0,rest,((a,d,1),))
    return sorted(core_set.items())


@lru_cache(maxsize=50000)
def power_bounds(v,p):
    if v in (0,1):
        return Q(v),Q(v)
    z=ball(v)**ball(p)*10**8
    return Q(int(z.lower().floor().fmpz()),10**8),Q(int(z.upper().ceil().fmpz()),10**8)


def coefficient_box(j,u,shape):
    return [(Q(1),Q(1))]+[(comb(j,k)*power_bounds(u[0],hi)[0],
                            comb(j,k)*power_bounds(u[1],lo)[1])
                           for k,(lo,hi) in enumerate(shape,1)]+[(Q(1),Q(1))]


def boxes(core_bounds,tail_bounds,tail,core_attenuation):
    result={j:[multiply_interval(a,b) for a,b in zip(core_bounds[j],coefficient_box(j,tail,tail_bounds[j]),strict=True)]
            for j in ORDERS}
    for k in range(1,9):
        upper=[]
        for distance,index in ((k,k),(9-k,k-1)):
            factor=core_attenuation[distance]*power_bounds(tail[1],Q(distance**2,96))[1]
            upper.append(Q(comb(9,k),comb(8,index))*result[8][index][1]*factor)
        lo,hi=result[9][k]
        result[9][k]=lo,min(hi,*upper)
    return result


def nominal_cosines(core,midpoint,shape):
    answers=[]
    for j in ORDERS:
        coeff=[arb(1)]
        for k in range(1,j):
            f=anchor_exponent(j,k)
            require(shape[j][k-1][0]<=f<=shape[j][k-1][1],("Gaussian anchor",j,k,shape[j][k-1]))
            coeff.append(comb(j,k)*ball(core[j][k])**ball(Q(1,j))*ball(midpoint)**ball(f))
        coeff.append(arb(1))
        try:
            roots=acb_poly(coeff).roots(maxprec=256)
        except ValueError:
            return None
        upper=[z.real/abs(z) for z in roots if z.imag>0]
        if len(upper)!=j//2 or sum(z.imag<0 for z in roots)!=j//2:
            return None
        answers.extend((cosine,j) for cosine in upper)
    answers.sort(key=lambda item:float(item[0].mid()),reverse=True)
    if tuple(j for _,j in answers)!=(9,8)*4:
        return None
    if any(not a>b for (a,_),(b,_) in zip(answers,answers[1:])):
        return None
    return [cosine for cosine,_ in answers]


def rays(angles,coefficient_boxes,models=None):
    positions = POSITIONS[:12] if models else POSITIONS[:5]
    count=0
    for left,right in zip(angles,angles[1:]):
        for p in positions:
            proposal=float(left.mid())*float(p)+float(right.mid())*float(1-p)
            c=Q(round(100000*proposal),100000)
            if not left>ball(c)>right:
                continue
            values=[]
            for j in ORDERS:
                value=radial_replay(j,coefficient_boxes[j],(c,c))
                if value is None and models and models[j] is not None:
                    value=paired_radial(j,rational(c),models[j])
                if value is None:
                    break
                values.append(value)
            if len(values)==2:
                count+=sum(values)
                break
        else:
            return None
    return count

def verify_core(rows,budget):
    powers={j:core_powers(j,rows) for j in ORDERS}
    core_bounds={j:[power_bounds(value,Q(1,j)) for value in powers[j]] for j in ORDERS}
    for j in ORDERS:
        cap=[comb(j,k)*hi for k,(_,hi) in enumerate(core_bounds[j])]
        require(sum(a*INNER**k for k,a in enumerate(cap) if k)<1,("inner",rows,j))
        require(sum(a*OUTER**k for k,a in enumerate(cap[:-1]))<OUTER**j,("outer",rows,j))

    full_cores={j:tuple((comb(j,k)*lo,comb(j,k)*hi)
                       for k,(lo,hi) in enumerate(core_bounds[j])) for j in ORDERS}
    lo=power_bounds(Q(8,15),budget)[0]
    hi=Q(1)
    curvature=CURVATURE_CAP-budget
    attenuation={k:power_bounds(Q(152,153),curvature*k*k)[1] for k in range(1,9)}
    pending=[((lo,hi),0)]
    leaves=cells=0
    while pending:
        tail,depth=pending.pop()
        midpoint=sum(tail)/2
        effective=((tail[0]**-2-1)/7)*2000
        t8=min(CAP,Q(max(1,-(-effective.numerator//effective.denominator)),2000))
        shape=exponents(t8)
        angles=nominal_cosines(powers,midpoint,shape)
        if angles is not None:

            values=boxes(core_bounds,shape,tail,attenuation)
            count=rays(angles,values)
            if count is None and depth>=4:
                models={j:paired_box(j,full_cores[j],tail,(t8,8*t8/(9-t8))) for j in ORDERS}
                count=rays(angles,values,models)
            if count is not None:
                leaves+=1
                cells+=count
                continue
        require(depth<11,("unresolved tail",rows,tuple(map(float,tail)),angles))
        pending.extend((((tail[0],midpoint), depth+1),((midpoint,tail[1]),depth+1)))
    return leaves,cells


def main():
    ctx.prec=128
    verify_convexity()
    require(Q(153,152)**96>Q(15,8),"tail attenuation from outer loss")
    retained = cores()
    require(len(retained) == 93, "retained list count")
    shard, shards = tuple(map(int, argv[1:])) if len(argv) > 1 else (0, 1)
    require(0 <= shard < shards, "shard index")
    count = leaves = cells = 0
    for rows, budget in retained[shard::shards]:
        row_leaves, row_cells = verify_core(rows, budget)
        count += 1
        leaves += row_leaves
        cells += row_cells
    print(f"Verified {count} retained lists, {leaves} outer cells and {cells} radial cells.")


if __name__ == "__main__":
    main()
