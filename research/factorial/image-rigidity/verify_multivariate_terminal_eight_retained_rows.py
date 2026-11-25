"""Replay the terminal-eight retained lists with a light shared remainder."""

from fractions import Fraction as Q
from math import comb
from flint import acb_poly, arb, ctx
from verify_multivariate_terminal_six_retained_rows import core_powers, ball, RAY_FRACTIONS
from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_private_angles import multiply_interval, require
from verify_multivariate_terminal_eight_shared_high_curvature import exponent_bounds
from functools import lru_cache

ORDERS=(7,8)
MAX_CURVATURE=Q(19,4)
TAIL_CAP=Q(1,8)
EXCEPTIONS = {((1, 2, -1),), ((2, 2, -1),), ((3, 3, -1),)}

def rows_above(cap,budget):
    for gap in range(1, 1+int(budget*(1+1/(7*cap)))):
        for lower in range(1, 1+int(Q(gap,7)/cap-Q(1,7))):
            s=Q(gap*gap,lower+gap)
            if s<=budget and Q(gap,7*lower+1)>cap:
                yield lower,gap,s

def cores():
    core_set={}
    for a,d,s in rows_above(TAIL_CAP,MAX_CURVATURE):
        if s>=MAX_CURVATURE:
            continue
        rest=MAX_CURVATURE-s
        candidates=tuple((aa,dd,o,ss) for aa,dd,ss in rows_above(TAIL_CAP,rest) for o in (-1,1))
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
    for k in range(1,8):
        upper=[]
        for distance,index in ((k,k),(8-k,k-1)):
            factor=core_attenuation[distance]*power_bounds(tail[1],Q(distance**2,75))[1]
            upper.append(Q(comb(8,k),comb(7,index))*result[7][index][1]*factor)
        lo,hi=result[8][k]
        result[8][k]=lo,min(hi,*upper)
    return result

def nominal_cosines(core,midpoint,shape):
    answers=[]
    for j in ORDERS:
        coeff=[arb(1)]
        for k in range(1,j):
            f=Q(k*(j-k),j-1)*(1 if j==7 else Q(49,48))
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
    if tuple(j for _,j in answers)!=(8,7,8,7,8,7,8):
        return None
    if any(not a>b for (a,_),(b,_) in zip(answers,answers[1:])):
        return None
    return [cosine for cosine,_ in answers]

def rays(angles,coefficient_boxes):
    count=0
    for left,right in zip(angles,angles[1:]):
        for p in RAY_FRACTIONS:
            proposal=float(left.mid())*float(p)+float(right.mid())*float(1-p)
            c=Q(round(100000*proposal),100000)
            if not left>ball(c)>right:
                continue
            n=radial_replay(7,coefficient_boxes[7],(c,c))
            if n is None:
                continue
            m=radial_replay(8,coefficient_boxes[8],(c,c))
            if m is None:
                continue
            count+=n+m
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
    lo=power_bounds(Q(7,13),budget)[0]
    hi = Q(1)
    curvature=MAX_CURVATURE-budget
    attenuation={k:power_bounds(Q(119,120),curvature*k*k)[1] for k in range(1,8)}
    pending=[((lo,hi),0)]
    leaves=cells=0
    while pending:
        tail,depth=pending.pop()
        midpoint=sum(tail)/2
        effective=((tail[0]**-2-1)/6)*2000
        t7=min(TAIL_CAP,Q(max(1,-(-effective.numerator//effective.denominator)),2000))
        shape=exponent_bounds(t7,7*t7/(8-t7))
        angles=nominal_cosines(powers,midpoint,shape)
        if angles is not None:
            count=rays(angles,boxes(core_bounds,shape,tail,attenuation))
            if count is not None:
                leaves+=1
                cells+=count
                continue
        require(depth<15,("unresolved tail",rows,tuple(map(float,tail)),angles))
        pending.extend((((tail[0],midpoint), depth+1),((midpoint,tail[1]),depth+1)))
    return leaves,cells

def main():
    ctx.prec=128
    require(Q(120,119)**75>Q(13,7),"tail attenuation from outer loss")
    retained=cores()
    print("retained",len(retained),flush=True)
    keys = {rows for rows, _ in retained}
    require(len(keys) == 245 and EXCEPTIONS <= keys, "retained list enumeration")
    leaves = cells = count = 0
    for rows,b in retained:
        if rows in EXCEPTIONS:
            continue
        row_leaves, row_cells = verify_core(rows, b)
        count += 1
        leaves += row_leaves
        cells += row_cells
    print(f"certified {count} lists, {leaves} tail intervals, and {cells} radial cells")

if __name__=="__main__":
    main()
