"""Exact grouped disk separation on terminal-fourteen low outer means."""

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from math import prod
from sys import argv

from flint import acb, acb_poly, arb, ctx

from verify_multivariate_terminal_fourteen_retained_profiles import configure_shared_model
from verify_multivariate_terminal_fourteen_shared_disks import disk_model, propose_centers, radius_at
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds

BANDS = {
    4: (Q(1, 8), Q(1, 8), True, 53),
    5: (Q(1, 8), Q(1, 8), False, 52),
    6: (Q(1, 8), Q(1, 8), False, 105),
    7: (Q(1, 8), Q(1, 8), False, 181),
    8: (Q(1, 8), Q(1, 8), True, 337),
    9: (Q(1, 8), Q(1, 8), False, 1234),
    10: (Q(1, 9), Q(1, 9), False, 347),
    11: (Q(1, 8), Q(1, 9), False, 914),
    12: (Q(1, 10), Q(1, 10), False, 295),
    13: (Q(1, 10), Q(1, 10), True, 304),
    14: (Q(1, 10), Q(1, 10), False, 1477),
    15: (Q(1, 11), Q(5, 66), False, 106),
    16: (Q(3, 40), Q(3, 40), False, 163),
    17: (Q(3, 40), Q(3, 40), False, 303),
    18: (Q(3, 40), Q(3, 40), True, 339),
    19: (Q(3, 40), Q(3, 40), False, 637),
    20: (Q(3, 40), Q(3, 40), False, 1285),
    21: (Q(1, 14), Q(1, 14), False, 703),
    22: (Q(1, 14), Q(1, 14), False, 809),
    23: (Q(1, 14), Q(1, 14), False, 831),
    24: (Q(1, 14), Q(1, 14), False, 3878),
}

ctx.prec=112
RADII=tuple(arb(k)/2**power for power in range(20,3,-1) for k in (4,5,6,7))

N,UPPER=13,14
INTEGER_FLOOR=Q(1,80)

@lru_cache(maxsize=8192)
def simple_ratio(cap):
    return 1+log_bounds(1/(1-cap/UPPER))[1]/log_bounds(1+(N-1)*cap)[0]

@lru_cache(maxsize=8192)
def row_square(base,gap):
    return prod(Q(N*base+h,N*base+(N-1)*gap+h) for h in range(1,gap+1))

@lru_cache(maxsize=8192)
def row_ratio(base,gap):
    lower=sum(log_bounds(1+Q((N-1)*gap,N*base+h))[0] for h in range(1,gap+1))
    upper=sum(log_bounds(1+Q(N*gap,UPPER*base+h))[1] for h in range(1,gap+1))
    return upper/lower

@lru_cache(maxsize=8192)
def envelope(cap,u):
    floor=min(cap,max(INTEGER_FLOOR,cap/2))
    upper=simple_ratio(floor)
    reciprocal=floor
    if floor==cap:
        return reciprocal,upper
    u_square=u**2
    minimum_factor=1+(N-1)*floor/(1+floor)
    gap=1
    weight=minimum_factor*u_square
    while weight<=1:
        cutoff=(gap/cap-1)/N
        lo=max(1,-(-cutoff.numerator//cutoff.denominator))
        cutoff=(gap/floor-1)/N
        hi=-(-cutoff.numerator//cutoff.denominator)-1
        if lo<=hi and row_square(hi,gap)>=u_square:
            while lo<hi:
                mid=(lo+hi)//2
                if row_square(mid,gap)>=u_square:
                    hi=mid
                else:
                    lo=mid+1
            reciprocal=max(reciprocal,Q(gap,N*lo+1))
            upper=max(upper,row_ratio(lo,gap))
        gap+=1
        weight*=minimum_factor
    return reciprocal,min(upper,simple_ratio(reciprocal))

@lru_cache(maxsize=8192)
def exact_cap(cap,u):
    return envelope(cap,u)[0]

@lru_cache(maxsize=8192)
def ordinary_ratio_upper(cap,u):
    return envelope(cap,u)[1]

@dataclass
class RootBox:
    center: object
    radius: object
    drift: tuple
    error: object

def isolate(center,groups,order,radii):
    shift=acb_poly([center,1])
    translated=tuple((weight,p(shift)) for weight,p in groups)
    coeff=sum((weight*p for weight,p in translated),acb_poly())
    sup=tuple(abs(coeff[k]).upper() for k in range(order+1))
    inf=abs(coeff[1]).lower()
    for r in radii:
        high=sum(r**k*sup[k] for k in range(2,order+1))
        if r*inf>sup[0]+high:
            slope=acb(coeff[1].real.mid(),coeff[1].imag.mid())
            if not abs(slope).lower()>0:
                return None
            error=((r*abs(coeff[1]-slope).upper()+high)/abs(slope).lower()).upper()
            return RootBox(center,r,tuple((weight,-p[0]/slope) for weight,p in translated),error)
    return None

def projected(box,axis):
    return sum(weight*(axis.conjugate()*value).imag for weight,value in box.drift)

def ordered(a,b):
    c,d=a.center,b.center
    angle=(c.conjugate()*d).imag
    if angle>abs(c)*b.radius+abs(d)*a.radius+a.radius*b.radius:
        return True
    margin=angle+projected(b,c)-projected(a,d)
    return margin>abs(c)*b.error+abs(d)*a.error+a.radius*b.radius

def check_rectangle(model,means,skew,cap,radii):
    models={}
    for order in model.ORDERS:
        ceiling=model.ordinary_ratio_upper(cap,means[0]) if order==model.UPPER else None
        raw=model.coefficient_model(order,cap,means,skew,True,precise=True,monotone=True,
                 coupled=True,ratio_ceiling=ceiling,loss_transport=True,
                 center_ratio=True,center_skew=True)
        if raw is None:
            return True
        if not raw:
            return False
        models[order]=disk_model(raw)
    centers=propose_centers({order: data[0] for order,data in models.items()})
    if centers is None:
        return False
    upper,(odd,negative)=centers
    boxes=[]
    for i,c in upper:
        item=isolate(c,models[i][1],i,radii)
        if item is None or not c.imag>item.radius:
            return False
        boxes.append(item)
    if not all(ordered(a,b) for a,b in zip(boxes,boxes[1:])):
        return False
    r=radius_at(negative,models[odd][1],odd,radii)
    return r is not None and negative.real < -r

def check_band(model, lo, hi, cap):
    cap = model.exact_cap(cap, lo)
    pending = [((lo, hi), (Q(0), Q(1)), 0, 0)]
    accepted = 0
    while pending:
        means, signed, mean_depth, skew_depth = pending.pop()
        if check_rectangle(model, means, signed, cap, RADII):
            accepted += 1
            continue
        if mean_depth + skew_depth >= 20:
            raise ArithmeticError(("unresolved rectangle", means, signed, cap))
        if mean_depth < skew_depth + 2:
            mid = sum(means) / 2
            pending.extend((
                ((means[0], mid), signed, mean_depth + 1, skew_depth),
                ((mid, means[1]), signed, mean_depth + 1, skew_depth),
            ))
        else:
            mid = sum(signed) / 2
            pending.extend((
                (means, (signed[0], mid), mean_depth, skew_depth + 1),
                (means, (mid, signed[1]), mean_depth, skew_depth + 1),
            ))
    return accepted



def main():
    selected={int(value) for value in argv[1:]}
    if selected-BANDS.keys():
        raise ValueError("Select slabs by their integer lower endpoint.")
    model=configure_shared_model()
    model.FLOOR=INTEGER_FLOOR
    original_cap,original_ratio=model.exact_cap,model.ordinary_ratio_upper
    total=0
    for lo,(declared,effective,original,expected) in BANDS.items():
        if selected and lo not in selected:
            continue
        model.exact_cap=original_cap if original else exact_cap
        model.ordinary_ratio_upper=original_ratio if original else ordinary_ratio_upper
        lo_u,hi_u=Q(lo,100),Q(lo+1,100)
        if model.exact_cap(declared,lo_u)!=effective:
            raise ArithmeticError(("reciprocal cap",lo,effective))
        checked=check_band(model,lo_u,hi_u,declared)
        if checked!=expected:
            raise ArithmeticError(("unexpected rectangle count",lo,checked,expected))
        total+=checked
        print(lo,lo+1,declared,effective,checked,flush=True)
        model.coefficient_model.cache_clear()
        model.ordinary_ratio_upper.cache_clear()
    print(f"Exact accepted rectangles: {total}.",flush=True)


if __name__=="__main__":
    main()
