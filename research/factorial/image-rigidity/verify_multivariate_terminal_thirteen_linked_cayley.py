"""Exact matched-skew Cayley rays for the terminal-thirteen seam."""

from sys import argv
from fractions import Fraction as Q

from verify_multivariate_terminal_thirteen_cayley import (
    GRID, I, ONE, REFERENCE_RATIO as R, accurate_projection, anchor_for, bisect,
    check_ray, height, interval, limits_for, log_bounds, moment_ratio, rational,
    ray_for, require, skew_model, strength_factor, validate_rays,
)

CAP=Q(1,73)
DISK=Q(13,5)
NODES=tuple(range(66,71))
SCALE=100000
RAYS = {
    66: (122510, 109688, 95343, 84255, 72225, 62424, 50988, 42391, 32038, 23346, 13664),
    67: (122592, 109752, 95388, 84290, 72250, 62443, 51001, 42401, 32044, 23351, 13667),
    68: (122674, 109817, 95433, 84325, 72274, 62462, 51014, 42411, 32051, 23355, 13669),
    69: (122756, 109881, 95478, 84360, 72299, 62481, 51027, 42421, 32057, 23360, 13672),
    70: (122839, 109946, 95524, 84396, 72324, 62500, 51040, 42431, 32064, 23364, 13674),
}

# The movement term changes the first ray only.
SKEW_RAY=Q(196,10000)

def transport_error(limits,left,right):
    c=limits.cap
    B=143*c/(24+143*c)
    E=11*(1+13*c)/(24+143*c)
    L=1+log_bounds(1/(1-c/13))[1]/log_bounds(1+11*c)[0]
    U=log_bounds(1/(1-Q(right,1000)))[1]
    q_min=Q(left,1000)*limits.ratio[0]
    value=R**2*Q(52,11)*B*E**2*(L-R)*U/q_min**2
    return Q(-(-value.numerator*10**7//value.denominator),10**7)

def linked_first(order,left,right,limits,rays):
    ratio=interval(*limits.ratio)
    lam=interval(Q(left,1000),Q(right,1000)).sqrt()
    C=strength_factor(limits.qmax)
    lower_radius=C*moment_ratio(12,left,right,limits)/2
    A=moment_ratio(order,left,right,limits)
    upper_radius=C*A/2
    error=transport_error(limits,left,right)
    pending=[(lam,interval(Q(0),DISK),ratio,lower_radius*Q(k,8),lower_radius*Q(k+1,8),0)
             for k in range(-8,8)]
    count=0
    while pending:
        lam,x,r,sklo,skhi,depth=pending.pop()
        skew=interval(sklo,skhi)
        p=lam.square()
        sr=r.sqrt() if order==13 else ONE
        scale=lam*sr
        q=scale.square()
        b=(ray_for(p,0,left,right,rays)+rational(SKEW_RAY)*skew.square())/sr
        if order==13:
            ps=skew*(rational(R)/r).square()+interval(-error,error)
            blo=max(-upper_radius,Q(ps.lo,GRID))
            bhi=min(upper_radius,Q(ps.hi,GRID))
        else:
            blo,bhi=sklo,skhi
        if blo>bhi:
            continue
        base,groups,box,ok=skew_model(order,scale.lo,scale.hi,A,limits,True,blo,bhi)
        if not ok:
            require(depth<32,("scaled disk",order,left,right,lam.shown(),r.shown(),sklo,skhi))
            pending.extend((half,x,r,sklo,skhi,depth+1) for half in bisect(lam))
            continue
        y=height(q,x,b)
        if accurate_projection(order,x,y,base,groups,box):
            count+=1
            continue
        require(depth<46,("linked ray",order,left,right,lam.shown(),x.shown(),r.shown(),
                          float(sklo),float(skhi),count))
        wl,wx,wr=10*(lam.hi-lam.lo),x.hi-x.lo,10*(r.hi-r.lo) if order==13 else 0
        ws=Q(GRID,55)*(skhi-sklo)
        if ws>max(wl,wx,wr) and x.hi<GRID//3:
            center=(sklo+skhi)/2
            pending.extend(((lam,x,r,sklo,center,depth+1),(lam,x,r,center,skhi,depth+1)))
        elif wr>max(wl,wx):
            pending.extend((lam,x,half,sklo,skhi,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,sklo,skhi,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,sklo,skhi,depth+1) for half in bisect(x))
    return count

EXPECTED_CELLS = {
    (12, 66): (1740, 45, 351, 48, 162, 44, 183, 36, 37, 28, 22),
    (13, 66): (1851, 57, 672, 83, 346, 79, 72, 44, 49, 32, 30),
    (12, 67): (2318, 46, 476, 53, 273, 45, 280, 37, 39, 31, 24),
    (13, 67): (26435, 65, 1029, 87, 527, 82, 240, 52, 350, 35, 55),
    (12, 68): (2330, 46, 476, 53, 266, 45, 280, 37, 39, 31, 24),
    (13, 68): (26961, 65, 1017, 87, 527, 82, 240, 52, 365, 35, 55),
    (12, 69): (2295, 47, 465, 53, 241, 45, 279, 37, 39, 31, 23),
    (13, 69): (19227, 63, 955, 87, 493, 83, 204, 50, 310, 35, 54),
}


def run(order,index,left,right):
    require(tuple(RAYS)==NODES and left in NODES and right in NODES and left<right,"slab range")
    require(order in (0,12,13) and -1<=index<11,"order or ray")
    anchor_for(limits_for(CAP,NODES[0],NODES[0],DISK,False,refined=True),RAYS,NODES[0])
    count=0
    for slab in range(left,right):
        lim=limits_for(CAP,slab,slab+1,DISK,False,refined=True)
        validate_rays(lim,RAYS,slab,slab+1)
        C=strength_factor(lim.qmax)
        radius=C*moment_ratio(12,slab,slab+1,lim)/2
        max_level=Q(max(RAYS[slab][0],RAYS[slab+1][0]),SCALE)+SKEW_RAY*radius**2
        require(0<max_level<DISK and max_level**2*lim.qmax<1,"linked ray half-plane")
        for j in (12,13):
            if order and order!=j:
                continue
            for k in range(11):
                if index>=0 and index!=k:
                    continue
                if k==0:
                    cells=linked_first(j,slab,slab+1,lim,RAYS)
                else:
                    cells=check_ray(j,k,slab,slab+1,lim,RAYS,True,conditioned=True)
                require(cells == EXPECTED_CELLS[j,slab][k],
                        ("replay cell count", j, k, slab, cells))
                count+=cells
                print(j,k,slab,slab+1,cells,flush=True)
    print(f"Verified {count} projection cells.",flush=True)

if __name__=="__main__":
    run(*(tuple(map(int,argv[1:])) if len(argv)>1 else (0,-1,NODES[0],NODES[-1])))
