"""Exact grouped root disks on the terminal-fourteen lower mean bands."""
from fractions import Fraction as Q
from sys import argv

from flint import arb, ctx
import verify_multivariate_terminal_fourteen_shared_disks as disks
from verify_multivariate_terminal_fourteen_retained_profiles import configure_shared_model

BANDS = {50: 129, 51: 128, 52: 128, 53: 128, 54: 157, 55: 166, 56: 222, 57: 353, 58: 386, 59: 401, 60: 353, 61: 408, 62: 416, 63: 288}

def main():
    levels=set(map(int,argv[1:]))
    if levels-BANDS.keys():
        raise ValueError("Select slabs by the integer lower endpoint.")
    ctx.prec=112
    model=configure_shared_model()
    model.FLOOR=Q(1,80)
    radii=tuple(arb(k)/2**p for p in range(20,4,-1) for k in (4,5,6,7))
    total=0
    for level,expected in BANDS.items():
        if levels and level not in levels:
            continue
        cap=Q(1,40) if level<63 else Q(1,44)
        count=disks.check_band(model,Q(level,100),Q(level+1,100),cap,radii)
        if count!=expected:
            raise ArithmeticError(("unexpected disk count", level, count, expected))
        print(level,level+1,cap,count,flush=True)
        total+=count
        model.coefficient_model.cache_clear()
        model.ordinary_ratio_upper.cache_clear()
    print(f"Exact accepted rectangles: {total}.",flush=True)

if __name__=="__main__":
    main()
