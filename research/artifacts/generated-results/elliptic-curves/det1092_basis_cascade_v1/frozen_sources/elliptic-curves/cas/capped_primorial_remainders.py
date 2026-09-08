"""Exact remainder trees whose products never need to grow beyond the modulus.

For positive a,b, min(min(a,C)*min(b,C),C)=min(a*b,C). Thus each
node stores its true product capped at C=M+1. Every propagated remainder
is at most M: reducing it modulo a product >=C leaves it unchanged, so
the cap preserves every leaf remainder exactly. The bit-length shortcut
only avoids a multiplication when its product is certainly larger than C.
"""

def residues(values, modulus):
    if not values or modulus<0 or any(v<=0 for v in values):
        raise ValueError('positive values and a nonnegative modulus are required')
    cap=modulus+1
    levels=[[min(v,cap) for v in values]]
    threshold=cap.bit_length()+2
    while len(levels[-1])>1:
        old=levels[-1];new=[]
        for i in range(0,len(old),2):
            if i+1==len(old):
                new.append(old[i]);continue
            a,b=old[i],old[i+1]
            if a==cap or b==cap or a.bit_length()+b.bit_length()>=threshold:
                new.append(cap)
            else:new.append(min(a*b,cap))
        levels.append(new)
    current=[modulus%levels[-1][0]]
    for level in reversed(levels[:-1]):
        current=[current[i//2]%v for i,v in enumerate(level)]
    return current
