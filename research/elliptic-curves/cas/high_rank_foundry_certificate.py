"""Select a smaller exact finite witness. This is not an independence prover."""
from copy import deepcopy


def compact_packet(packet):
    out=deepcopy(packet);rank=packet['rank_lower_bound'];pivots={};selected=[]
    for signature in sorted(packet['proof']['signatures'],key=lambda s:s['prime']):
        before=len(pivots)
        for row in signature['rows']:
            if len(row)!=rank or any(x not in (0,1) for x in row):
                raise ArithmeticError('invalid binary certificate row')
            v=sum(x<<i for i,x in enumerate(row))
            while v:
                bit=v.bit_length()-1
                if bit not in pivots:pivots[bit]=v;break
                v^=pivots[bit]
        if len(pivots)>before:selected.append(signature)
        if len(pivots)==rank:break
    if len(pivots)!=rank:raise ArithmeticError('finite certificate does not span the point columns')
    out['proof']['signatures']=selected
    out['status']='PASS_FOUNDRY_COMPACT_CERTIFICATE_PENDING_REPLAY'
    out['certificate_compaction']={'original_prime_count':len(packet['proof']['signatures']),
        'selected_prime_count':len(selected),'rule':'Keep rank-increasing good-prime row blocks in prime order. Both exact finite implementations must replay this subset before promotion.'}
    return out
