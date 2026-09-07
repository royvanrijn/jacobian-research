"""Frozen, target-free shortlist metrics. No artifact reads or orbit IDs."""
import hashlib
from fractions import Fraction

def coverage_key(fingerprint):
    return hashlib.sha256(('visibility-v3-coset:'+str(int(fingerprint))).encode()).hexdigest()

def tied_take(indices, metric, count, reverse=False):
    order=sorted(indices,key=metric,reverse=reverse)
    if len(order)<=count:return set(order)
    edge=metric(order[count-1])
    return {i for i in order if (metric(i)>=edge if reverse else metric(i)<=edge)}

def stratified(indices, metric, strata, count):
    order=sorted(indices,key=metric)
    selected=set()
    for j in range(strata):
        # Score-quantile centres. Exact integer distances, all boundary ties.
        centre=metric(order[min(len(order)-1,(2*j+1)*len(order)//(2*strata))])
        selected |= tied_take(indices,lambda i:abs(metric(i)-centre),count)
    return selected

def cheap_shortlist(norms, keys):
    ids=list(range(len(norms)))
    selected=tied_take(ids,lambda i:int(norms[i]),4,True)
    selected |= tied_take(ids,lambda i:int(norms[i]),4)
    selected |= stratified(ids,lambda i:int(norms[i]),16,2)
    selected |= tied_take(ids,lambda i:coverage_key(keys[i]),8)
    median=sorted(int(n) for n in norms)[len(norms)//2]
    poor=[i for i in ids if int(norms[i])<=median]
    selected |= tied_take(poor,lambda i:coverage_key(keys[i]),8)
    return sorted(selected,key=lambda i:(-int(norms[i]),int(keys[i])))

def chart_profile(mapping):
    values=[Fraction(str(x)) for x in mapping['discriminant_quartic']]
    sizes=[abs(x.numerator).bit_length()+x.denominator.bit_length() for x in values]
    return {'quartic_bits':sum(sizes),'quartic_max_bits':max(sizes)}

def final_shortlist(rows):
    ids=list(range(len(rows)))
    if not ids:return []
    selected=stratified(ids,lambda i:rows[i]['metric_norm'],8,1)
    selected |= stratified(ids,lambda i:rows[i]['quartic_bits'],8,1)
    selected |= tied_take(ids,lambda i:coverage_key(rows[i]['fingerprint']),4)
    # Two Pareto frontiers: shallow and deep exact centres are both useful
    # structural arms; small chart coefficients and high multiplicity preferred.
    for sign in (-1,1):
        def score(i):
            r=rows[i]
            return (sign*r['metric_norm'],r['quartic_bits'],r['quartic_max_bits'],-r['multiplicity'])
        scores={i:score(i) for i in ids}
        selected |= {i for i in ids if not any(all(a<=b for a,b in zip(scores[j],scores[i])) and scores[j]!=scores[i] for j in ids)}
    return [rows[i] for i in sorted(selected,key=lambda i:(-rows[i]['metric_norm'],rows[i]['fingerprint']))]
