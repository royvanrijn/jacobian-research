"""Exact pair triage and outcome-independent scoring policy (no Sage dependency)."""
from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt, log


def require(condition, message):
    if not condition:
        raise ValueError(message)


def square_root(value):
    value = Q(str(value))
    if value < 0:
        return None
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    return Q(a, b) if a*a == value.numerator and b*b == value.denominator else None


def same_quadratic_class(first, second):
    """Return c with q_second=c^2*q_first, or None; retain nonsquare constants."""
    a, b = list(map(lambda x: Q(str(x)), first)), list(map(lambda x: Q(str(x)), second))
    while a and not a[-1]: a.pop()
    while b and not b[-1]: b.pop()
    if not a or len(a) != len(b): return None
    ratio = b[-1]/a[-1]
    return square_root(ratio) if all(y == ratio*x for x,y in zip(a,b)) else None


def dot(first, gram, second):
    require(len(first) == len(second) == len(gram), 'incompatible lattice dimensions')
    return sum(int(a)*int(g)*int(b) for a,row in zip(first,gram) for g,b in zip(row,second))


def pair_order(rows, gram, limit):
    """Same-cover attempts first, then odd intersections; no outcome labels read."""
    require(limit > 0, 'pair limit must be positive')
    require(len({r['mask'] for r in rows}) == len(rows), 'duplicate orbit masks')
    for r in rows:
        require(dot(r['word'], gram, r['word']) == 10, 'not a norm-ten trace')
    same, odd = [], []
    for i,a in enumerate(rows):
        for b in rows[i+1:]:
            scale = same_quadratic_class(a['q'], b['q'])
            intersection = 8-dot(a['word'], gram, b['word'])
            row = {'masks':[a['mask'],b['mask']], 'intersection':intersection,
                   'coefficient_bits':a['coefficient_bits']+b['coefficient_bits']}
            if scale is not None:
                same.append({**row,'kind':'same-cover','scale':str(scale)})
            elif intersection > 0 and intersection % 2:
                odd.append({**row,'kind':'odd-intersection'})
    key = lambda r:(r['coefficient_bits'],abs(r['intersection']),r['masks'])
    # A same-cover prefix must not consume the fallback's entire finite budget.
    return sorted(same,key=key)[:limit] + sorted(odd,key=key)[:limit]


def hash_order(row, salt='x1092-elkies-controls-v1'):
    return sha256((salt+'/'+str(row['id'])).encode()).hexdigest(), str(row['id'])


def effective_control_count(eligible, requested):
    """Preserve at least one ranked fibre whenever any eligible fibre exists."""
    require(eligible >= 0 and requested >= 0, 'negative control population')
    return min(requested, max(0, eligible-1))


def carrier_quality_summary(height_bits, primary_cap):
    """Outcome-free carrier quality from induced parameter heights only."""
    require(primary_cap > 0, 'nonpositive carrier height cap')
    values=sorted(int(v) for v in height_bits)
    require(all(v >= 0 for v in values), 'negative parameter height')
    extended=4*primary_cap
    if values:
        median=values[len(values)//2]
        worst=values[-1]
    else:
        median=worst=10**18
    return {
        'sampled':len(values),
        'within_primary':sum(v <= primary_cap for v in values),
        'within_extended':sum(v <= extended for v in values),
        'primary_cap':primary_cap,
        'extended_cap':extended,
        'median_bits':median,
        'max_bits':worst,
    }


def carrier_quality_key(summary, same_cover=False):
    """Prefer usable yield first; same-cover is only a tie-break, never a blocker."""
    return (
        int(summary['within_primary']),
        int(summary['within_extended']),
        int(bool(same_cover)),
        -int(summary['median_bits']),
        -int(summary['max_bits']),
        int(summary['sampled']),
    )


def split_controls(rows, count):
    require(0 <= count < len(rows), 'controls must leave a nonempty ranked arm')
    require(len({r['id'] for r in rows}) == len(rows), 'duplicate candidate IDs')
    selected = sorted(rows,key=hash_order)[:count]
    ids = {r['id'] for r in selected}
    return selected, [r for r in rows if r['id'] not in ids]


def shortlist(rows, controls, keep):
    require(keep > 0, 'deep shortlist must be nonempty')
    ids = {r['id'] for r in controls}
    eligible = [r for r in rows if r['id'] not in ids]
    return sorted(eligible,key=lambda r:(-r['shallow']['score_units'],r['j_bits'],str(r['id'])))[:keep]


def final_selection(deep_rows, controls, ranked):
    require(ranked > 0, 'ranked count must be positive')
    ids = {r['id'] for r in controls}
    pool = [r for r in deep_rows if r['id'] not in ids]
    best = sorted(pool,key=lambda r:(-r['deep']['score_units'],r['j_bits'],str(r['id'])))[:ranked]
    # Interleave; controls never acquire a rank-based refill or promotion rule.
    result=[]
    for i in range(max(len(best),len(controls))):
        if i < len(best): result.append({**best[i],'control':False})
        if i < len(controls): result.append({**controls[i],'control':True})
    return result


def mestre_score(traces):
    """Only a scheduling score. None is a missing term, never trace/rank zero."""
    rows=list(traces)
    require(len({p for p,a in rows}) == len(rows), 'duplicate scoring primes')
    used=[(int(p),int(a)) for p,a in rows if a is not None]
    for p,a in used:
        require(p >= 5 and p+1-a > 0 and a*a <= 4*p, 'invalid good-fibre trace')
    return {'score_units':sum(round(10**12*(2-a)*log(p)/(p+1-a)) for p,a in used),
            'used_primes':len(used),'missing_primes':[int(p) for p,a in rows if a is None],
            'rank_bound':None}


def rational_height(value):
    value=Q(str(value))
    return max(abs(value.numerator).bit_length(),value.denominator.bit_length())


def nth_root(value, n):
    value=Q(value)
    if value < 0: return None
    def root(k):
        lo,hi=0,1 << ((k.bit_length()+n-1)//n)
        while lo < hi:
            mid=(lo+hi+1)//2
            if mid**n <= k: lo=mid
            else: hi=mid-1
        return lo if lo**n == k else None
    a,b=root(value.numerator),root(value.denominator)
    return None if a is None or b is None else Q(a,b)


def short_isomorphic(a,b,c,d):
    """Exact Q-isomorphism of nonsingular short models; equal j is insufficient."""
    a,b,c,d=map(lambda x:Q(str(x)),(a,b,c,d))
    require(4*a**3+27*b*b != 0 and 4*c**3+27*d*d != 0, 'singular model')
    if not a or not c:
        return not a and not c and nth_root(d/b,6) is not None
    if not b or not d:
        return not b and not d and nth_root(c/a,4) is not None
    u2=(d/b)/(c/a)
    return square_root(u2) is not None and u2*u2 == c/a and u2**3 == d/b
