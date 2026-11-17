"""Exact bounded controls for the mixed Eisenstein ramified-prime calculation."""
import argparse
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path


def lattice_counts(limit):
    bound = isqrt(4 * limit // 3)
    counts = [0] * (limit + 1)
    for x in range(-bound, bound + 1):
        for y in range(-bound, bound + 1):
            value = x*x - x*y + y*y
            if value <= limit:
                counts[value] += 1
    return counts, bound


def reduced_and_omega(n):
    remaining = n
    omega = 0
    d = 2
    while d*d <= remaining:
        if remaining % d == 0:
            if d % 3 != 2:
                return False, None
            while remaining % d == 0:
                remaining //= d
                omega += 1
        d += 1
    if remaining > 1:
        if remaining % 3 != 2:
            return False, None
        omega += 1
    return True, omega


def sigma(n):
    return sum(d + (n//d if d*d != n else 0)
               for d in range(1, isqrt(n)+1) if n % d == 0)


def direct_counts(n, r):
    total = sum(r[a]*r[n-a] for a in range(n+1))
    odd = sum(r[a]*r[n-a] for a in range(1,n+1,2))
    bad = sum(r[a]*r[n-a] for a in range(n+1)
              if a % 9 == 0 or (n-a) % 9 == 0)
    bad_odd = sum(r[a]*r[n-a] for a in range(1,n+1,2)
                  if a % 9 == 0 or (n-a) % 9 == 0)
    return total, odd, bad, bad_odd


def verify_equal(actual, expected):
    if actual != expected:
        raise ValueError(f'count mismatch: {actual} != {expected}')


def group_control():
    def rmul(z,u):
        a,b=z; c,d=u
        return ((a*c-b*d)%3,(a*d+b*c-b*d)%3)
    def bar(z):
        a,b=z
        return ((a-b)%3,(-b)%3)
    def add(z,w):
        return tuple((a+b)%3 for a,b in zip(z,w))
    def neg(z):
        return tuple((-a)%3 for a in z)
    def mul(q,r):
        z,w=q[:2],q[2:]; u,v=r[:2],r[2:]
        return add(rmul(z,u),neg(rmul(w,bar(v))))+add(rmul(z,v),rmul(w,bar(u)))
    def norm(q):
        a,b,c,d=q
        return (a*a-a*b+b*b+c*c-c*d+d*d)%3
    def canon(q):
        return min(q,neg(q))
    units=[q for q in product(range(3),repeat=4) if norm(q)]
    group={canon(q) for q in units}
    runits=[(1,0),(2,0),(0,1),(0,2),(2,2),(1,1)]
    gamma={canon(z+(0,0)) for z in runits}|{canon((0,0)+z) for z in runits}
    def coset(q):
        return min(canon(mul(g,q)) for g in gamma)
    cosets={coset(q) for q in group}
    verify_equal((len(units),len(group),len(gamma),len(cosets)),(72,36,6,6))
    for q in group:
        for r in group:
            verify_equal(norm(mul(q,r)),norm(q)*norm(r)%3)
    opposite={q for q in cosets if norm(q)==2}
    verify_equal(len(opposite),3)
    omega=(0,1,0,0)
    for q in opposite:
        orbit={coset(q),coset(mul(q,omega)),coset(mul(mul(q,omega),omega))}
        verify_equal(orbit,opposite)
    return {'units':72,'projective_group':36,'unit_image':6,'cosets':6,
            'norm_classes':[3,3],'right_omega_cycles_opposite_class':True,
            'multiplicativity_checks':len(group)**2}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--limit',type=int,default=3000)
    parser.add_argument('--output',type=Path,default=Path('ramified_counts_receipt.json'))
    args=parser.parse_args()
    if args.limit<16:
        parser.error('limit must be at least16')
    r,bound=lattice_counts(args.limit)
    cases=0; squares=0; summaries=[]
    for n in range(4,args.limit+1,6):
        reduced,omega=reduced_and_omega(n)
        if not reduced:
            continue
        cases+=1
        lam=0
        if isqrt(n)**2==n:
            squares+=1
            lam=(-1)**(omega//2)*isqrt(n)
        S=sigma(n)-(sigma(n//4) if n%4==0 else 0)
        actual=direct_counts(n,r)
        expected=(12*sigma(n),12*S,4*sigma(n)+8*lam,4*S+12*lam)
        verify_equal(actual,expected)
        if actual[3]*3*n > actual[1]*(n+2*isqrt(n)) and lam:
            # Compare only the exact square case here; nonsquares have fraction1/3.
            raise ValueError('Ramified odd bound failed')
        if n<=100:
            summaries.append({'n':n,'lambda':lam,'S':S,'counts':actual})
    wrong_sign_rejected=False
    try:
        verify_equal(direct_counts(4,r)[2],4*sigma(4)-8*(-2))
    except ValueError:
        wrong_sign_rejected=True
    if not wrong_sign_rejected:
        raise ValueError('Wrong square sign was accepted')
    receipt={'scope':'Exact finite controls only; no unrestricted identity proved by enumeration',
             'limit':args.limit,'coordinate_bound':bound,
             'coordinate_box_complete_reason':'Q(x,y)>=3*x*x/4 and>=3*y*y/4',
             'reduced_even_n1mod3_cases':cases,'square_cases':squares,
             'all_four_counts_match':True,'small_examples':summaries,
             'wrong_lambda_sign_rejected_n4':True,'group_control':group_control(),
             'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='small_examples'},indent=2))

if __name__=='__main__':
    main()
