"""Fraction-only replay of the two integral, globally minimal neighbours."""
from fractions import Fraction as F
from functools import reduce
import hashlib
import json
from pathlib import Path
from operator import mul as product

from verify_pointed_height_bounds import add, mul, affine, homogeneous, require


def invariants(P,Q):
    d=add([4*c for c in P],mul(Q,Q));d += [F(0)]*(5-len(d))
    a,b,c,e,f=d[4],d[3],d[2],d[1],d[0]
    return 12*a*f-3*b*e+c*c,(72*a*c*f+9*b*c*e-27*a*e*e-27*b*b*f-2*c**3)/2


def matrix_mul(M,N):
    a,b,c,d=M;e,f,g,h=N
    return [a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h]


def padded(v,n):return v+[F(0)]*(n-len(v))


def verify(folder,root):
    base=json.loads((folder/'preconditioned_full-input.json').read_text())
    m=base['mapping'];P=list(map(F,m['reduced_P']));Q=list(map(F,m['reduced_Q']))
    A,B=map(F,base['curve'][3:]);c4,c6=-48*A/6**4,-864*B/6**6
    require(F(A+27,1296).denominator==1 and F(B+15*A-8289,46656).denominator==1,'nonintegral elliptic transport')
    source=root/base['prime_certificate']['path'];raw=source.read_bytes()
    require(hashlib.sha256(raw).hexdigest()==base['prime_certificate']['sha256'],'prime packet changed')
    proof=json.loads(raw);factors=proof['records']['302']['factorizations']['DISCRIMINANT_FACTORIZATION']
    require((c4**3-c6**2)/1728==reduce(product,(int(p)**e for p,e in factors),1),'minimal discriminant product differs')
    # Nonminimal integral Weierstrass equations require v(c4)>=4 and
    # v(Delta)>=12. Only2 reaches the latter threshold; c4 is odd there.
    require(all(e<12 or (p=='2' and c4.numerator%2) for p,e in factors),'minimality gate missing')
    require(invariants(P,Q)==(c4,c6),'base quartic is not level zero')
    rows=[];keys=[]
    for file in sorted((folder/'neighbours').glob('neighbour-*-input.json')):
        data=json.loads(file.read_text());n=data['mapping'];w=data['neighbour_witness']
        p,r,s=map(int,(w['prime'],w['abscissa_residue'],w['ordinate_residue']))
        require(p>2 and str(p) in proof['proved_primes'] and 0<=r<p and 0<=s<p,'unproved neighbour prime')
        newP=[c/p**2 for c in add(affine(P,r,p),add([-s*c for c in affine(Q,r,p)],[-s*s]))]
        newQ=[c/p for c in add(affine(Q,r,p),[2*s])]
        require(all(c.denominator==1 for c in newP+newQ),'nonintegral local neighbour')
        require(padded(newP,5)==list(map(F,w['unreduced_P'])) and padded(newQ,3)==list(map(F,w['unreduced_Q'])),'neighbour equation differs')
        require(invariants(newP,newQ)==(c4,c6),'neighbour invariants changed')
        finalP=list(map(F,n['reduced_P']));finalQ=list(map(F,n['reduced_Q']))
        require(all(c.denominator==1 for c in finalP+finalQ) and invariants(finalP,finalQ)==(c4,c6),'final model not integral level zero')
        M=matrix_mul(list(map(F,m['matrix'])),matrix_mul([p,r,0,1],list(map(F,w['reduction_matrix']))))
        require(M==list(map(F,n['matrix'])),'composed parameter transport changed')
        ratio=F(n['square_ratio']);from math import isqrt
        disc=padded(add([4*c for c in finalP],mul(finalQ,finalQ)),5)
        require(ratio>0 and isqrt(ratio.numerator)**2==ratio.numerator and isqrt(ratio.denominator)**2==ratio.denominator,'ordinate scale not a square')
        require(homogeneous(list(map(F,m['raw_coefficients'])),M)==[ratio*c for c in disc],'quartic equivalence failed')
        # Reimplement the four height-preserving projective signed permutations.
        def key(M):
            a,b,c,d=M;variants=[M,[-a,b,-c,d],[b,a,d,c],[-b,a,-d,c]]
            return min(tuple(c/next(z for z in v if z) for c in v) for v in variants)
        k=key(M);require(k!=key(list(map(F,m['matrix']))) and k not in keys,'duplicate bounded-box key')
        keys.append(k);rows.append({'file':file.name,'prime':p,'residue':r,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    require(len(rows)==2,'incomplete frozen pair')
    return {'status':'PASS_FRACTION_MINIMALITY_AND_TRANSPORT_REPLAY','models':rows,
        'boundary':'Global minimality follows from integral equations with the same discriminant as the verified minimal elliptic model. Different signed-permutation box keys do not prove different global GL2(Z) orbits.'}


if __name__=='__main__':
    import argparse,time
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();require(not args.output.exists(),'preserve replay')
    start=time.process_time();r=verify(args.folder,Path(__file__).resolve().parents[2])
    r.update(cpu_seconds=time.process_time()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
