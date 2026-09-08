"""Finite-field six-pair recognition without scanning every projective pole.

Each two disjoint equal-gap pairs determines a quadratic in the Mobius pole.
Six disjoint pairs give fifteen pair-of-pair coincidences with the same pole
and unsigned gap. Count these keys, then verify an actual six-edge matching.
This is a necessary finite-field test, never a rational construction proof.
"""
from functools import lru_cache
from itertools import combinations
import numpy as np


@lru_cache(maxsize=24)
def equations(size):
    ids=[]
    for i,j,k,l in combinations(range(size),4):
        for a,b,c,d in [(i,j,k,l),(i,k,j,l),(i,l,j,k)]:
            for sign in [-1,1]:ids.append((a,b,c,d,sign))
    return np.asarray(ids,dtype=np.int64).T


def matching(edges,count,chosen=()):
    if count==0:return chosen
    if len(edges)<count:return None
    for pos,(i,j) in enumerate(edges):
        rest=[(k,l) for k,l in edges[pos+1:] if i not in (k,l) and j not in (k,l)]
        found=matching(rest,count-1,chosen+((i,j),))
        if found is not None:return found
    return None


class PairCollisionFilter:
    def __init__(self,prime):
        p=int(prime)
        if p<5 or p>2000000:raise ValueError('Odd prime in 5..2000000 required')
        # The caller pins and verifies primality; inversion detects composites.
        inv=[0,1]+[0]*(p-2)
        for a in range(2,p):inv[a]=(-(p//a)*inv[p%a])%p
        self.inv=np.asarray(inv,dtype=np.int64)
        assert np.all((np.arange(1,p,dtype=np.int64)*self.inv[1:])%p==1)
        roots=np.full(p,-1,dtype=np.int64);v=np.arange((p+1)//2,dtype=np.int64);roots[v*v%p]=v
        self.sqrt=roots;self.p=p

    def __call__(self,z,witness=False):
        if z is None:return None
        n=len(z);p=self.p
        if n<12 or n>=p:raise ValueError('Need 12 <= number of points < p')
        if len(set(z))!=n:return None
        offset=0
        while offset in z:offset+=1
        w=np.asarray([0 if a is None else self.inv[(int(a)-offset)%p] for a in z],dtype=np.int64)
        i,j,k,l,sign=equations(n);a,b,c,d=w[i],w[j],w[k],w[l]
        ab=(a-b)%p;cd=sign*(c-d)%p
        A=(ab-cd)%p;B=(-ab*(c+d)+cd*(a+b))%p
        C=((ab*c%p)*d-(cd*a%p)*b)%p
        keys=[]
        def retain(poles,indices):
            aa=a[indices];bb=b[indices]
            den=(aa-poles)*(bb-poles)%p;good=den!=0
            poles=poles[good];indices=indices[good];den=den[good]
            gap=ab[indices]*self.inv[den]%p;gap=np.minimum(gap,p-gap)
            keys.append(poles*p+gap)
        disc=(B*B-4*(A*C%p))%p;root=self.sqrt[disc]
        good=np.flatnonzero((A!=0)&(root>=0))
        for sign in [-1,1]:retain(((-B[good]+sign*root[good])*self.inv[2*A[good]%p])%p,good)
        linear=np.flatnonzero((A==0)&(B!=0))
        retain((-C[linear]*self.inv[B[linear]])%p,linear)
        infinity=np.flatnonzero(A==0);gap=np.minimum(ab[infinity],p-ab[infinity]);keys.append(p*p+gap)
        values,counts=np.unique(np.concatenate(keys),return_counts=True)
        ii,jj=np.triu_indices(n,1);hits=[]
        for key in values[counts>=15]:
            pole,gap=divmod(int(key),p)
            if pole==p:xs=w;valid=np.ones(n,dtype=bool)
            else:
                delta=(w-pole)%p;valid=delta!=0;xs=self.inv[delta]
            differences=(xs[ii]-xs[jj])%p;differences=np.minimum(differences,p-differences)
            selected=np.flatnonzero(valid[ii]&valid[jj]&(differences==gap))
            found=matching([(int(ii[q]),int(jj[q])) for q in selected],6)
            if found is not None:
                if not witness:return True
                hits.append({'finite_chart_offset':offset,'pole':None if pole==p else pole,'gap':gap,'pairs':[list(edge) for edge in found]})
        return hits if witness else False
