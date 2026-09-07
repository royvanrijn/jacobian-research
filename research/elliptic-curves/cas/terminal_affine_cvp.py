"""Exact k-nearest affine lattice enumeration across all parity classes."""
import heapq
from math import isqrt
from gmpy2 import mpq as Q


class AffineCVP:
    def __init__(self, gram):
        self.g = [[int(x) for x in row] for row in gram]
        self.n = n = len(self.g)
        self.mu = [[Q(0) for _ in range(n)] for _ in range(n)]
        self.d = []
        for i in range(n):
            for j in range(i):
                self.mu[i][j] = (self.g[i][j]-sum(self.mu[i][k]*self.mu[j][k]*self.d[k] for k in range(j)))/self.d[j]
            self.d.append(Q(self.g[i][i])-sum(self.mu[i][k]**2*self.d[k] for k in range(i)))
            assert self.d[-1]>0

    def nearest(self, target, count=64, node_limit=20000000):
        t = [Q(str(x)) for x in target]
        shifts = [t[i]+sum(self.mu[j][i]*t[j] for j in range(i+1,self.n)) for i in range(self.n)]
        def cost(w):
            return sum(self.d[i]*(w[i]+sum(self.mu[j][i]*w[j] for j in range(i+1,self.n))-shifts[i])**2 for i in range(self.n))
        seed = [0]*self.n
        for i in range(self.n-1,-1,-1):
            seed[i] = int(round(shifts[i]-sum(self.mu[j][i]*seed[j] for j in range(i+1,self.n))))
        seeds = {tuple(seed)}
        step = 0
        while len(seeds)<count:
            step += 1
            for i in range(self.n):
                for sign in (-1,1):
                    w=seed[:];w[i]+=sign*step;seeds.add(tuple(w))
        heap=[];seen={}
        def add(w,value):
            w=tuple(w)
            if w in seen: return
            seen[w]=value
            if len(heap)<count: heapq.heappush(heap,(-value,w))
            elif value < -heap[0][0]: heapq.heapreplace(heap,(-value,w))
        for w in seeds: add(w,cost(w))
        initial = str(-heap[0][0])
        w=[0]*self.n;nodes=0
        def visit(i,used):
            nonlocal nodes
            nodes+=1
            if nodes>node_limit: raise RuntimeError('affine CVP exact node limit exhausted')
            remaining=-heap[0][0]-used
            if remaining<0:return
            if i<0:
                add(w,used);return
            shift=sum(self.mu[j][i]*w[j] for j in range(i+1,self.n))-shifts[i]
            a,b=int(shift.numerator),int(shift.denominator)
            radius=remaining/self.d[i]*b*b
            r=isqrt(int(radius.numerator//radius.denominator))
            lo=-((r+a)//b);hi=(r-a)//b
            for value in sorted(range(lo,hi+1),key=lambda x:abs(Q(x)+shift)):
                w[i]=value;visit(i-1,used+self.d[i]*(Q(value)+shift)**2)
        visit(self.n-1,Q(0))
        bound=-heap[0][0]
        rows=sorted((v,w) for w,v in seen.items() if v<=bound)
        return {'rows':[{'distance':str(v),'word':list(w)} for v,w in rows],
            'initial_radius_squared':initial,'final_radius_squared':str(bound),'nodes':nodes,
            'requested_count':count,'all_boundary_ties_retained':True,
            'coverage':'Complete rational closed ellipsoid over Z^r, without any parity restriction. Certifies k-nearest rounded-metric words, not a global chart-height optimum.'}
