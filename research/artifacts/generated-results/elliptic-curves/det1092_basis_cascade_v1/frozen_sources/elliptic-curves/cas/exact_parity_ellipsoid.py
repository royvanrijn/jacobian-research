"""Exhaustive rational LDL enumeration in an integral parity coset."""
from fractions import Fraction as F
from math import isqrt


def ldl(gram):
    n=len(gram)
    if any(len(r)!=n for r in gram) or any(gram[i][j]!=gram[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Gram must be square and symmetric')
    L=[[F(i==j) for j in range(n)] for i in range(n)];D=[]
    for j in range(n):
        d=F(gram[j][j])-sum(L[j][k]**2*D[k] for k in range(j))
        if d<=0:raise ValueError('Gram is not positive definite')
        D.append(d)
        for i in range(j+1,n):
            L[i][j]=(F(gram[i][j])-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    if any(sum(L[i][k]*D[k]*L[j][k] for k in range(n))!=gram[i][j]
           for i in range(n) for j in range(n)):
        raise ArithmeticError('exact LDL reconstruction failed')
    return L,D


def enumerate_coset(gram,residue,bound):
    """Return all norm counts and a minimum witness within the closed ellipsoid.

    At descending coordinate i, LDL gives exactly D_i*(z_i-centre)^2 plus
    the fixed nonnegative tail norm. A conservative integral interval contains
    every allowed z_i; exact inequalities and parity filter it. No float, CVP
    oracle, pruning score, or unproved lattice reduction property is used.
    """
    n=len(gram);L,D=ldl(gram);bound=F(bound)
    if len(residue)!=n or any(r not in (0,1) for r in residue) or bound<0:
        raise ValueError('invalid parity residue or radius')
    z=[0]*n;nodes=[0]*n;counts={};best=None;witness=None
    def visit(i,partial):
        nonlocal best,witness
        if i<0:
            if partial.denominator!=1:raise ArithmeticError('integral Gram produced fractional norm')
            q=int(partial);counts[q]=counts.get(q,0)+1
            if best is None or q<best:best=q;witness=list(z)
            return
        centre=-sum(L[j][i]*z[j] for j in range(i+1,n))
        radius_squared=(bound-partial)/D[i]
        if radius_squared<0:return
        radius=isqrt(radius_squared.numerator//radius_squared.denominator)+1
        lo=-((-(centre-radius).numerator)//(centre-radius).denominator)
        hi=(centre+radius).numerator//(centre+radius).denominator
        lo+=(residue[i]-lo)%2
        for value in range(lo,hi+1,2):
            nodes[i]+=1;total=partial+D[i]*(F(value)-centre)**2
            if total<=bound:z[i]=value;visit(i-1,total)
    visit(n-1,F(0))
    return {'bound':str(bound),'norm_counts':{str(k):counts[k] for k in sorted(counts)},
            'minimum_within_bound':best,'minimum_witness':witness,'visited_candidates_by_depth':nodes}
