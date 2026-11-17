"""All unfilled peripheral rows for a finite iterated-knot cover.

Longitude monodromy must be identity. Every meridian cycle then gives
one boundary torus with its actual positive normal meridian and lifted
preferred longitude. No word expansion or branched filling is performed.
"""
from pathlib import Path
import sys
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_loop_first_arrival as K
E=K.E
def reduced_h1(rows,extra=()):
    """Cancel exact unit relations before the generic Smith computation.

    Each Schur step is an integral unimodular row/column elimination.
    Removing its pivot relation and generator leaves the presented group
    unchanged. This avoids generic Smith coefficient growth on long covers.
    """
    original=rows+list(extra)
    # HNF of the transposed relation matrix is an exact basis for the
    # SAME row lattice, eliminating redundant relations before Smith.
    block=[list(map(int,row)) for row in hermite_normal_form(S.Matrix(original).T).T.tolist()]
    if not block:return dict(free=len(original[0]),torsion=[])
    width=len(block[0]);assert all(len(row)==width for row in block)
    while block and width:
        rowcounts=[sum(x!=0 for x in row) for row in block]
        colcounts=[sum(row[j]!=0 for row in block) for j in range(width)]
        pivots=[(rowcounts[i]+colcounts[j],i,j) for i,row in enumerate(block) for j,x in enumerate(row) if abs(x)==1]
        if not pivots:break
        _,i,j=min(pivots);pivot=block[i];unit=pivot[j]
        block=[[row[k]-row[j]*unit*pivot[k] for k in range(width) if k!=j] for h,row in enumerate(block) if h!=i]
        width-=1
    if not width:return dict(free=0,torsion=[])
    if not block:return dict(free=width,torsion=[])
    return E.h1(block)
def peripheral_cover(p,q,cables,rho,reverse=False):
    degree=len(rho[0]);count=len(rho);width=degree*count;ident=tuple(range(degree))
    assert count==2+len(cables)
    assert all(sorted(pi)==list(ident) for pi in rho)
    unit=(ident,[[0]*width for _ in range(degree)])
    def compose(a,b):
        da,fa=a;db,fb=b
        return tuple(db[da[s]] for s in range(degree)),[[x+y for x,y in zip(fa[s],fb[da[s]])] for s in range(degree)]
    def invert(a):
        dest,flows=a;back=E.inv(dest)
        return back,[[-x for x in flows[back[s]]] for s in range(degree)]
    def power(a,k):
        if k<0:return power(invert(a),-k)
        out=unit
        while k:
            if k%2:out=compose(out,a)
            a=compose(a,a);k//=2
        return out
    def prod(*parts):
        out=unit
        for a in parts:out=compose(out,a)
        return out
    generators=[]
    for g,pi in enumerate(rho):
        flows=[[0]*width for _ in range(degree)]
        for s in range(degree):flows[s][g*degree+s]=1
        generators.append((E.inv(pi),flows))
    u,v,gg=S.gcdex(q,p);assert gg==1
    meridian=prod(power(generators[0],int(u)),power(generators[1],int(v)))
    longitude=prod(power(generators[0],p),power(meridian,-p*q))
    relations=[prod(power(generators[0],p),power(generators[1],-q))]
    for j,(r,k) in enumerate(cables,2):
        u,v,gg=S.gcdex(r,k);assert gg==1;u=int(u);v=-int(v)
        fibre=prod(power(meridian,k),power(longitude,r))
        relations.append(prod(power(generators[j],r),invert(fibre)))
        meridian=prod(power(meridian,u),power(longitude,v),power(generators[j],-v))
        longitude=prod(power(meridian,-r*k),fibre)
    assert all(a[0]==ident for a in relations)
    assert longitude[0]==ident,'requires actual longitude monodromy1'
    edges=[(g,s) for g in range(count) for s in range(degree)]
    parent=list(range(degree));tree=set()
    def root(s):
        while parent[s]!=s:s=parent[s]
        return s
    for g,s in (list(reversed(edges)) if reverse else edges):
        t=E.inv(rho[g])[s]
        if root(s)!=root(t):parent[root(s)]=root(t);tree.add((g,s))
    assert len(tree)==degree-1,'requires connected cover'
    keep=[e for e in edges if e not in tree];indices=[g*degree+s for g,s in keep]
    project=lambda row:[row[j] for j in indices]
    rows=[project(a[1][s]) for a in relations for s in range(degree)]
    seen=set();boundaries=[]
    for s in range(degree):
        if s in seen:continue
        cycle=[];t=s
        while t not in cycle:cycle.append(t);seen.add(t);t=meridian[0][t]
        normal=power(meridian,len(cycle));assert normal[0][s]==s
        boundaries.append(dict(sheets=cycle,degree=len(cycle),meridian=project(normal[1][s]),longitude=project(longitude[1][s])))
    return dict(degree=degree,relation_rows=rows,boundaries=boundaries,keep=keep,
        meridian_monodromy=E.inv(meridian[0]),longitude_monodromy=longitude[0],
        word_expansion_used=False)

def integral_null_core_linking(data):
    """Preferred framing matrix, requiring ALL filling cores integrally zero.

    This is a reconstruction from actual unfilled peripheral rows, not a
    general rational-linking algorithm for nonzero torsion core classes.
    """
    rows=data['relation_rows'];boundaries=data['boundaries'];count=len(boundaries)
    meridians=[b['meridian'] for b in boundaries];longitudes=[b['longitude'] for b in boundaries]
    closed=rows+meridians;closed_h1=reduced_h1(closed)
    assert all(reduced_h1(closed,[lon])==closed_h1 for lon in longitudes),'requires INTEGRAL zero for every filling core'
    basis=S.Matrix.hstack(*S.Matrix(rows).nullspace())
    mu=S.Matrix(meridians)*basis;lon=S.Matrix(longitudes)*basis
    assert mu.rank()==count,'requires primitive independent meridians'
    coefficients=[]
    for i in range(count):
        solutions=list(S.linsolve((mu.T,lon.row(i).T)));assert len(solutions)==1
        coefficient=solutions[0];assert all(x.is_Integer for x in coefficient)
        coefficients.append([int(x) for x in coefficient])
    matrix=S.Matrix(coefficients);relations=[]
    for i,longitude in enumerate(longitudes):
        relations.append([longitude[j]-sum(coefficients[i][k]*meridians[k][j] for k in range(count)) for j in range(len(longitude))])
    assert reduced_h1(rows,relations)==reduced_h1(rows),'free calculation must also hold integrally through torsion'
    assert matrix==matrix.T
    assert matrix*S.Matrix([b['degree'] for b in boundaries])==S.zeros(count,1)
    return dict(matrix=coefficients,closed_H1=closed_h1,all_cores_integrally_zero=True,
        preferred_diagonal_framings_retained=True)
