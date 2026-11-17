"""Exact two-curvette blowups using QQ polynomial rational functions.

This is a finite control implementation of the retained local_cluster
algorithm. No series truncation, floating arithmetic or general termination
claim is used. Its user-supplied blowup guard is an implementation limit.
"""
from pathlib import Path
import sys
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_trefoil_complementary_fusion as F
C=F.C;J=F.J
def local_cluster(sigma,support,max_steps=300):
    z=S.Symbol('z');support=dict(support);p=max(support)
    assert sigma>=1 and all(j>=sigma+1 for j in support)
    def pair(num,den):
        nn=S.Poly(num,z,domain=S.QQ);dd=S.Poly(den,z,domain=S.QQ);assert not dd.is_zero
        gg=S.gcd(nn,dd);nn=nn.exquo(gg);dd=dd.exquo(gg);lc=dd.LC()
        return nn.mul_ground(1/lc),dd.mul_ground(1/lc)
    def order(a):
        nn,dd=a
        if nn.is_zero:return (10**9,S.Integer(0))
        ni=min(m[0] for m,c in nn.terms());di=min(m[0] for m,c in dd.terms())
        return ni-di,nn.nth(ni)/dd.nth(di)
    def divide(a,b):return pair(a[0]*b[1],a[1]*b[0])
    def subtract(a,c):return pair(a[0]-a[1].mul_ground(c),a[1])
    aa=[pair(z**p,sum(c*z**(p-j) for j,c in support.items())+param*z**(p+1-sigma)) for param in [3,7]]
    bb=[pair(z**sigma,1) for _ in [3,7]]
    ba,bc='Hx',None;path=('finite','0');rows=[]
    for step in range(max_steps):
        oa,_=order(aa[0]);ob,_=order(bb[0]);assert oa>0 and ob>0
        key=path;rows.append(dict(key=key,parents=[p for p in [ba,bc] if p is not None],multiplicity=min(oa,ob)))
        if ob<oa:aa,bb=bb,aa;ba,bc=bc,ba;direction='b'
        else:direction='a'
        bb=[divide(b,a) for a,b in zip(aa,bb)];ba,bc=key,bc
        ro=[order(b) for b in bb];residues=[c if n==0 else 0 for n,c in ro]
        if any(n<0 for n,c in ro) or residues[0]!=residues[1]:return dict(rows=rows,terminal=key)
        cb=residues[0]
        if cb!=0:bb=[subtract(b,cb) for b in bb];bc=None
        path=path+((direction,str(cb)),)
    raise AssertionError('finite control exceeded its blowup guard; no general termination claim')
def graph(sigma,support):
    """Full integral rho2 complementary source control, rooted at R."""
    z,T=S.symbols('z T');rr=J.retag(C.local_cluster(z**-1+T,z**2,0,z,T),'R','Hx')
    dd=J.retag(local_cluster(sigma,support),'D','Hy');out=C.assemble([rr,dd]);mat=out['full']
    vertices=set(range(mat.rows))-set(out['index']);arrow=out['labels'][rr['rows'][-1]['parents'][0]]
    spine=J.rooted_nodes(mat,vertices,arrow);adj={i:[j for j in vertices if i!=j and mat[i,j]] for i in vertices}
    nadj={i:[b if a==i else a for a,b in spine['edges'] if i in (a,b)] for i in spine['nodes']}
    ordernodes=[];last=None;cur=spine['root']
    while True:
        ordernodes.append(cur);nexts=[x for x in nadj[cur] if x!=last]
        if not nexts:break
        assert len(nexts)==1;last,cur=cur,nexts[0]
    assert set(ordernodes)==set(spine['nodes'])
    bridges=[]
    for a,b in zip(ordernodes,ordernodes[1:]):
        todo=[(a,[a])];seen={a};path=None
        while todo:
            cur,path0=todo.pop()
            if cur==b:path=path0;break
            for nxt in adj[cur]:
                if nxt not in seen:seen.add(nxt);todo.append((nxt,path0+[nxt]))
        assert path is not None;intern=path[1:-1]
        block=-mat.extract(intern,intern)
        bridges.append(dict(nodes=[a,b],D=int(block.det()),weights=[int(mat[i,i]) for i in intern]))
    legs={i:[dict(order=x['order'],cofactor=int(abs(mat.extract(x['vertices'][1:],x['vertices'][1:]).det())),weights=[int(mat[j,j]) for j in x['vertices']]) for x in spine['legs'][i]] for i in ordernodes}
    return dict(sigma=sigma,support=support,source=[[int(x) for x in row] for row in out['source'].tolist()],
        full_matrix=[[int(x) for x in row] for row in mat.tolist()],omitted_indices=out['index'],
        root=spine['root'],ordered_nodes=ordernodes,legs=legs,bridges=bridges,
        truncated_series_used=False)
