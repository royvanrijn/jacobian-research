"""Certify a univariate eliminant for the normalized ratio ell/v."""

from collections import defaultdict

from sage.all import GF, QQ, PolynomialRing, block_diagonal_matrix, matrix, vector

from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows,weighted_exponents

base=PolynomialRing(QQ,"s"); K=base;s=K.gen()
P=PolynomialRing(K,names=("r","h","w","x","t"))
r,h,w,x,t=P.gens()
def homrows(rows):
    homfs=[]
    for k,row in enumerate(rows,7):
        deg=(k-1)//2
        co=defaultdict(K.zero)
        for (ri,li,hi,vi,wi),coef in row.dict().items():
            assert (li+vi)%2==(k-1)%2
            xi=(li+vi)//2
            ti=deg-ri-hi-wi-xi
            assert ti>=0
            co[(ri,hi,wi,xi,ti)] += K(coef)*s**li
        homfs.append(P(co))
    return homfs

def localrows(fs,chart,order):
    B=PolynomialRing(K,names=("a","b","c","d"))
    a,b,c,d=B.gens()
    if chart==1: cmap=P.hom([a,b,c+s*(8-s)/17,B(1),d],B); pivots=(0,2)
    elif chart==2: cmap=P.hom([a,b+17,B(1),c,d],B); pivots=(0,)
    local=[jet(cmap(f),order) for f in fs]; assert all(not f.constant_coefficient() for f in local)
    return local,pivots

def wsum(e):
    return sum(e) if len(e)!=2 else e[0]+2*e[1]

def jet(f,d,hom=False):
    check=(lambda v:v==d) if hom else (lambda v:v<=d)
    return f.parent()({e:c for e,c in f.dict().items() if check(wsum(e))})

def products(pows,exponents,cut):
    a=pows[0][0]
    for i,n in enumerate(exponents):
        if n:
            a=jet(a*pows[i][n],cut)
        if not a: break
    return a

def substitute(fs,vals,cut,B):
    maxs=[max([max((e[i] for e in f.dict()),default=0) for f in fs]) for i in range(len(vals))]
    powers=[[B.one()] for _ in vals]
    for val,m,pow in zip(vals,maxs,powers):
        for j in range(m): pow.append(jet(pow[-1]*val,cut))
    return [sum(c*products(powers,e,cut) for e,c in f.dict().items()) for f in fs]

def local_eliminate(fs,pivots,order):
    orig=fs[0].parent(); free=[j for j in range(orig.ngens()) if j not in pivots]
    B=PolynomialRing(K,names=tuple(f"u{i}" for i in free)); gens=B.gens()
    J=matrix(K,[[f.monomial_coefficient(orig.gen(j)) for j in pivots] for f in fs])
    chosen=J.change_ring(K.fraction_field()).transpose().pivots()
    assert len(chosen)==len(pivots)
    jac=J.matrix_from_rows(chosen); d_jac=jac.det(); adj=jac.adjugate()
    dvals=[B.zero()]*4
    for j,u in zip(free,gens,strict=True):dvals[j]=d_jac**2*u
    pivot_fs=[fs[j] for j in chosen]
    for q in range(1,order+1):
        residual = vector(B, [jet(f,q,hom=True) for f in substitute(pivot_fs,dvals,q,B)])
        corrections=adj.change_ring(B)*residual
        for j,a in zip(pivots,corrections,strict=True):
            corrected={}
            for e,c in a.dict().items():
                qx,rem=c.quo_rem(d_jac); assert not rem
                corrected[e]=qx
            dvals[j]-=B(corrected)
    raw=substitute(fs,dvals,order,B)
    assert all(not raw[j] for j in chosen)
    return dvals,raw

def concatenate_J(raw,degree):
    R=raw[0].parent()
    n=R.ngens(); weights=(1,2) if n==2 else (1,)*n; mon=list(weighted_exponents(weights,degree)); ix={e:i for i,e in enumerate(mon)}
    parts=[]
    for f in raw:
        if not f: continue
        d=min(wsum(e) for e in f.dict())
        for a in weighted_exponents(weights,degree-d): parts.append((a,f.dict()))
    A=matrix(K,len(parts),len(mon))
    for i,(a,items) in enumerate(parts):
        for e,c in items.items():
            m=tuple(x+y for x,y in zip(a,e))
            if wsum(m)<=degree:A[i,ix[m]]=c
    return A,mon


def mod_eval(coef, value, finite):
    return finite(coef(QQ(value)))

def local_restriction_matrix(chart, vals, mon, mon9, finite, s0):
    B=PolynomialRing(finite,names=vals[0].parent().variable_names())
    dvals=[B({e:mod_eval(c,s0,finite) for e,c in g.dict().items()}) for g in vals]
    a,b,c,d=dvals
    if chart==1: mapvals=[a,b,c+finite(s0*(8-s0))/17,B.one(),d]
    if chart==2: mapvals=[a,b+17,B.one(),c,d]
    S=PolynomialRing(finite, names=P.variable_names())
    mons=[S({e:finite.one()}) for e in mon9]
    cut=max(wsum(e) for e in mon)
    return matrix(finite, [list(B(g).monomial_coefficient(B({e:finite.one()})) for e in mon)
                            for g in substitute(mons,mapvals,cut,B)])

def certify(rows):
    fs=homrows(rows)
    finite=GF(65521); s0=1
    mon9=[e for e in weighted_exponents((1,)*5,9) if sum(e)==9]
    assert len(mon9)==715
    isom={}
    for chart,order,expected in [(2,4,20),(1,8,10)]:
        loc,piv=localrows(fs,chart,order)
        vals,raw=local_eliminate(loc,piv,order)
        A,mon=concatenate_J(raw,order)
        rank=A.change_ring(K.fraction_field()).rank()
        assert rank==expected and A.ncols()-rank==15
        mod=matrix(finite,A.nrows(),A.ncols(),{
            k:mod_eval(c,s0,finite) for k,c in A.dict().items()})
        assert mod.rank()==expected
        eig=local_restriction_matrix(chart,vals,mon,mon9,finite,s0)
        isom[chart]=(mod,eig)
    target=(0,0,0,0,9)
    mod1,img1=isom[1]; mod2,img2=isom[2]
    rel=block_diagonal_matrix(mod1,mod2)
    img=img1.augment(img2)
    assert not any(img[mon9.index(target)])
    assert rel.stack(img).rank()==60
    print("Two quotient restriction map ranks total thirty",flush=True)

    pos={e:i for i,e in enumerate(mon9)}
    labels=[]
    mo=[]
    sd=[]
    for i,f in enumerate(fs):
        deg=(i+6)//2
        assert all(sum(e)==deg for e in f.dict())
        for exp in weighted_exponents((1,)*5,9-deg):
            if sum(exp)!=9-deg: continue
            labels.append((i,exp))
        mo.append([(e,mod_eval(c,s0,finite)) for e,c in f.dict().items()])
        sd.append(max(c.degree() for c in f.dict().values()))
    M=matrix(finite,len(labels),len(mon9))
    for i,(k,a) in enumerate(labels):
        for ex,c in mo[k]:
            j=pos[tuple(x+y for x,y in zip(a,ex))]
            M[i,j]=c
    I=M.transpose().pivots()
    J=M.pivots()
    assert len(I)==len(J)==685
    d0=M.matrix_from_rows_and_columns(I,J).det()
    assert d0 == 63052
    deg_bound=2*(sum(mon9[j][3] for j in J)-sum(labels[i][1][3] for i in I))+sum(labels[i][0]%2 for i in I)
    assert deg_bound == 1262 <= sum(sd[labels[i][0]] for i in I)
    print("Certified projection pivot minor",d0,"degree bound",deg_bound,flush=True)

    # Sort by the bounds for contributions to the degree in s.
    # Independent pivot sets on each side still give a nonsingular minor.
    rows_order = sorted(range(len(labels)), key=lambda i: labels[i][0] % 2 - 2*labels[i][1][3])
    cols_order = sorted(range(len(mon9)), key=lambda j: mon9[j][3])
    rows = [rows_order[i] for i in M.transpose().matrix_from_columns(rows_order).pivots()]
    cols = [cols_order[j] for j in M.matrix_from_columns(cols_order).pivots()]
    d1 = M.matrix_from_rows_and_columns(rows, cols).det()
    bound = sum(2*mon9[j][3] for j in cols) + sum(labels[i][0] % 2 - 2*labels[i][1][3] for i in rows)
    assert d1 == 12473 and bound == 835
    print("Degree-ordered minor", d1, "degree bound", bound, flush=True)

if __name__=="__main__":
    certify(reduced_residual_rows(5,16,even_linear=True,primitive=False))
